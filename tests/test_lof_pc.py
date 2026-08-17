import csv
import gzip
import importlib
import json
import logging
import math
from pathlib import Path

import numpy as np
import pytest


def lof_pc_module():
    return importlib.import_module("rare_variant_enrichment.lof_pc")


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("ENSG000001.12", "ENSG000001"),
        ("ENSG000001.12.3", "ENSG000001.12"),
        ("ENSG000001.alpha", "ENSG000001.alpha"),
        ("ENSG000001.", "ENSG000001."),
        ("001.2", "001"),
    ],
)
def test_normalize_ensembl_id_strips_exactly_one_terminal_numeric_version(
    raw: str, expected: str
):
    assert lof_pc_module().normalize_ensembl_id(raw) == expected


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("ENSG00000000419.14", "ENSG00000000419"),
        ("A0JNW5_ENSG00000111647.13", "ENSG00000111647"),
        (
            "chr20:50941209:50942031:clu_63027_-:ENSG00000000419.14",
            "ENSG00000000419",
        ),
        (" ENSG1 ", "ENSG1"),
    ],
)
def test_normalize_molecular_phenotype_id(raw: str, expected: str):
    assert lof_pc_module().normalize_molecular_phenotype_id(raw, 7) == expected


@pytest.mark.parametrize("raw", ["", "A0JNW5", "ENSG1_ENSG2"])
def test_normalize_molecular_phenotype_id_rejects_invalid_ids(raw: str):
    with pytest.raises(ValueError, match="Line 7"):
        lof_pc_module().normalize_molecular_phenotype_id(raw, 7)


def test_collapse_gene_expression_rows_uses_minimum_finite_z_score():
    rows = [
        ("ENSG1", np.array([1.0, np.nan, -2.0])),
        ("ENSG1", np.array([0.5, -3.0, np.nan])),
        ("ENSG2", np.array([4.0, 5.0, 6.0])),
    ]

    collapsed = lof_pc_module()._collapse_gene_expression_rows(rows)

    assert [gene_id for gene_id, _ in collapsed] == ["ENSG1", "ENSG2"]
    np.testing.assert_allclose(
        collapsed[0][1], np.array([0.5, -3.0, -2.0]), equal_nan=True
    )
    np.testing.assert_allclose(collapsed[1][1], np.array([4.0, 5.0, 6.0]))


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_prepare_protein_coding_genes_streams_gtf_and_writes_sorted_unique_ids(
    tmp_path: Path, compressed: bool
):
    gtf_text = (
        "##description test\n"
        "chr1\tsrc\texon\t1\t2\t.\t+\t.\tgene_id \"ENSG_EXON.1\"; gene_type \"protein_coding\";\n"
        "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG2.7\"; gene_type \"protein_coding\";\n"
        "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG1.2.3\"; gene_type \"protein_coding\";\n"
        "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG2.8\"; gene_type \"protein_coding\";\n"
        "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG3.1\"; gene_type \"lncRNA\";\n"
        "chr1\ttoo\tfew\tfields\n"
    )
    gtf = tmp_path / ("genes.gtf.gz" if compressed else "genes.gtf")
    if compressed:
        with gzip.open(gtf, "wt", encoding="utf-8") as handle:
            handle.write(gtf_text)
    else:
        gtf.write_text(gtf_text)
    genes = tmp_path / "coding.tsv"
    qc = tmp_path / "coding.json"

    lof_pc_module().prepare_protein_coding_genes(gtf, genes, qc)

    assert genes.read_text().splitlines() == ["gene_id", "ENSG1.2", "ENSG2"]
    payload = json.loads(qc.read_text())
    assert payload == {
        "duplicate_normalized_gene_records": 1,
        "nine_field_gene_records": 4,
        "protein_coding_gene_records": 3,
        "total_lines": 7,
        "unique_protein_coding_gene_count": 2,
    }


def test_prepare_protein_coding_genes_rejects_gtf_without_coding_genes(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"ENSG1.1\"; gene_type \"lncRNA\";\n"
    )

    with pytest.raises(ValueError, match="No protein-coding genes found in GTF"):
        lof_pc_module().prepare_protein_coding_genes(
            gtf, tmp_path / "genes.tsv", tmp_path / "qc.json"
        )


@pytest.mark.parametrize(
    ("available", "expected"),
    [
        (3, [0, 1, 2]),
        (12, [*range(0, 11), 11]),
        (105, [*range(0, 11), *range(20, 101, 10), 104]),
        (
            603,
            [
                *range(0, 11),
                *range(20, 101, 10),
                *range(150, 501, 50),
                600,
                602,
            ],
        ),
    ],
)
def test_adaptive_pc_grid_spans_all_required_intervals(
    available: int, expected: list[int]
):
    assert lof_pc_module().build_pc_grid([], available) == expected


@pytest.mark.parametrize(
    "requested",
    [[0, 0], [1, 0], [-1], [0, 3], [True]],
    ids=["duplicate", "decreasing", "negative", "unavailable", "boolean"],
)
def test_explicit_pc_grid_is_strictly_increasing_unique_nonnegative_and_available(
    requested: list[int],
):
    with pytest.raises(ValueError, match="PC counts"):
        lof_pc_module().build_pc_grid(requested, 3)


def test_build_pc_chunks_partitions_explicit_grid_with_short_final_chunk():
    assert lof_pc_module().build_pc_chunks([0, 1, 10, 20, 30], 31, 2) == [
        [0, 1],
        [10, 20],
        [30],
    ]


def test_build_pc_chunks_uses_adaptive_grid_when_requested_grid_is_empty():
    assert lof_pc_module().build_pc_chunks([], 25, 3) == [
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8],
        [9, 10, 20],
        [24],
    ]


def test_explicit_pc_grid_rejects_the_full_available_pc_count():
    with pytest.raises(ValueError, match="PC counts"):
        lof_pc_module().build_pc_grid([0, 3], 3)


@pytest.mark.parametrize("available", [0, 1])
def test_adaptive_pc_grid_keeps_intercept_only_model_when_no_pc_is_allowed(available: int):
    assert lof_pc_module().build_pc_grid([], available) == [0]


def test_build_pc_chunks_rejects_nonpositive_chunk_size():
    with pytest.raises(ValueError, match="positive"):
        lof_pc_module().build_pc_chunks([0], 1, 0)


def test_read_principal_component_header_does_not_parse_data_rows(tmp_path: Path):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\tPC2\nS1\tnot-a-number\tinf\n")

    assert lof_pc_module().read_principal_component_header(pcs) == 2


def test_read_principal_components_preserves_string_ids_and_all_finite_values(
    tmp_path: Path,
):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\tPC2\n001\t0.5\t-2\nA-2\t1\t3.25\n")

    matrix = lof_pc_module().read_principal_components(pcs)

    assert matrix.sample_ids == ("001", "A-2")
    assert matrix.available_pc_count == 2
    np.testing.assert_allclose(matrix.values, [[0.5, -2.0], [1.0, 3.25]])


def test_read_principal_components_rejects_nonfinite_trailing_unselected_pc(
    tmp_path: Path,
):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\tPC2\nS1\t0\tinf\nS2\t1\t2\n")

    with pytest.raises(ValueError, match="finite"):
        lof_pc_module().read_principal_components(pcs)


def test_read_principal_components_rejects_header_only_file(tmp_path: Path):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\n")

    with pytest.raises(ValueError, match="at least one sample"):
        lof_pc_module().read_principal_components(pcs)


def test_read_additional_covariates_accepts_sample_id_as_final_column(
    tmp_path: Path,
):
    path = tmp_path / "genetic-pcs.tsv"
    path.write_text(
        "GENETICPC1\tGENETICPC2\tsample_id\n"
        "0.1\t-0.2\t001\n"
        "0.3\t0.4\t1000291\n"
    )

    matrix = lof_pc_module().read_additional_covariates(path)

    assert matrix.sample_ids == ("001", "1000291")
    assert matrix.names == ("GENETICPC1", "GENETICPC2")
    assert matrix.sample_count == 2
    assert matrix.covariate_count == 2
    np.testing.assert_allclose(matrix.values, [[0.1, -0.2], [0.3, 0.4]])


def test_read_additional_covariates_accepts_id_column_in_any_position(
    tmp_path: Path,
):
    path = tmp_path / "covariates.tsv"
    path.write_text("ID\tGENETICPC1\n001\t0.1\n")

    matrix = lof_pc_module().read_additional_covariates(path)

    assert matrix.sample_ids == ("001",)
    assert matrix.names == ("GENETICPC1",)


@pytest.mark.parametrize(
    "text",
    [
        "GENETICPC1\tsample_id\n0.1\tS1\n0.2\tS1\n",
        "GENETICPC1\tother\n0.1\tS1\n",
        "sample_id\tGENETICPC1\nS1\tNA\n",
        "sample_id\tGENETICPC1\nS1\tinf\n",
        "sample_id\tGENETICPC1\n\t0.1\n",
        "sample_id\tsample_id\nS1\t0.1\n",
        "sample_id\nS1\n",
    ],
)
def test_read_additional_covariates_rejects_invalid_schema_or_values(
    tmp_path: Path, text: str
):
    path = tmp_path / "invalid.tsv"
    path.write_text(text)

    with pytest.raises(ValueError):
        lof_pc_module().read_additional_covariates(path)


def test_lof_pc_sample_id_readers_preserve_numeric_looking_ids_as_strings(
    tmp_path: Path,
):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\n001\t0.1\n1000291\t0.2\n")

    assert lof_pc_module().read_principal_components(pcs).sample_ids == (
        "001",
        "1000291",
    )


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("sample\tPC1\nS1\t0\n", "header"),
        ("ID\tPC1\tPC3\nS1\t0\t1\n", "consecutive"),
        ("ID\tPC1\nS1\t0\nS1\t1\n", "Duplicate PC sample ID"),
        ("ID\tPC1\n\t0\n", "empty sample ID"),
    ],
)
def test_read_principal_components_rejects_schema_and_id_errors(
    tmp_path: Path, text: str, message: str
):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text(text)

    with pytest.raises(ValueError, match=message):
        lof_pc_module().read_principal_components(pcs)


def test_read_lof_carriers_collapses_normalized_pairs_and_class_tokens(tmp_path: Path):
    carriers = tmp_path / "lof.tsv"
    carriers.write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "001\tENSG1.2\tG1\tYES\t1\tv1\t hc , LC \n"
        "001\tENSG1.3\tG1\ttrue\t1\tv2\tlc\n"
        "002\tENSG2\tG2\t1\t1\tv3\tunknown\n"
        "003\tENSG3\tG3\tno\t1\tv4\tHC\n"
        "004\tENSG4\tG4\t0\t0\t\t\n"
    )

    parsed = lof_pc_module().read_lof_carriers(carriers)

    assert parsed.pairs_by_definition == {
        "any_lof": {("001", "ENSG1"), ("002", "ENSG2")},
        "HC": {("001", "ENSG1")},
        "HC_or_LC": {("001", "ENSG1")},
    }
    assert parsed.qc == {
        "input_row_count": 5,
        "truthy_row_count": 3,
        "unique_any_lof_pair_count": 2,
        "unique_hc_or_lc_pair_count": 1,
        "unique_hc_pair_count": 1,
    }


def test_read_lof_carriers_rejects_unknown_boolean(tmp_path: Path):
    carriers = tmp_path / "lof.tsv"
    carriers.write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "S1\tENSG1\tG1\tmaybe\t1\tv1\tHC\n"
    )

    with pytest.raises(ValueError, match="true/false/1/0/yes/no"):
        lof_pc_module().read_lof_carriers(carriers)


def test_residualize_expression_uses_intercept_selected_pcs_and_population_sd():
    fit = lof_pc_module().residualize_expression(
        np.array([3.0, 1.0, 4.0, 6.0, 8.0]),
        np.array([[0.0], [0.0], [1.0], [1.0], [2.0]]),
        1,
    )

    assert fit.exclusion_reason is None
    assert fit.usable_sample_count == 5
    assert fit.rank == 2
    assert fit.residual_mean == pytest.approx(0.0, abs=1e-14)
    assert fit.residual_sd == pytest.approx(math.sqrt(4.0 / 5.0))
    np.testing.assert_allclose(
        fit.z_scores,
        [math.sqrt(5) / 2, -math.sqrt(5) / 2, -math.sqrt(5) / 2, math.sqrt(5) / 2, 0],
        atol=1e-12,
    )


def test_residualize_expression_keeps_missing_observations_nan():
    fit = lof_pc_module().residualize_expression(
        np.array([0.0, np.nan, 0.0, 2.0, 2.0]),
        np.zeros((5, 0)),
        0,
    )

    assert fit.exclusion_reason is None
    assert fit.usable_sample_count == 4
    assert np.isnan(fit.z_scores[1])
    np.testing.assert_allclose(fit.z_scores[[0, 2, 3, 4]], [-1, -1, 1, 1])


def test_complete_data_projection_matches_legacy_residuals():
    expression = np.array(
        [
            [3.0, 2.0],
            [1.0, 4.0],
            [4.0, 1.0],
            [6.0, 3.0],
            [8.0, 7.0],
        ]
    )
    pcs = np.array(
        [[-2.0, 1.0], [-1.0, -1.0], [0.0, 0.0], [1.0, -1.0], [2.0, 1.0]]
    )
    module = lof_pc_module()
    state = module.prepare_complete_data_projection(expression, pcs, [0, 1])
    prediction = np.zeros_like(expression)
    previous_pc_count = 0

    for pc_count in [0, 1]:
        prediction = state.advance_prediction(
            previous_pc_count, pc_count, prediction
        )
        actual = state.z_scores(prediction)
        expected = np.column_stack(
            [
                module.residualize_expression(expression[:, column], pcs, pc_count).z_scores
                for column in range(expression.shape[1])
            ]
        )
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-10)
        previous_pc_count = pc_count


def test_complete_data_projection_matches_legacy_for_nonorthogonal_pc_prefixes():
    expression = np.array(
        [
            [1.0, 6.0],
            [2.0, 4.0],
            [4.0, 7.0],
            [3.0, 3.0],
            [7.0, 5.0],
            [8.0, 9.0],
        ]
    )
    pcs = np.array(
        [[0.0, 0.0], [1.0, 0.0], [2.0, 1.0], [3.0, 1.0], [4.0, 2.0], [5.0, 3.0]]
    )
    module = lof_pc_module()
    state = module.prepare_complete_data_projection(expression, pcs, [0, 1, 2])
    prediction = np.zeros_like(expression)
    previous_pc_count = 0

    for pc_count in [0, 1, 2]:
        prediction = state.advance_prediction(
            previous_pc_count, pc_count, prediction
        )
        expected = np.column_stack(
            [
                module.residualize_expression(expression[:, column], pcs, pc_count).z_scores
                for column in range(expression.shape[1])
            ]
        )
        np.testing.assert_allclose(
            state.z_scores(prediction), expected, atol=1e-10, rtol=1e-10
        )
        previous_pc_count = pc_count


def test_complete_data_projection_marks_zero_variance_expression_as_nan():
    state = lof_pc_module().prepare_complete_data_projection(
        np.array([[5.0], [5.0], [5.0]]),
        np.array([[-1.0], [0.0], [1.0]]),
        [0, 1],
    )

    prediction = state.advance_prediction(0, 1, np.zeros((3, 1)))

    assert np.isnan(state.z_scores(prediction)).all()


def test_complete_data_projection_rejects_missing_expression_values():
    with pytest.raises(ValueError, match="complete finite expression"):
        lof_pc_module().prepare_complete_data_projection(
            np.array([[1.0], [np.nan], [3.0]]),
            np.array([[-1.0], [0.0], [1.0]]),
            [0, 1],
        )


@pytest.mark.parametrize(
    "thresholds", [[], [0.0], [1.0], [-2.0, -2.0], [-2.0, float("nan")]]
)
def test_negative_z_thresholds_must_be_finite_unique_and_negative(
    thresholds: list[float],
):
    with pytest.raises(ValueError, match="negative z-score thresholds"):
        lof_pc_module().validate_negative_z_thresholds(thresholds)


def _write_analysis_fixture(tmp_path: Path) -> dict[str, Path]:
    phenotype = tmp_path / "phenotypes.bed"
    phenotype.write_text(
        "#chr\tstart\tend\tfeature_id\tS1\tS2\tS3\tS4\tS5\tS6\n"
        "chr1\t0\t1\tENSG1.4\t0\t1\t2\t3\t4\t5\n"
        "chr1\t1\t2\tENSG2.8\t5\t4\t3\t2\t1\t0\n"
        "chr1\t2\t3\tENSG3.1\t0\t1\t2\t3\t4\t5\n"
    )
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text(
        "ID\tPC1\tPC2\n"
        "S1\t-1\t0\nS2\t1\t1\nS3\t-1\t2\n"
        "S4\t1\t3\nS5\t-1\t4\nS6\t1\t5\nPC_ONLY\t0\t6\n"
    )
    genes = tmp_path / "genes.tsv"
    genes.write_text("gene_id\nENSG1\nENSG2\n")
    carriers = tmp_path / "lof.tsv"
    carriers.write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "S1\tENSG1.1\tG1\tyes\t1\tv1\tHC\n"
        "S3\tENSG1.2\tG1\ttrue\t1\tv2\tLC\n"
        "S4\tENSG1\tG1\t1\t1\tv3\tunknown\n"
        "S2\tENSG2.1\tG2\ttrue\t1\tv4\tHC\n"
        "S5\tENSG2\tG2\ttrue\t1\tv5\tlc\n"
        "SX\tENSG1\tG1\ttrue\t1\tv6\tHC\n"
        "S1\tENSGX\tGX\ttrue\t1\tv7\tHC\n"
    )
    return {
        "phenotype": phenotype,
        "pcs": pcs,
        "genes": genes,
        "carriers": carriers,
    }


def _run_analysis(
    tmp_path: Path,
    inputs: dict[str, Path],
    *,
    thresholds: list[float],
    pc_counts: list[int],
) -> dict[str, Path]:
    outputs = {
        "results": tmp_path / "results.tsv",
        "summary": tmp_path / "summary.json",
        "gene_qc": tmp_path / "gene-pc-qc.tsv.gz",
        "analysis_qc": tmp_path / "analysis-qc.json",
    }
    lof_pc_module().calculate_lof_pc_enrichment(
        inputs["phenotype"],
        inputs["carriers"],
        inputs["pcs"],
        inputs["genes"],
        thresholds,
        pc_counts,
        outputs["results"],
        outputs["summary"],
        outputs["gene_qc"],
        outputs["analysis_qc"],
    )
    return outputs


def test_lof_pc_enrichment_accepts_susie_ids_and_collapses_duplicate_features(
    tmp_path: Path,
):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["phenotype"].write_text(
        "#chr\tstart\tend\tfeature_id\tS1\tS2\tS3\tS4\tS5\tS6\n"
        "chr1\t0\t1\tA0JNW5_ENSG1.7\t0\t1\t2\t3\t4\t5\n"
        "chr1\t0\t1\tchr1:100:101:clu_1_+:ENSG1.9\t0\t-1\t-2\t-3\t-4\t-5\n"
        "chr1\t1\t2\tchr1:200:201:clu_2_-:ENSG2.3\t5\t4\t3\t2\t1\tNA\n"
    )

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0])

    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    assert analysis_qc["bed_feature_count"] == 3
    assert analysis_qc["bed_gene_count"] == 2
    assert analysis_qc["duplicate_feature_count"] == 1
    assert analysis_qc["protein_coding_bed_feature_count"] == 3
    assert analysis_qc["protein_coding_bed_gene_count"] == 2
    assert analysis_qc["protein_coding_duplicate_feature_count"] == 1
    assert analysis_qc["per_pc"]["0"]["total_observations"] == 11

    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        gene_qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {(row["gene_id"], row["pc_count"]) for row in gene_qc_rows} == {
        ("ENSG1", "0"),
        ("ENSG2", "0"),
    }


def test_calculate_lof_pc_enrichment_preserves_supplied_adaptive_grid_mode(
    tmp_path: Path,
):
    """A chunked adaptive WDL shard must not be relabeled as an explicit grid."""
    inputs = _write_analysis_fixture(tmp_path)
    outputs = _run_analysis(
        tmp_path, inputs, thresholds=[-0.8], pc_counts=[0, 1]
    )

    lof_pc_module().calculate_lof_pc_enrichment(
        inputs["phenotype"],
        inputs["carriers"],
        inputs["pcs"],
        inputs["genes"],
        [-0.8],
        [0, 1],
        outputs["results"],
        outputs["summary"],
        outputs["gene_qc"],
        outputs["analysis_qc"],
        pc_grid_mode="adaptive",
    )

    assert json.loads(outputs["summary"].read_text())["pc_grid_mode"] == "adaptive"


def _write_merge_shard(
    tmp_path: Path,
    name: str,
    *,
    pc_counts: list[int],
    p_values: list[float],
) -> dict[str, Path]:
    carrier_definitions = list(lof_pc_module().CARRIER_DEFINITIONS)
    assert len(p_values) == len(carrier_definitions)
    shard_dir = tmp_path / name
    shard_dir.mkdir()
    outputs = {
        "results": shard_dir / "results.tsv",
        "summary": shard_dir / "summary.json",
        "gene_qc": shard_dir / "gene-pc-qc.tsv.gz",
        "analysis_qc": shard_dir / "analysis-qc.json",
    }
    result_header = lof_pc_module().RESULT_HEADER
    with outputs["results"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=result_header, delimiter="\t")
        writer.writeheader()
        for pc_count in pc_counts:
            for index, p_value in enumerate(p_values):
                writer.writerow(
                    {
                        "pc_count": pc_count,
                        "z_threshold": -2.0,
                        "carrier_definition": carrier_definitions[index],
                        "eligible_gene_count": 1,
                        "total_observations": 4,
                        "outlier_observations": 1,
                        "carrier_observations": 1,
                        "n11": 1,
                        "n10": 0,
                        "n01": 0,
                        "n00": 3,
                        "outlier_carrier_rate": "1.0",
                        "nonoutlier_carrier_rate": "0.0",
                        "carrier_rate_ratio": "NA",
                        "odds_ratio": "NA",
                        "odds_ratio_corrected_0_5": "21.0",
                        "fisher_p_value": p_value,
                        "fisher_fdr_bh": "1.0",
                    }
                )
    with gzip.open(outputs["gene_qc"], "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(lof_pc_module().GENE_PC_QC_HEADER)
        for pc_count in pc_counts:
            writer.writerow([f"GENE_{pc_count}", pc_count, 4, 1, 0.0, 1.0, "included", ""])
    summary = {
        "available_pc_count": 2,
        "carrier_definitions": carrier_definitions,
        "emitted_result_rows": len(pc_counts) * len(p_values),
        "fdr_scope": "global_across_all_emitted_rows",
        "negative_z_thresholds": [-2.0],
        "observation_unit": "eligible sample-gene residual",
        "pc_grid_mode": "explicit",
        "provenance": {
            "input_files": {
                "lof_carriers": "lof.tsv",
                "phenotype_bed": "input.bed",
                "principal_components": "pcs.tsv",
                "protein_coding_genes": "genes.tsv",
            },
            "software_versions": {
                "numpy": "2.0.0",
                "python": "3.12.0",
                "rare_variant_enrichment": "0.3.0",
            },
        },
        "residualization": {
            "design": "intercept plus first k principal components",
            "outlier_rule": "residual_z <= z_threshold",
            "residual_standard_deviation_ddof": 0,
        },
        "selected_pc_counts": pc_counts,
        "statistical_limitation": "screening",
    }
    outputs["summary"].write_text(json.dumps(summary))
    analysis_qc = {
        "bed_gene_count": 1,
        "bed_sample_count": 4,
        "pc_sample_count": 4,
        "protein_coding_bed_gene_count": 1,
        "protein_coding_gene_count": 1,
        "shared_bed_pc_sample_count": 4,
        "pre_join_carrier_pair_counts": {"any_lof": 1, "HC": 1, "HC_or_LC": 1},
        "lof_carrier_table": {"input_row_count": 1},
        "per_pc": {
            str(pc_count): {
                "eligible_gene_count": 1,
                "total_observations": 4,
                "carrier_observations": {"any_lof": 1, "HC": 1, "HC_or_LC": 1},
                "exclusion_counts": {
                    "insufficient_dof": 0,
                    "invalid_or_zero_residual_sd": 0,
                    "other": 0,
                    "rank_deficiency": 0,
                },
            }
            for pc_count in pc_counts
        },
    }
    outputs["analysis_qc"].write_text(json.dumps(analysis_qc))
    return outputs


def _read_merge_result_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_merge_result_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=lof_pc_module().RESULT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_gene_pc_qc_rows(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def _write_gene_pc_qc_rows(path: Path, rows: list[list[str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def _rewrite_merge_summary(path: Path, **changes: object) -> None:
    summary = json.loads(path.read_text())
    summary.update(changes)
    path.write_text(json.dumps(summary))


def _rewrite_analysis_qc(path: Path, **changes: object) -> None:
    analysis_qc = json.loads(path.read_text())
    analysis_qc.update(changes)
    path.write_text(json.dumps(analysis_qc))


def _merge_shards(tmp_path: Path, shards: list[dict[str, Path]]) -> dict[str, Path]:
    outputs = {
        "results": tmp_path / "merged-results.tsv",
        "summary": tmp_path / "merged-summary.json",
        "gene_qc": tmp_path / "merged-gene-pc-qc.tsv.gz",
        "analysis_qc": tmp_path / "merged-analysis-qc.json",
    }
    lof_pc_module().merge_lof_pc_enrichment(
        [shard["results"] for shard in shards],
        [shard["summary"] for shard in shards],
        [shard["gene_qc"] for shard in shards],
        [shard["analysis_qc"] for shard in shards],
        outputs["results"],
        outputs["summary"],
        outputs["gene_qc"],
        outputs["analysis_qc"],
    )
    return outputs


def test_merge_lof_pc_enrichment_recomputes_global_fdr_and_combines_qc(tmp_path: Path):
    shard_one = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    shard_two = _write_merge_shard(
        tmp_path, "two", pc_counts=[1], p_values=[0.02, 0.50, 0.60]
    )

    outputs = _merge_shards(tmp_path, [shard_one, shard_two])

    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert [(row["pc_count"], row["z_threshold"], row["carrier_definition"]) for row in rows] == [
        ("0", "-2.0", "any_lof"),
        ("0", "-2.0", "HC"),
        ("0", "-2.0", "HC_or_LC"),
        ("1", "-2.0", "any_lof"),
        ("1", "-2.0", "HC"),
        ("1", "-2.0", "HC_or_LC"),
    ]
    assert [row["fisher_fdr_bh"] for row in rows] == [
        "0.06",
        "0.4000000000000001",
        "0.44999999999999996",
        "0.06",
        "0.6",
        "0.6",
    ]
    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 2
    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    assert set(analysis_qc["per_pc"]) == {"0", "1"}
    assert analysis_qc["bed_gene_count"] == 1
    summary = json.loads(outputs["summary"].read_text())
    assert summary["selected_pc_counts"] == [0, 1]
    assert summary["emitted_result_rows"] == 6
    assert summary["fdr_scope"] == "global_across_all_emitted_rows"


def test_merge_lof_pc_enrichment_rejects_duplicate_pc_counts(tmp_path: Path):
    shard_one = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    shard_two = _write_merge_shard(
        tmp_path, "two", pc_counts=[0], p_values=[0.02, 0.50, 0.60]
    )

    with pytest.raises(ValueError, match="Duplicate PC count"):
        _merge_shards(tmp_path, [shard_one, shard_two])


def test_merge_lof_pc_enrichment_rejects_results_missing_a_selected_pc(tmp_path: Path):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    shard["results"].write_text("\t".join(lof_pc_module().RESULT_HEADER) + "\n")

    with pytest.raises(ValueError, match="do not match the shard summary"):
        _merge_shards(tmp_path, [shard])


def test_merge_lof_pc_enrichment_rejects_fisher_p_value_outside_probability_range(
    tmp_path: Path,
):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    rows = _read_merge_result_rows(shard["results"])
    rows[0]["fisher_p_value"] = "1.01"
    _write_merge_result_rows(shard["results"], rows)

    with pytest.raises(ValueError, match="fisher_p_value"):
        _merge_shards(tmp_path, [shard])


@pytest.mark.parametrize("row_index", [0, -1], ids=["duplicate", "missing"])
def test_merge_lof_pc_enrichment_rejects_duplicate_or_missing_result_combinations(
    tmp_path: Path, row_index: int
):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    rows = _read_merge_result_rows(shard["results"])
    if row_index == 0:
        rows[1]["carrier_definition"] = rows[0]["carrier_definition"]
    else:
        rows.pop()
    _write_merge_result_rows(shard["results"], rows)

    with pytest.raises(ValueError, match="result combinations"):
        _merge_shards(tmp_path, [shard])


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("duplicate", "duplicate gene-PC QC row"),
        ("missing", "expected gene-PC QC rows"),
        ("wrong_pc", "unselected PC count"),
        ("excluded", "included row count"),
        ("wrong_usable_count", "total_observations"),
    ],
)
def test_merge_lof_pc_enrichment_rejects_incompatible_gene_pc_qc(
    tmp_path: Path, mutation: str, message: str
):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    rows = _read_gene_pc_qc_rows(shard["gene_qc"])
    if mutation == "duplicate":
        rows.append(rows[1])
    elif mutation == "missing":
        rows.pop()
    elif mutation == "wrong_pc":
        rows[1][1] = "1"
    elif mutation == "excluded":
        rows[1][6:] = ["excluded", "rank_deficiency"]
    elif mutation == "wrong_usable_count":
        rows[1][2] = "3"
    else:
        raise AssertionError(f"Unhandled mutation: {mutation}")
    _write_gene_pc_qc_rows(shard["gene_qc"], rows)

    with pytest.raises(ValueError, match=message):
        _merge_shards(tmp_path, [shard])


def test_merge_lof_pc_enrichment_requires_identical_analysis_qc_metadata(
    tmp_path: Path,
):
    shard_one = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    shard_two = _write_merge_shard(
        tmp_path, "two", pc_counts=[1], p_values=[0.02, 0.50, 0.60]
    )
    _rewrite_analysis_qc(shard_two["analysis_qc"], bed_sample_count=5)

    with pytest.raises(ValueError, match="top-level metadata"):
        _merge_shards(tmp_path, [shard_one, shard_two])


def test_merge_lof_pc_enrichment_rejects_invalid_available_pc_count(tmp_path: Path):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    _rewrite_merge_summary(shard["summary"], available_pc_count=True)

    with pytest.raises(ValueError, match="available_pc_count"):
        _merge_shards(tmp_path, [shard])


def test_merge_lof_pc_enrichment_rejects_malformed_provenance(tmp_path: Path):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    _rewrite_merge_summary(shard["summary"], provenance={"input_files": {}})

    with pytest.raises(ValueError, match="provenance"):
        _merge_shards(tmp_path, [shard])


def test_merge_lof_pc_enrichment_forces_global_fdr_scope(tmp_path: Path):
    shard = _write_merge_shard(
        tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20, 0.30]
    )
    _rewrite_merge_summary(shard["summary"], fdr_scope="shard_local")

    outputs = _merge_shards(tmp_path, [shard])

    assert json.loads(outputs["summary"].read_text())["fdr_scope"] == (
        "global_across_all_emitted_rows"
    )


def test_analysis_logs_configuration_progress_and_outputs(tmp_path: Path, caplog):
    inputs = _write_analysis_fixture(tmp_path)

    with caplog.at_level(logging.INFO, logger="rare_variant_enrichment.lof_pc"):
        outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0])

    messages = "\n".join(caplog.messages)
    assert "Starting LoF/PC enrichment" in messages
    assert "shared BED/PC samples" in messages
    assert "Completed PC count 0" in messages
    assert "Wrote LoF/PC enrichment outputs" in messages
    assert str(outputs["results"]) in messages


def test_pc_major_complete_data_analysis_uses_incremental_projection_and_logs_in_order(
    tmp_path: Path, caplog, monkeypatch: pytest.MonkeyPatch
):
    inputs = _write_analysis_fixture(tmp_path)
    module = lof_pc_module()
    original_prepare = module.prepare_complete_data_projection
    original_advance = module.CompleteDataProjection.advance_prediction
    original_residualize = module.residualize_expression
    prepared_shapes: list[tuple[tuple[int, int], tuple[int, int], tuple[int, ...]]] = []
    advance_calls: list[tuple[int, int]] = []
    residualize_calls: list[int] = []

    def prepare(expression, principal_components, requested_pc_counts):
        prepared_shapes.append(
            (expression.shape, principal_components.shape, tuple(requested_pc_counts))
        )
        return original_prepare(expression, principal_components, requested_pc_counts)

    def advance(self, previous_pc_count, pc_count, prediction):
        advance_calls.append((previous_pc_count, pc_count))
        return original_advance(self, previous_pc_count, pc_count, prediction)

    def residualize(expression, principal_components, pc_count):
        residualize_calls.append(pc_count)
        return original_residualize(expression, principal_components, pc_count)

    monkeypatch.setattr(module, "prepare_complete_data_projection", prepare)
    monkeypatch.setattr(module.CompleteDataProjection, "advance_prediction", advance)
    monkeypatch.setattr(module, "residualize_expression", residualize)

    with caplog.at_level(logging.INFO, logger="rare_variant_enrichment.lof_pc"):
        outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0, 1])

    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert {row["pc_count"] for row in rows} == {"0", "1"}
    assert prepared_shapes == [((6, 2), (6, 2), (0, 1))]
    assert advance_calls == [(0, 0), (0, 1)]
    assert residualize_calls == []
    completion_indexes = {
        pc_count: next(
            index
            for index, message in enumerate(caplog.messages)
            if f"Completed PC count {pc_count}" in message
        )
        for pc_count in (0, 1)
    }
    assert completion_indexes[0] < completion_indexes[1]


def test_vectorized_multi_pc_nonorthogonal_analysis_matches_legacy_results_and_qc(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """The complete-data optimization must retain legacy outputs for nonorthogonal PCs."""
    inputs = _write_analysis_fixture(tmp_path)
    vectorized_directory = tmp_path / "vectorized"
    vectorized_directory.mkdir()
    vectorized_outputs = _run_analysis(
        vectorized_directory, inputs, thresholds=[-0.8], pc_counts=[0, 1]
    )

    module = lof_pc_module()
    original_matrix_rank = module.np.linalg.matrix_rank
    rank_calls = 0

    def force_initial_rank_deficiency(*arguments, **kwargs):
        nonlocal rank_calls
        rank_calls += 1
        if rank_calls == 1:
            return 0
        return original_matrix_rank(*arguments, **kwargs)

    monkeypatch.setattr(module.np.linalg, "matrix_rank", force_initial_rank_deficiency)
    legacy_directory = tmp_path / "legacy"
    legacy_directory.mkdir()
    legacy_outputs = _run_analysis(
        legacy_directory, inputs, thresholds=[-0.8], pc_counts=[0, 1]
    )

    assert _read_merge_result_rows(vectorized_outputs["results"]) == _read_merge_result_rows(
        legacy_outputs["results"]
    )
    vectorized_qc_rows = _read_gene_pc_qc_rows(vectorized_outputs["gene_qc"])
    legacy_qc_rows = _read_gene_pc_qc_rows(legacy_outputs["gene_qc"])
    assert [row[:4] + row[6:] for row in vectorized_qc_rows] == [
        row[:4] + row[6:] for row in legacy_qc_rows
    ]
    for vectorized_row, legacy_row in zip(
        vectorized_qc_rows[1:], legacy_qc_rows[1:], strict=True
    ):
        for column in (4, 5):
            if vectorized_row[column] == "NA":
                assert legacy_row[column] == "NA"
            else:
                assert float(vectorized_row[column]) == pytest.approx(
                    float(legacy_row[column]), abs=1e-12
                )
    assert json.loads(vectorized_outputs["analysis_qc"].read_text()) == json.loads(
        legacy_outputs["analysis_qc"].read_text()
    )


def test_missing_expression_fallback_retains_legacy_results_and_qc_exclusions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["phenotype"].write_text(
        "#chr\tstart\tend\tfeature_id\tS1\tS2\tS3\tS4\tS5\tS6\n"
        "chr1\t0\t1\tENSG1.4\t0\t1\tNA\tNA\tNA\tNA\n"
        "chr1\t1\t2\tENSG2.8\t5\t4\t3\t2\t1\t0\n"
    )
    module = lof_pc_module()

    def fail_if_complete_data_path_is_used(*args, **kwargs):
        raise AssertionError("missing expression must select the legacy route")

    monkeypatch.setattr(
        module, "prepare_complete_data_projection", fail_if_complete_data_path_is_used
    )

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0, 1])

    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert {row["pc_count"] for row in rows} == {"0", "1"}
    assert {int(row["eligible_gene_count"]) for row in rows} == {1}
    assert {int(row["total_observations"]) for row in rows} == {6}

    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    for pc_count in ("0", "1"):
        assert analysis_qc["per_pc"][pc_count]["exclusion_counts"] == {
            "insufficient_dof": 1,
            "invalid_or_zero_residual_sd": 0,
            "other": 0,
            "rank_deficiency": 0,
        }

    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        gene_qc = list(csv.DictReader(handle, delimiter="\t"))
    by_gene_pc = {(row["gene_id"], row["pc_count"]): row for row in gene_qc}
    for pc_count in ("0", "1"):
        assert by_gene_pc[("ENSG1", pc_count)]["status"] == "excluded"
        assert (
            by_gene_pc[("ENSG1", pc_count)]["exclusion_reason"]
            == "insufficient_dof"
        )
        assert by_gene_pc[("ENSG2", pc_count)]["status"] == "included"


def test_rank_deficient_complete_data_preserves_legacy_qc_exclusions(tmp_path: Path):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["pcs"].write_text(
        "ID\tPC1\tPC2\n"
        "S1\t0\t0\nS2\t0\t1\nS3\t0\t2\nS4\t0\t3\nS5\t0\t4\nS6\t0\t5\nPC_ONLY\t0\t6\n"
    )

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0, 1])

    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert {row["pc_count"] for row in rows} == {"0", "1"}
    assert {int(row["eligible_gene_count"]) for row in rows if row["pc_count"] == "0"} == {
        2
    }
    assert {int(row["eligible_gene_count"]) for row in rows if row["pc_count"] == "1"} == {
        0
    }

    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    assert analysis_qc["per_pc"]["0"]["exclusion_counts"] == {
        "insufficient_dof": 0,
        "invalid_or_zero_residual_sd": 0,
        "other": 0,
        "rank_deficiency": 0,
    }
    assert analysis_qc["per_pc"]["1"]["exclusion_counts"] == {
        "insufficient_dof": 0,
        "invalid_or_zero_residual_sd": 0,
        "other": 0,
        "rank_deficiency": 2,
    }

    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        gene_qc = list(csv.DictReader(handle, delimiter="\t"))
    by_gene_pc = {(row["gene_id"], row["pc_count"]): row for row in gene_qc}
    for gene_id in ("ENSG1", "ENSG2"):
        assert by_gene_pc[(gene_id, "0")]["status"] == "included"
        assert by_gene_pc[(gene_id, "1")]["status"] == "excluded"
        assert by_gene_pc[(gene_id, "1")]["rank"] == "1"
        assert by_gene_pc[(gene_id, "1")]["exclusion_reason"] == "rank_deficiency"


def test_lof_pc_enrichment_emits_hand_checked_pooled_fisher_cells_and_global_fdr(
    tmp_path: Path,
):
    inputs = _write_analysis_fixture(tmp_path)

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0])

    with outputs["results"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    by_definition = {row["carrier_definition"]: row for row in rows}
    assert {
        definition: tuple(
            int(by_definition[definition][cell]) for cell in ("n11", "n10", "n01", "n00")
        )
        for definition in by_definition
    } == {
        "any_lof": (2, 3, 2, 5),
        "HC": (1, 1, 3, 7),
        "HC_or_LC": (2, 2, 2, 6),
    }
    assert {int(row["eligible_gene_count"]) for row in rows} == {2}
    assert {int(row["total_observations"]) for row in rows} == {12}
    assert {int(row["outlier_observations"]) for row in rows} == {4}
    assert all(0.0 <= float(row["fisher_fdr_bh"]) <= 1.0 for row in rows)
    assert float(by_definition["any_lof"]["outlier_carrier_rate"]) == pytest.approx(0.5)
    assert float(by_definition["any_lof"]["nonoutlier_carrier_rate"]) == pytest.approx(3 / 8)

    summary = json.loads(outputs["summary"].read_text())
    assert summary["selected_pc_counts"] == [0]
    assert summary["pc_grid_mode"] == "explicit"
    assert summary["available_pc_count"] == 2
    assert summary["negative_z_thresholds"] == [-0.8]
    assert summary["fdr_scope"] == "global_across_all_emitted_rows"
    assert set(summary["provenance"]) >= {"input_files", "software_versions"}

    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    assert analysis_qc["pre_join_carrier_pair_counts"] == {
        "HC": 4,
        "HC_or_LC": 6,
        "any_lof": 7,
    }
    assert analysis_qc["per_pc"]["0"]["carrier_observations"] == {
        "HC": 2,
        "HC_or_LC": 4,
        "any_lof": 5,
    }
    for row in rows:
        assert int(row["carrier_observations"]) == int(row["n11"]) + int(row["n10"])
        assert int(row["carrier_observations"]) == analysis_qc["per_pc"]["0"][
            "carrier_observations"
        ][row["carrier_definition"]]

    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        gene_qc = list(csv.DictReader(handle, delimiter="\t"))
    assert {(row["gene_id"], row["pc_count"], row["status"]) for row in gene_qc} == {
        ("ENSG1", "0", "included"),
        ("ENSG2", "0", "included"),
    }


def test_negative_threshold_is_inclusive(tmp_path: Path):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["phenotype"].write_text(
        "#chr\tstart\tend\tfeature_id\tS1\tS2\tS3\tS4\n"
        "chr1\t0\t1\tENSG1\t0\t0\t2\t2\n"
    )
    inputs["pcs"].write_text("ID\tPC1\nS1\t0\nS2\t1\nS3\t2\nS4\t3\n")
    inputs["genes"].write_text("gene_id\nENSG1\n")
    inputs["carriers"].write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "S1\tENSG1\tG1\ttrue\t1\tv1\tHC\n"
    )

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-1.0], pc_counts=[0])

    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert {int(row["outlier_observations"]) for row in rows} == {2}
    hc = next(row for row in rows if row["carrier_definition"] == "HC")
    assert (int(hc["n11"]), int(hc["n10"])) == (1, 0)


def test_analysis_qc_reports_pc_specific_exclusion_reasons_and_carriers(tmp_path: Path):
    phenotype = tmp_path / "phenotypes.bed"
    phenotype.write_text(
        "#chr\tstart\tend\tfeature_id\tS1\tS2\tS3\tS4\tS5\n"
        "chr1\t0\t1\tENSG10\t1\t2\t4\t8\t16\n"
        "chr1\t1\t2\tENSG11\t1\t2\tNA\tNA\tNA\n"
        "chr1\t2\t3\tENSG12\t1\t2\t4\tNA\tNA\n"
        "chr1\t3\t4\tENSG13\t5\t5\t5\t5\t5\n"
    )
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\tPC2\nS1\t0\t0\nS2\t0\t1\nS3\t0\t2\nS4\t1\t3\nS5\t2\t4\n")
    genes = tmp_path / "genes.tsv"
    genes.write_text("gene_id\nENSG10\nENSG11\nENSG12\nENSG13\n")
    carriers = tmp_path / "lof.tsv"
    carriers.write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "S1\tENSG10\tGOOD\ttrue\t1\tv1\tHC\n"
        "S1\tENSG12\tRANK\ttrue\t1\tv2\tLC\n"
    )
    inputs = {"phenotype": phenotype, "pcs": pcs, "genes": genes, "carriers": carriers}

    outputs = _run_analysis(tmp_path, inputs, thresholds=[-1.0], pc_counts=[0, 1])

    analysis_qc = json.loads(outputs["analysis_qc"].read_text())
    assert analysis_qc["per_pc"]["0"] == {
        "carrier_observations": {"HC": 1, "HC_or_LC": 2, "any_lof": 2},
        "eligible_gene_count": 2,
        "exclusion_counts": {
            "insufficient_dof": 1,
            "invalid_or_zero_residual_sd": 1,
            "other": 0,
            "rank_deficiency": 0,
        },
        "total_observations": 8,
    }
    assert analysis_qc["per_pc"]["1"] == {
        "carrier_observations": {"HC": 1, "HC_or_LC": 1, "any_lof": 1},
        "eligible_gene_count": 1,
        "exclusion_counts": {
            "insufficient_dof": 1,
            "invalid_or_zero_residual_sd": 1,
            "other": 0,
            "rank_deficiency": 1,
        },
        "total_observations": 5,
    }

    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    reasons = {(row["gene_id"], row["pc_count"]): row["exclusion_reason"] for row in rows}
    assert reasons[("ENSG11", "0")] == "insufficient_dof"
    assert reasons[("ENSG11", "1")] == "insufficient_dof"
    assert reasons[("ENSG12", "1")] == "rank_deficiency"
    assert reasons[("ENSG13", "0")] == "invalid_or_zero_residual_sd"
    assert reasons[("ENSG13", "1")] == "invalid_or_zero_residual_sd"


def test_analysis_rejects_trailing_infinite_pc_even_when_zero_pcs_selected(tmp_path: Path):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["pcs"].write_text("ID\tPC1\tPC2\nS1\t0\tinf\nS2\t1\t2\n")

    with pytest.raises(ValueError, match="finite"):
        _run_analysis(tmp_path, inputs, thresholds=[-1.0], pc_counts=[0])
