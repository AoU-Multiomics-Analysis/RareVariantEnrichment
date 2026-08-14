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
        (3, [0, 1, 2, 3]),
        (12, [*range(0, 11), 12]),
        (105, [*range(0, 11), *range(20, 101, 10), 105]),
        (
            603,
            [
                *range(0, 11),
                *range(20, 101, 10),
                *range(150, 501, 50),
                600,
                603,
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
    [[0, 0], [1, 0], [-1], [0, 4], [True]],
    ids=["duplicate", "decreasing", "negative", "unavailable", "boolean"],
)
def test_explicit_pc_grid_is_strictly_increasing_unique_nonnegative_and_available(
    requested: list[int],
):
    with pytest.raises(ValueError, match="PC counts"):
        lof_pc_module().build_pc_grid(requested, 3)


def test_build_pc_chunks_partitions_explicit_grid_with_short_final_chunk():
    assert lof_pc_module().build_pc_chunks([0, 1, 10, 20, 30], 30, 2) == [
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
        [25],
    ]


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
        "ID\tPC1\n"
        "S1\t0\nS2\t0\nS3\t0\nS4\t0\nS5\t0\nS6\t0\nPC_ONLY\t0\n"
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
        "chr1\t0\t1\tGOOD\t1\t2\t4\t8\t16\n"
        "chr1\t1\t2\tSPARSE\t1\t2\tNA\tNA\tNA\n"
        "chr1\t2\t3\tRANK\t1\t2\t4\tNA\tNA\n"
        "chr1\t3\t4\tZERO\t5\t5\t5\t5\t5\n"
    )
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\nS1\t0\nS2\t0\nS3\t0\nS4\t1\nS5\t2\n")
    genes = tmp_path / "genes.tsv"
    genes.write_text("gene_id\nGOOD\nSPARSE\nRANK\nZERO\n")
    carriers = tmp_path / "lof.tsv"
    carriers.write_text(
        "sample_id\tgene_id\tgene_symbol\thas_lof_variant\tn_lof_variants\tvariant_ids\tlof_classes\n"
        "S1\tGOOD\tGOOD\ttrue\t1\tv1\tHC\n"
        "S1\tRANK\tRANK\ttrue\t1\tv2\tLC\n"
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
    assert reasons[("SPARSE", "0")] == "insufficient_dof"
    assert reasons[("SPARSE", "1")] == "insufficient_dof"
    assert reasons[("RANK", "1")] == "rank_deficiency"
    assert reasons[("ZERO", "0")] == "invalid_or_zero_residual_sd"
    assert reasons[("ZERO", "1")] == "invalid_or_zero_residual_sd"


def test_analysis_rejects_trailing_infinite_pc_even_when_zero_pcs_selected(tmp_path: Path):
    inputs = _write_analysis_fixture(tmp_path)
    inputs["pcs"].write_text("ID\tPC1\tPC2\nS1\t0\tinf\nS2\t1\t2\n")

    with pytest.raises(ValueError, match="finite"):
        _run_analysis(tmp_path, inputs, thresholds=[-1.0], pc_counts=[0])
