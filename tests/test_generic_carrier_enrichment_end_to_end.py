import csv
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from rare_variant_enrichment.artifacts import file_artifact
from rare_variant_enrichment.carrier_definitions import (
    CARRIER_DEFINITION_HEADER,
    build_carrier_definitions,
)
from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER
from rare_variant_enrichment.lof_pc import (
    calculate_carrier_pc_enrichment,
    merge_carrier_pc_enrichment,
    prepare_protein_coding_genes,
)
from rare_variant_enrichment.pc_selection import analyze_carrier_pc_enrichment


FIXTURES = Path(__file__).parent / "fixtures"
DEFINITIONS = [
    "lof_hc",
    "lof_hc_or_lc",
    "missense",
    "missense_revel_ge_0_75",
    "splice_core",
    "splice_region",
    "splice_any",
    "missense_revel_ge_1",
]


def _audit_row(
    sample_id: str,
    gene_id: str,
    variant_id: str,
    variant_classes: str,
    *,
    consequence: str,
    loftee: str = "",
    revel: str = "",
) -> dict[str, str]:
    chrom, pos, ref, alt = variant_id.split(":")
    return {
        "sample_id": sample_id,
        "gene_id": gene_id,
        "gene_symbol": "G1" if gene_id == "ENSG1" else "G2",
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_id": variant_id,
        "variant_ac": "1",
        "variant_af": "0.001",
        "sample_alt_allele_count": "1",
        "most_severe_consequence": consequence,
        "all_consequences": consequence,
        "unknown_consequences": "",
        "loftee": loftee,
        "revel": revel,
        "gvs_max_af": "0.001",
        "variant_classes": variant_classes,
    }


def _write_materialization_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    rows = [
        _audit_row(
            "S1",
            "ENSG1",
            "chr1:100:A:C",
            "lof_hc,lof_hc_or_lc",
            consequence="stop_gained",
            loftee="HC",
        ),
        _audit_row(
            "S4",
            "ENSG2",
            "chr1:110:A:G",
            "lof_hc,lof_hc_or_lc",
            consequence="frameshift_variant",
            loftee="HC",
        ),
        _audit_row(
            "S6",
            "ENSG2",
            "chr1:115:A:T",
            "lof_hc_or_lc",
            consequence="stop_gained",
            loftee="LC",
        ),
        _audit_row(
            "S3",
            "ENSG1",
            "chr1:120:C:T",
            "missense",
            consequence="missense_variant",
            revel="0.74",
        ),
        _audit_row(
            "S2",
            "ENSG1",
            "chr1:125:C:G",
            "missense",
            consequence="missense_variant",
            revel="0.75",
        ),
        _audit_row(
            "S2",
            "ENSG1",
            "chr1:130:G:A",
            "missense",
            consequence="missense_variant",
            revel="0.90",
        ),
        _audit_row(
            "S4",
            "ENSG1",
            "chr1:135:G:T",
            "missense",
            consequence="missense_variant",
        ),
        _audit_row(
            "S5",
            "ENSG1",
            "chr1:140:T:C",
            "missense",
            consequence="missense_variant",
            revel="0.80",
        ),
        _audit_row(
            "S1",
            "ENSG2",
            "chr1:145:T:G",
            "splice_core",
            consequence="splice_donor_variant",
        ),
        _audit_row(
            "S6",
            "ENSG2",
            "chr1:150:A:T",
            "splice_region",
            consequence="splice_region_variant",
        ),
        _audit_row(
            "S2",
            "ENSG2",
            "chr1:155:A:G",
            "splice_core,splice_region",
            consequence="splice_donor_variant",
        ),
    ]
    audit = tmp_path / "variant_carrier_audit.tsv.gz"
    with gzip.open(audit, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=AUDIT_HEADER, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    extraction_qc = tmp_path / "variant_carriers.qc.json"
    extraction_qc.write_text(
        json.dumps(
            {
                "audit_artifact": file_artifact(
                    audit,
                    "variant_carrier_audit.tsv.gz",
                    AUDIT_HEADER,
                    len(rows),
                ),
                "vcf_index_provenance": "supplied",
                "transcript_index_provenance": "generated",
                "quality_or_frequency_filters_applied": False,
            }
        ),
        encoding="utf-8",
    )

    config = tmp_path / "carrier_definitions.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "lof_hc", "variant_classes": ["lof_hc"]},
                    {
                        "name": "lof_hc_or_lc",
                        "variant_classes": ["lof_hc_or_lc"],
                    },
                    {"name": "missense", "variant_classes": ["missense"]},
                    {
                        "name": "missense_revel_ge_0_75",
                        "variant_classes": ["missense"],
                        "minimum_revel": 0.75,
                    },
                    {
                        "name": "splice_core",
                        "variant_classes": ["splice_core"],
                    },
                    {
                        "name": "splice_region",
                        "variant_classes": ["splice_region"],
                    },
                    {
                        "name": "splice_any",
                        "variant_classes": ["splice_core", "splice_region"],
                    },
                    {
                        "name": "missense_revel_ge_1",
                        "variant_classes": ["missense"],
                        "minimum_revel": 1.0,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    return audit, extraction_qc, config


def _read_carriers(path: Path) -> list[tuple[str, ...]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == CARRIER_DEFINITION_HEADER
        return [tuple(row[column] for column in CARRIER_DEFINITION_HEADER) for row in reader]


def _calculate_shard(
    directory: Path,
    carrier_table: Path,
    manifest: Path,
    genes: Path,
    pc_count: int,
) -> dict[str, Path]:
    directory.mkdir()
    outputs = {
        "results": directory / "results.tsv",
        "summary": directory / "summary.json",
        "gene_qc": directory / "gene-pc-qc.tsv.gz",
        "analysis_qc": directory / "analysis-qc.json",
    }
    calculate_carrier_pc_enrichment(
        FIXTURES / "lof_pc_phenotypes.bed",
        carrier_table,
        manifest,
        FIXTURES / "principal_components.tsv",
        genes,
        [-0.8],
        [pc_count],
        outputs["results"],
        outputs["summary"],
        outputs["gene_qc"],
        outputs["analysis_qc"],
        pc_grid_mode="explicit",
    )
    return outputs


def _bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [1.0] * len(values)
    running = 1.0
    for rank_index in range(len(values) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * len(values) / rank)
        adjusted[original_index] = running
    return adjusted


def test_generic_carrier_enrichment_runs_from_audit_through_selection(tmp_path: Path):
    audit, extraction_qc, config = _write_materialization_inputs(tmp_path)
    carrier_table = tmp_path / "carrier_definitions.tsv.gz"
    manifest_path = tmp_path / "carrier_definitions.qc.json"
    build_carrier_definitions(
        audit,
        extraction_qc,
        config,
        carrier_table,
        manifest_path,
        container_image="example.invalid/enrichment@sha256:0123456789abcdef",
    )

    assert _read_carriers(carrier_table) == [
        ("S1", "ENSG1", "G1", "lof_hc", "1", "chr1:100:A:C"),
        ("S1", "ENSG1", "G1", "lof_hc_or_lc", "1", "chr1:100:A:C"),
        ("S1", "ENSG2", "G2", "splice_core", "1", "chr1:145:T:G"),
        ("S1", "ENSG2", "G2", "splice_any", "1", "chr1:145:T:G"),
        (
            "S2",
            "ENSG1",
            "G1",
            "missense",
            "2",
            "chr1:125:C:G,chr1:130:G:A",
        ),
        (
            "S2",
            "ENSG1",
            "G1",
            "missense_revel_ge_0_75",
            "2",
            "chr1:125:C:G,chr1:130:G:A",
        ),
        ("S2", "ENSG2", "G2", "splice_core", "1", "chr1:155:A:G"),
        ("S2", "ENSG2", "G2", "splice_region", "1", "chr1:155:A:G"),
        ("S2", "ENSG2", "G2", "splice_any", "1", "chr1:155:A:G"),
        ("S3", "ENSG1", "G1", "missense", "1", "chr1:120:C:T"),
        ("S4", "ENSG1", "G1", "missense", "1", "chr1:135:G:T"),
        ("S4", "ENSG2", "G2", "lof_hc", "1", "chr1:110:A:G"),
        ("S4", "ENSG2", "G2", "lof_hc_or_lc", "1", "chr1:110:A:G"),
        ("S5", "ENSG1", "G1", "missense", "1", "chr1:140:T:C"),
        (
            "S5",
            "ENSG1",
            "G1",
            "missense_revel_ge_0_75",
            "1",
            "chr1:140:T:C",
        ),
        ("S6", "ENSG2", "G2", "lof_hc_or_lc", "1", "chr1:115:A:T"),
        ("S6", "ENSG2", "G2", "splice_region", "1", "chr1:150:A:T"),
        ("S6", "ENSG2", "G2", "splice_any", "1", "chr1:150:A:T"),
    ]

    manifest = json.loads(manifest_path.read_text())
    assert manifest["definition_order"] == DEFINITIONS
    assert manifest["input_artifacts"]["audit"] == file_artifact(
        audit, "variant_carrier_audit.tsv.gz", AUDIT_HEADER, 11
    )
    assert manifest["output_artifact"] == file_artifact(
        carrier_table,
        "carrier_definitions.tsv.gz",
        CARRIER_DEFINITION_HEADER,
        18,
    )
    assert manifest["definition_counts"]["missense"] == {
        "matched_audit_rows": 5,
        "distinct_variants": 5,
        "distinct_sample_gene_pairs": 4,
        "output_rows": 4,
    }
    assert manifest["definition_counts"]["missense_revel_ge_0_75"] == {
        "matched_audit_rows": 3,
        "distinct_variants": 3,
        "distinct_sample_gene_pairs": 2,
        "output_rows": 2,
    }
    assert manifest["definition_counts"]["missense_revel_ge_1"] == {
        "matched_audit_rows": 0,
        "distinct_variants": 0,
        "distinct_sample_gene_pairs": 0,
        "output_rows": 0,
    }

    genes = tmp_path / "protein_coding_genes.tsv"
    genes_qc = tmp_path / "protein_coding_genes.qc.json"
    prepare_protein_coding_genes(
        FIXTURES / "gene_annotation.gtf", genes, genes_qc
    )
    shards = [
        _calculate_shard(tmp_path / "pc0", carrier_table, manifest_path, genes, 0),
        _calculate_shard(tmp_path / "pc1", carrier_table, manifest_path, genes, 1),
    ]

    merged = {
        "results": tmp_path / "carrier_pc_enrichment.tsv",
        "summary": tmp_path / "carrier_pc_enrichment.summary.json",
        "gene_qc": tmp_path / "carrier_pc_enrichment.gene_pc_qc.tsv.gz",
        "analysis_qc": tmp_path / "carrier_pc_enrichment.analysis_qc.json",
    }
    merge_carrier_pc_enrichment(
        [item["results"] for item in shards],
        [item["summary"] for item in shards],
        [item["gene_qc"] for item in shards],
        [item["analysis_qc"] for item in shards],
        merged["results"],
        merged["summary"],
        merged["gene_qc"],
        merged["analysis_qc"],
    )

    with merged["results"].open(encoding="utf-8", newline="") as handle:
        results = list(csv.DictReader(handle, delimiter="\t"))
    assert len(results) == 2 * len(DEFINITIONS)
    assert [row["carrier_definition"] for row in results] == DEFINITIONS * 2
    expected_cells = {
        "lof_hc": (1, 1, 3, 7),
        "lof_hc_or_lc": (2, 1, 2, 7),
        "missense": (1, 3, 3, 5),
        "missense_revel_ge_0_75": (1, 1, 3, 7),
        "splice_core": (0, 2, 4, 6),
        "splice_region": (1, 1, 3, 7),
        "splice_any": (1, 2, 3, 6),
        "missense_revel_ge_1": (0, 0, 4, 8),
    }
    for row in results:
        assert tuple(
            int(row[column]) for column in ("n11", "n10", "n01", "n00")
        ) == expected_cells[row["carrier_definition"]]

    expected_fdr = _bh_adjust([float(row["fisher_p_value"]) for row in results])
    assert [float(row["fisher_fdr_bh"]) for row in results] == pytest.approx(
        expected_fdr
    )
    summary = json.loads(merged["summary"].read_text())
    assert summary["carrier_definitions"] == DEFINITIONS
    assert summary["selected_pc_counts"] == [0, 1]
    assert summary["fdr_scope"] == "global_across_all_emitted_rows"
    assert summary["carrier_definition_materialization"] == {
        "definitions": manifest["definitions"],
        "manifest_artifact": {
            "logical_name": "carrier_definitions.qc.json",
            "size_bytes": manifest_path.stat().st_size,
            "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        },
        "materializer": manifest["provenance"],
    }

    analysis_qc = json.loads(merged["analysis_qc"].read_text())
    expected_carrier_counts = {
        "lof_hc": 2,
        "lof_hc_or_lc": 3,
        "missense": 4,
        "missense_revel_ge_0_75": 2,
        "splice_core": 2,
        "splice_region": 2,
        "splice_any": 3,
        "missense_revel_ge_1": 0,
    }
    assert list(analysis_qc["per_pc"]) == ["0", "1"]
    assert analysis_qc["carrier_definitions"] == DEFINITIONS
    assert analysis_qc["carrier_pair_counts"] == {
        definition: {
            "input_sample_gene_pairs": count,
            "sample_intersection_pairs": count,
            "protein_coding_gene_pairs": count,
            "complete_analysis_universe_pairs": count,
        }
        for definition, count in expected_carrier_counts.items()
    }
    for pc_count in ("0", "1"):
        assert analysis_qc["per_pc"][pc_count]["total_observations"] == 12
        assert analysis_qc["per_pc"][pc_count]["carrier_observations"] == (
            expected_carrier_counts
        )
    assert all(
        int(row["carrier_observations"]) == int(row["n11"]) + int(row["n10"])
        for row in results
    )
    assert merged["gene_qc"].read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(merged["gene_qc"], "rt", encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 4

    selection_path = tmp_path / "carrier_pc_selection.json"
    plot_path = tmp_path / "carrier_pc_enrichment.svg"
    analyze_carrier_pc_enrichment(
        merged["results"],
        selection_path,
        plot_path,
        carrier_definitions=DEFINITIONS,
        selection_z_thresholds=[-0.8],
        plateau_fraction=0.95,
    )
    selection = json.loads(selection_path.read_text())["selection"]
    assert selection["carrier_definitions"] == DEFINITIONS
    assert selection["selected_pc_count"] == 0
    assert selection["excluded_definitions"] == {
        "splice_core": "incomplete_finite_curve",
        "missense_revel_ge_1": "zero_carriers",
    }
    assert selection["estimable_carrier_definitions"] == [
        "lof_hc",
        "lof_hc_or_lc",
        "missense",
        "missense_revel_ge_0_75",
        "splice_region",
        "splice_any",
    ]
    svg = plot_path.read_text()
    assert svg.endswith("</svg>\n")
    for definition in DEFINITIONS:
        assert f'data-carrier-definition="{definition}"' in svg
    assert 'data-exclusion-reason="zero_carriers"' in svg
