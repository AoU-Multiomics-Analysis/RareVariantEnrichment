import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from rare_variant_enrichment.artifacts import file_artifact
from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER


WORKFLOW = Path("workflows/carrier_enrichment.wdl").resolve()
FIXTURES = Path(__file__).parent / "fixtures"
TEST_IMAGE = os.environ.get(
    "RARE_VARIANT_ENRICHMENT_TEST_IMAGE", "rare-variant-enrichment:test"
)


def _unavailable(message: str) -> None:
    if os.environ.get("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME") == "1":
        pytest.fail(message)
    pytest.skip(message)


def _require_wdl_runtime() -> str:
    missing = [name for name in ("miniwdl", "docker") if not shutil.which(name)]
    if missing:
        _unavailable("WDL runtime prerequisites are missing: " + ", ".join(missing))
    if subprocess.run(
        ["docker", "info"], text=True, capture_output=True, check=False
    ).returncode != 0:
        _unavailable("Docker daemon is unavailable")
    if subprocess.run(
        ["docker", "image", "inspect", TEST_IMAGE],
        text=True,
        capture_output=True,
        check=False,
    ).returncode != 0:
        _unavailable(f"Build {TEST_IMAGE} before WDL runtime tests")
    return shutil.which("miniwdl") or "miniwdl"


def _audit_row(
    sample_id: str,
    gene_id: str,
    variant_id: str,
    variant_classes: str,
    *,
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
        "most_severe_consequence": "missense_variant",
        "all_consequences": "missense_variant",
        "unknown_consequences": "",
        "loftee": "HC" if variant_classes == "lof_hc" else "",
        "revel": revel,
        "gvs_max_af": "0.001",
        "variant_classes": variant_classes,
    }


def _write_materialization_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    rows = [
        _audit_row("S1", "ENSG1", "chr1:100:A:C", "lof_hc"),
        _audit_row("S2", "ENSG2", "chr1:110:A:G", "lof_hc"),
        _audit_row("S3", "ENSG1", "chr1:120:C:T", "missense", revel="0.8"),
        _audit_row("S4", "ENSG1", "chr1:130:G:A", "missense", revel="0.7"),
        _audit_row("S5", "ENSG2", "chr1:140:T:C", "missense", revel="0.9"),
        _audit_row("S1", "ENSG2", "chr1:150:A:T", "splice_core"),
    ]
    audit = tmp_path / "variant_carrier_audit.tsv.gz"
    with gzip.open(audit, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t")
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
                "transcript_index_provenance": "supplied",
                "quality_or_frequency_filters_applied": False,
            }
        )
    )
    definitions = tmp_path / "carrier_definitions.json"
    definitions.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "lof_hc", "variant_classes": ["lof_hc"]},
                    {"name": "missense", "variant_classes": ["missense"]},
                    {
                        "name": "splice_any",
                        "variant_classes": ["splice_core", "splice_region"],
                    },
                ],
            }
        )
    )
    return audit, extraction_qc, definitions


def test_carrier_enrichment_wdl_runs_hand_checked_generic_fixture(tmp_path: Path):
    miniwdl = _require_wdl_runtime()
    audit, extraction_qc, definitions = _write_materialization_inputs(tmp_path)
    inputs = {
        "CarrierEnrichment.variant_carrier_audit_tsv_gz": str(audit.resolve()),
        "CarrierEnrichment.variant_carriers_qc_json": str(extraction_qc.resolve()),
        "CarrierEnrichment.carrier_definitions_json": str(definitions.resolve()),
        "CarrierEnrichment.phenotype_bed": str(
            (FIXTURES / "lof_pc_phenotypes.bed").resolve()
        ),
        "CarrierEnrichment.principal_components_tsv": str(
            (FIXTURES / "principal_components.tsv").resolve()
        ),
        "CarrierEnrichment.gene_annotation_gtf": str(
            (FIXTURES / "gene_annotation.gtf").resolve()
        ),
        "CarrierEnrichment.negative_z_thresholds": [-0.8],
        "CarrierEnrichment.selection_z_thresholds": [-0.8],
        "CarrierEnrichment.pc_counts": [0],
        "CarrierEnrichment.pc_counts_per_job": 1,
        "CarrierEnrichment.docker_image": TEST_IMAGE,
        "CarrierEnrichment.prepare_cpu": 1,
        "CarrierEnrichment.prepare_memory_gb": 1,
        "CarrierEnrichment.prepare_disk_gb": 1,
        "CarrierEnrichment.analysis_cpu": 1,
        "CarrierEnrichment.analysis_memory_gb": 1,
        "CarrierEnrichment.analysis_disk_gb": 1,
        "CarrierEnrichment.max_retries": 1,
        "CarrierEnrichment.pc_preemptible": 0,
    }
    inputs_path = tmp_path / "inputs.json"
    outputs_path = tmp_path / "outputs.json"
    inputs_path.write_text(json.dumps(inputs))

    result = subprocess.run(
        [
            miniwdl,
            "run",
            str(WORKFLOW),
            "-i",
            str(inputs_path),
            "-d",
            str(tmp_path / "miniwdl"),
            "-o",
            str(outputs_path),
            "--no-cache",
        ],
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    outputs = json.loads(outputs_path.read_text())["outputs"]
    assert set(outputs) == {
        f"CarrierEnrichment.{name}"
        for name in (
            "carrier_definitions_tsv_gz",
            "carrier_definitions_qc_json",
            "results_tsv",
            "summary_json",
            "gene_pc_qc_tsv_gz",
            "analysis_qc_json",
            "pc_selection_json",
            "enrichment_plot_svg",
            "pc_sweep_qc_summary_tsv",
            "pc_sweep_qc_plot_png",
            "protein_coding_genes_tsv",
            "protein_coding_genes_qc_json",
        )
    }
    with Path(outputs["CarrierEnrichment.results_tsv"]).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {
        row["carrier_definition"]: tuple(
            int(row[cell]) for cell in ("n11", "n10", "n01", "n00")
        )
        for row in rows
    } == {
        "lof_hc": (1, 1, 3, 7),
        "missense": (1, 2, 3, 6),
        "splice_any": (0, 1, 4, 7),
    }
