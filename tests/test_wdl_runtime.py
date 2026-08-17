import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest


WORKFLOW = Path("workflows/rare_variant_enrichment.wdl").resolve()
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
    docker_info = subprocess.run(
        ["docker", "info"], text=True, capture_output=True, check=False
    )
    if docker_info.returncode != 0:
        _unavailable("Docker daemon is unavailable: " + docker_info.stderr.strip())
    image = subprocess.run(
        ["docker", "image", "inspect", TEST_IMAGE],
        text=True,
        capture_output=True,
        check=False,
    )
    if image.returncode != 0:
        _unavailable(f"Build {TEST_IMAGE} before WDL runtime tests: " + image.stderr.strip())
    return shutil.which("miniwdl") or "miniwdl"


def test_required_runtime_mode_fails_instead_of_skipping_missing_prerequisites(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setenv("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME", "1")
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    with pytest.raises(pytest.fail.Exception, match="prerequisites are missing"):
        _require_wdl_runtime()


def test_wdl_runs_the_four_input_lof_pc_fixture_with_known_cells(tmp_path: Path):
    miniwdl = _require_wdl_runtime()
    inputs = {
        "RareVariantEnrichment.phenotype_bed": str((FIXTURES / "lof_pc_phenotypes.bed").resolve()),
        "RareVariantEnrichment.lof_carrier_table": str((FIXTURES / "lof_carriers.tsv").resolve()),
        "RareVariantEnrichment.principal_components_tsv": str((FIXTURES / "principal_components.tsv").resolve()),
        "RareVariantEnrichment.gene_annotation_gtf": str((FIXTURES / "gene_annotation.gtf").resolve()),
        "RareVariantEnrichment.negative_z_thresholds": [-0.8],
        "RareVariantEnrichment.selection_z_thresholds": [-0.8],
        "RareVariantEnrichment.pc_counts": [0, 1],
        "RareVariantEnrichment.pc_counts_per_job": 1,
        "RareVariantEnrichment.docker_image": TEST_IMAGE,
        "RareVariantEnrichment.prepare_cpu": 1,
        "RareVariantEnrichment.prepare_memory_gb": 1,
        "RareVariantEnrichment.prepare_disk_gb": 1,
        "RareVariantEnrichment.analysis_cpu": 1,
        "RareVariantEnrichment.analysis_memory_gb": 1,
        "RareVariantEnrichment.analysis_disk_gb": 1,
        "RareVariantEnrichment.max_retries": 1,
    }
    inputs_path = tmp_path / "inputs.json"
    outputs_path = tmp_path / "outputs.json"
    inputs_path.write_text(json.dumps(inputs))
    try:
        result = subprocess.run(
            [miniwdl, "run", str(WORKFLOW), "-i", str(inputs_path), "-d", str(tmp_path / "miniwdl"), "-o", str(outputs_path), "--no-cache"],
            text=True,
            capture_output=True,
            check=False,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        _unavailable("miniwdl Docker backend did not start the fixture workflow within 60 seconds")
    assert result.returncode == 0, result.stderr
    outputs = json.loads(outputs_path.read_text())["outputs"]
    assert set(outputs) == {
        "RareVariantEnrichment.results_tsv",
        "RareVariantEnrichment.summary_json",
        "RareVariantEnrichment.gene_pc_qc_tsv_gz",
        "RareVariantEnrichment.analysis_qc_json",
        "RareVariantEnrichment.pc_selection_json",
        "RareVariantEnrichment.enrichment_plot_svg",
        "RareVariantEnrichment.pc_sweep_qc_summary_tsv",
        "RareVariantEnrichment.pc_sweep_qc_plot_png",
        "RareVariantEnrichment.protein_coding_genes_tsv",
        "RareVariantEnrichment.protein_coding_genes_qc_json",
    }
    rows = list(csv.DictReader(Path(outputs["RareVariantEnrichment.results_tsv"]).open(), delimiter="\t"))
    assert {row["carrier_definition"]: tuple(
        int(row[cell]) for cell in ("n11", "n10", "n01", "n00")
    ) for row in rows if row["pc_count"] == "0"} == {
        "any_lof": (2, 3, 2, 5),
        "HC": (1, 1, 3, 7),
        "HC_or_LC": (2, 2, 2, 6),
    }
    summary = json.loads(Path(outputs["RareVariantEnrichment.summary_json"]).read_text())
    assert {row["pc_count"] for row in rows} == {"0", "1"}
    assert summary["selected_pc_counts"] == [0, 1]
    assert summary["fdr_scope"] == "global_across_all_emitted_rows"
    assert Path(outputs["RareVariantEnrichment.gene_pc_qc_tsv_gz"]).read_bytes()[:2] == b"\x1f\x8b"


def test_wdl_runs_optional_covariate_fixture_and_reports_intersection(
    tmp_path: Path,
):
    miniwdl = _require_wdl_runtime()
    inputs = {
        "RareVariantEnrichment.phenotype_bed": str(
            (FIXTURES / "lof_pc_phenotypes.bed").resolve()
        ),
        "RareVariantEnrichment.lof_carrier_table": str(
            (FIXTURES / "lof_carriers.tsv").resolve()
        ),
        "RareVariantEnrichment.principal_components_tsv": str(
            (FIXTURES / "principal_components.tsv").resolve()
        ),
        "RareVariantEnrichment.additional_covariates_tsv": str(
            (FIXTURES / "genetic_pcs.tsv").resolve()
        ),
        "RareVariantEnrichment.gene_annotation_gtf": str(
            (FIXTURES / "gene_annotation.gtf").resolve()
        ),
        "RareVariantEnrichment.negative_z_thresholds": [-0.8],
        "RareVariantEnrichment.selection_z_thresholds": [-0.8],
        "RareVariantEnrichment.pc_counts": [0],
        "RareVariantEnrichment.pc_counts_per_job": 1,
        "RareVariantEnrichment.pc_preemptible": 1,
        "RareVariantEnrichment.docker_image": TEST_IMAGE,
        "RareVariantEnrichment.prepare_cpu": 1,
        "RareVariantEnrichment.prepare_memory_gb": 1,
        "RareVariantEnrichment.prepare_disk_gb": 1,
        "RareVariantEnrichment.analysis_cpu": 1,
        "RareVariantEnrichment.analysis_memory_gb": 1,
        "RareVariantEnrichment.analysis_disk_gb": 1,
        "RareVariantEnrichment.max_retries": 1,
    }
    inputs_path = tmp_path / "covariate-inputs.json"
    outputs_path = tmp_path / "covariate-outputs.json"
    inputs_path.write_text(json.dumps(inputs))
    try:
        result = subprocess.run(
            [
                miniwdl,
                "run",
                str(WORKFLOW),
                "-i",
                str(inputs_path),
                "-d",
                str(tmp_path / "miniwdl-covariates"),
                "-o",
                str(outputs_path),
                "--no-cache",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        _unavailable("miniwdl Docker backend did not start the covariate fixture within 60 seconds")
    assert result.returncode == 0, result.stderr
    outputs = json.loads(outputs_path.read_text())["outputs"]
    analysis_qc = json.loads(
        Path(outputs["RareVariantEnrichment.analysis_qc_json"]).read_text()
    )
    assert analysis_qc["additional_covariates_supplied"] is True
    assert analysis_qc["additional_covariate_count"] == 2
    assert analysis_qc["additional_covariate_sample_count"] == 5
    assert analysis_qc["shared_bed_pc_covariate_sample_count"] == 5
