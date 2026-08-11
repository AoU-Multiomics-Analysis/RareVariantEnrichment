import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest


WORKFLOW = Path("workflows/rare_variant_enrichment.wdl").resolve()
TEST_IMAGE = os.environ.get(
    "RARE_VARIANT_ENRICHMENT_TEST_IMAGE", "rare-variant-enrichment:test"
)


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
        _unavailable(
            f"Build {TEST_IMAGE} before running WDL integration tests: "
            + image.stderr.strip()
        )
    return shutil.which("miniwdl") or "miniwdl"


def _unavailable(message: str) -> None:
    if os.environ.get("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME") == "1":
        pytest.fail(message)
    pytest.skip(message)


def test_required_runtime_mode_fails_instead_of_skipping_missing_prerequisites(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setenv("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME", "1")
    monkeypatch.setattr(shutil, "which", lambda _name: None)

    with pytest.raises(pytest.fail.Exception, match="prerequisites are missing"):
        _require_wdl_runtime()


@pytest.mark.parametrize(
    ("exact_ac", "cumulative_ac_max", "expected_kind", "supply_index"),
    [
        ([1, 2], [], "exact", False),
        ([], [1, 2], "cumulative", True),
    ],
    ids=["exact-only-generated-index", "cumulative-only-supplied-index"],
)
def test_wdl_runs_single_ac_family_and_optional_index_modes(
    prepared_fixture,
    tmp_path: Path,
    exact_ac: list[int],
    cumulative_ac_max: list[int],
    expected_kind: str,
    supply_index: bool,
):
    miniwdl = _require_wdl_runtime()
    inputs = {
        "RareVariantEnrichment.phenotype_bed": str(prepared_fixture.bed.resolve()),
        "RareVariantEnrichment.rare_variant_vcf": str(prepared_fixture.vcf_gz.resolve()),
        "RareVariantEnrichment.chromosomes": ["chr1"],
        "RareVariantEnrichment.z_thresholds": [2.0, 3.0],
        "RareVariantEnrichment.exact_allele_counts": exact_ac,
        "RareVariantEnrichment.cumulative_allele_count_maxima": cumulative_ac_max,
        "RareVariantEnrichment.distance_thresholds_bp": [10, 100],
        "RareVariantEnrichment.outlier_tail": "absolute",
        "RareVariantEnrichment.docker_image": TEST_IMAGE,
        "RareVariantEnrichment.prepare_cpu": 1,
        "RareVariantEnrichment.prepare_memory_gb": 1,
        "RareVariantEnrichment.prepare_disk_gb": 1,
        "RareVariantEnrichment.scatter_cpu": 1,
        "RareVariantEnrichment.scatter_memory_gb": 1,
        "RareVariantEnrichment.scatter_disk_gb": 1,
        "RareVariantEnrichment.gather_cpu": 1,
        "RareVariantEnrichment.gather_memory_gb": 1,
        "RareVariantEnrichment.gather_disk_gb": 1,
        "RareVariantEnrichment.max_retries": 2,
        "RareVariantEnrichment.publish_carrier_audit": supply_index,
    }
    if supply_index:
        inputs["RareVariantEnrichment.rare_variant_vcf_tbi"] = str(
            prepared_fixture.vcf_tbi.resolve()
        )

    inputs_path = tmp_path / "inputs.json"
    inputs_path.write_text(json.dumps(inputs))
    outputs_path = tmp_path / "outputs.json"
    run_directory = tmp_path / "miniwdl-run"
    run_directory.mkdir()
    result = subprocess.run(
        [
            miniwdl,
            "run",
            str(WORKFLOW),
            "-i",
            str(inputs_path),
            "-d",
            f"{run_directory}/.",
            "-o",
            str(outputs_path),
            "--no-cache",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    outputs = json.loads(outputs_path.read_text())["outputs"]
    enrichment_path = Path(outputs["RareVariantEnrichment.enrichment_tsv"])
    rows = list(csv.DictReader(enrichment_path.open(), delimiter="\t"))
    assert len(rows) == 2 * 2 * 2
    assert {row["z_threshold"] for row in rows} == {"2.0", "3.0"}
    assert {row["distance_bp"] for row in rows} == {"10", "100"}
    assert {row["ac_kind"] for row in rows} == {expected_kind}
    assert {row["ac_class"] for row in rows} == {
        f"AC{'=' if expected_kind == 'exact' else '<='}{value}" for value in (1, 2)
    }
    assert Path(
        outputs["RareVariantEnrichment.generated_or_validated_vcf_tbi"]
    ).is_file()
    summary = json.loads(
        Path(outputs["RareVariantEnrichment.enrichment_json"]).read_text()
    )
    assert summary["provenance"]["vcf_index"] == (
        "supplied" if supply_index else "generated"
    )
    assert summary["provenance"]["max_retries"] == 2
    assert summary["provenance"]["selected_chromosomes"] == ["chr1"]
    phenotype_qc = json.loads(
        Path(outputs["RareVariantEnrichment.phenotype_qc_json"]).read_text()
    )
    assert phenotype_qc["selected_feature_count"] == 4
    audit_output = outputs["RareVariantEnrichment.carrier_minimum_distances_tsv"]
    assert (audit_output is not None) is supply_index

    by_key = {
        (row["z_threshold"], row["ac_class"], row["distance_bp"]): row
        for row in rows
    }
    expected_class = "AC=1" if expected_kind == "exact" else "AC<=1"
    hand_checked = by_key[("2.0", expected_class, "10")]
    assert (
        hand_checked["total_observations"],
        hand_checked["n11"],
        hand_checked["n10"],
        hand_checked["n01"],
        hand_checked["n00"],
    ) == ("11", "3", "0", "0", "8")
