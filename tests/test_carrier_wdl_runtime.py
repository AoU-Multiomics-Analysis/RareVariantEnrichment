import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest


WORKFLOW = Path("workflows/extract_variant_carriers.wdl").resolve()
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


def test_wdl_runs_gene_matched_carrier_fixture(prepared_fixture, tmp_path: Path):
    miniwdl = _require_wdl_runtime()
    inputs = {
        "ExtractVariantCarriers.filtered_vcf": str(prepared_fixture.vcf_gz.resolve()),
        "ExtractVariantCarriers.filtered_vcf_tbi": str(prepared_fixture.vcf_tbi.resolve()),
        "ExtractVariantCarriers.transcript_annotations": str(
            prepared_fixture.vat_bgz.resolve()
        ),
        "ExtractVariantCarriers.chromosomes": ["chr1"],
        "ExtractVariantCarriers.annotation_chunk_size_bp": 250,
        "ExtractVariantCarriers.docker_image": TEST_IMAGE,
        "ExtractVariantCarriers.prepare_cpu": 1,
        "ExtractVariantCarriers.prepare_memory_gb": 1,
        "ExtractVariantCarriers.prepare_disk_gb": 1,
        "ExtractVariantCarriers.scatter_cpu": 1,
        "ExtractVariantCarriers.scatter_memory_gb": 1,
        "ExtractVariantCarriers.scatter_disk_gb": 1,
        "ExtractVariantCarriers.gather_cpu": 1,
        "ExtractVariantCarriers.gather_memory_gb": 1,
        "ExtractVariantCarriers.gather_disk_gb": 1,
        "ExtractVariantCarriers.max_retries": 1,
        "ExtractVariantCarriers.scatter_preemptible": 0,
    }
    inputs_path = tmp_path / "inputs.json"
    outputs_path = tmp_path / "outputs.json"
    inputs_path.write_text(json.dumps(inputs))
    result = subprocess.run(
        [
            miniwdl, "run", str(WORKFLOW), "-i", str(inputs_path),
            "-d", str(tmp_path / "miniwdl"), "-o", str(outputs_path), "--no-cache",
        ],
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    outputs = json.loads(outputs_path.read_text())["outputs"]
    assert set(outputs) == {
        "ExtractVariantCarriers.variant_carrier_audit_tsv_gz",
        "ExtractVariantCarriers.variant_carriers_tsv_gz",
        "ExtractVariantCarriers.variant_carriers_qc_json",
        "ExtractVariantCarriers.chromosome_qc_jsons",
        "ExtractVariantCarriers.transcript_schema_json",
        "ExtractVariantCarriers.generated_or_validated_transcript_annotations_tbi",
        "ExtractVariantCarriers.transcript_index_provenance",
    }
    carrier_path = Path(outputs["ExtractVariantCarriers.variant_carriers_tsv_gz"])
    assert carrier_path.read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(carrier_path, "rt", newline="") as handle:
        classes = {row["variant_class"] for row in csv.DictReader(handle, delimiter="\t")}
    assert {"lof_hc", "lof_hc_or_lc", "missense"} <= classes
    assert outputs["ExtractVariantCarriers.transcript_index_provenance"] == "generated"
    qc = json.loads(
        Path(outputs["ExtractVariantCarriers.variant_carriers_qc_json"]).read_text()
    )
    assert qc["quality_or_frequency_filters_applied"] is False
    assert qc["audit_row_count"] == qc["chromosome_qc"][0]["carrier_audit_rows"]
