import subprocess
import sys
from pathlib import Path


def test_cli_lists_workflow_subcommands():
    result = subprocess.run(
        [sys.executable, "-m", "rare_variant_enrichment.cli", "--help"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    for command in ("prepare-phenotypes", "classify-chromosome", "gather", "calculate"):
        assert command in result.stdout


def test_gather_cli_writes_aggregated_outputs(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text("sample_id\tfeature_id\tac_class\tminimum_distance_bp\nS1\tGENE1\tAC=1\t4\n")
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1","extracted_records":1}')
    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "qc.tsv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rare_variant_enrichment.cli",
            "gather",
            "--carrier-input",
            str(carrier),
            "--qc-input",
            str(qc),
            "--carrier-output",
            str(carrier_output),
            "--qc-output",
            str(qc_output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "S1\tGENE1\tAC=1\t4" in carrier_output.read_text().splitlines()
    assert qc_output.read_text().splitlines() == ["chromosome\textracted_records", "chr1\t1"]
