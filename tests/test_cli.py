import subprocess
import sys
from pathlib import Path

import pytest

from rare_variant_enrichment.cli import build_parser, parse_csv_ints


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


def test_parse_csv_ints_decodes_empty_wdl_array():
    assert parse_csv_ints("") == []


@pytest.mark.parametrize(
    ("command", "exact_ac", "cumulative_ac_max"),
    [
        ("classify-chromosome", "", "1"),
        ("classify-chromosome", "1", ""),
        ("calculate", "", "1"),
        ("calculate", "1", ""),
    ],
)
def test_ac_cli_options_allow_one_empty_family(
    command: str, exact_ac: str, cumulative_ac_max: str
):
    common = [
        "--exact-ac",
        exact_ac,
        "--cumulative-ac-max",
        cumulative_ac_max,
    ]
    if command == "classify-chromosome":
        arguments = [
            command,
            "--vcf",
            "variants.vcf.gz",
            "--features",
            "features.tsv",
            "--shared-samples",
            "shared_samples.txt",
            "--chromosome",
            "chr1",
            *common,
            "--max-distance",
            "1000",
            "--carrier-output",
            "carriers.tsv",
            "--regions-output",
            "regions.bed",
            "--qc-output",
            "qc.json",
        ]
    else:
        arguments = [
            command,
            "--phenotype-bed",
            "phenotypes.bed",
            "--shared-samples",
            "shared_samples.txt",
            "--carriers",
            "carriers.tsv",
            "--features",
            "features.tsv",
            *common,
            "--z-thresholds",
            "2.0",
            "--distance-thresholds",
            "1000",
            "--tail",
            "absolute",
            "--output-tsv",
            "enrichment.tsv",
            "--output-json",
            "enrichment.json",
        ]

    parsed = build_parser().parse_args(arguments)

    assert parsed.exact_ac == ([] if exact_ac == "" else [1])
    assert parsed.cumulative_ac_max == ([] if cumulative_ac_max == "" else [1])
