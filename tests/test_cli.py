import subprocess
import sys


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
