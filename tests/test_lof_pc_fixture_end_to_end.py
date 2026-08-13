import csv
import json
from pathlib import Path
import subprocess
import sys


FIXTURES = Path(__file__).parent / "fixtures"


def test_new_cli_commands_run_four_input_fixture_with_known_fisher_cells(tmp_path: Path):
    genes = tmp_path / "protein-coding-genes.tsv"
    genes_qc = tmp_path / "protein-coding-genes.qc.json"
    prepare = subprocess.run(
        [
            sys.executable,
            "-m",
            "rare_variant_enrichment.cli",
            "prepare-protein-coding-genes",
            "--gtf",
            str(FIXTURES / "gene_annotation.gtf"),
            "--genes-output",
            str(genes),
            "--qc-output",
            str(genes_qc),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert prepare.returncode == 0, prepare.stderr
    assert genes.read_text().splitlines() == ["gene_id", "ENSG1", "ENSG2"]

    results = tmp_path / "results.tsv"
    summary = tmp_path / "summary.json"
    gene_qc = tmp_path / "gene-pc-qc.tsv.gz"
    analysis_qc = tmp_path / "analysis-qc.json"
    calculate = subprocess.run(
        [
            sys.executable,
            "-m",
            "rare_variant_enrichment.cli",
            "lof-pc-enrichment",
            "--phenotype-bed",
            str(FIXTURES / "lof_pc_phenotypes.bed"),
            "--lof-carriers",
            str(FIXTURES / "lof_carriers.tsv"),
            "--principal-components",
            str(FIXTURES / "principal_components.tsv"),
            "--protein-coding-genes",
            str(genes),
            "--negative-z-thresholds",
            "-0.8",
            "--pc-counts",
            "0",
            "--results-output",
            str(results),
            "--summary-output",
            str(summary),
            "--gene-pc-qc-output",
            str(gene_qc),
            "--analysis-qc-output",
            str(analysis_qc),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert calculate.returncode == 0, calculate.stderr

    rows = list(csv.DictReader(results.open(), delimiter="\t"))
    cells = {
        row["carrier_definition"]: tuple(
            int(row[column]) for column in ("n11", "n10", "n01", "n00")
        )
        for row in rows
    }
    assert cells == {
        "any_lof": (2, 3, 2, 5),
        "HC": (1, 1, 3, 7),
        "HC_or_LC": (2, 2, 2, 6),
    }
    assert json.loads(summary.read_text())["fdr_scope"] == "global_across_all_emitted_rows"
    assert gene_qc.read_bytes()[:2] == b"\x1f\x8b"
    assert json.loads(analysis_qc.read_text())["per_pc"]["0"][
        "carrier_observations"
    ] == {"HC": 2, "HC_or_LC": 4, "any_lof": 5}
