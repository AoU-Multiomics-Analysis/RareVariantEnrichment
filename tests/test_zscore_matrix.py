import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


FIXTURES = Path(__file__).parent / "fixtures"


def run_export(tmp_path, *, bed=None, pcs=None, covariates=None, selected=1):
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps({"selection": {"selected_pc_count": selected}}))
    command = [
        sys.executable, "-m", "rare_variant_enrichment.cli", "export-zscore-matrix",
        "--phenotype-bed", str(bed or FIXTURES / "lof_pc_phenotypes.bed"),
        "--principal-components", str(pcs or FIXTURES / "principal_components.tsv"),
        "--selection-input", str(selection),
        "--matrix-output", str(tmp_path / "matrix.tsv.gz"),
        "--gene-qc-output", str(tmp_path / "gene_qc.tsv.gz"),
        "--summary-output", str(tmp_path / "summary.json"),
    ]
    if covariates is not None:
        command.extend(["--additional-covariates", str(covariates)])
    return subprocess.run(command, text=True, capture_output=True, check=False)


def read_matrix(tmp_path):
    with gzip.open(tmp_path / "matrix.tsv.gz", "rt") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    values = np.array([
        [np.nan if value == "NA" else float(value) for value in row[1:]]
        for row in rows[1:]
    ])
    return rows[0], [row[0] for row in rows[1:]], values


@pytest.mark.parametrize("selected", [0, 1])
def test_export_uses_selected_count_and_includes_noncoding_genes(tmp_path, selected):
    result = run_export(tmp_path, selected=selected)
    assert result.returncode == 0, result.stderr
    header, genes, values = read_matrix(tmp_path)
    assert header == ["gene_id", "S1", "S2", "S3", "S4", "S5", "S6"]
    assert genes == ["ENSG1", "ENSG2", "ENSG3"]
    # PC1 removes the alternating offset from the linear expression sequence.
    expected = (
        np.array([-2, -2, 0, 0, 2, 2]) / np.sqrt(8 / 3)
        if selected else np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]) / np.sqrt(35 / 12)
    )
    np.testing.assert_allclose(values, [expected, -expected, expected], atol=1e-12)
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["selected_pc_count"] == selected
    assert summary["gene_count"] == 3
    assert summary["included_gene_count"] == 3
    assert summary["sample_count"] == 6
    assert summary["residual_standard_deviation_ddof"] == 0


def test_export_aligns_covariates_collapses_features_and_preserves_missing_values(tmp_path):
    bed = tmp_path / "phenotypes.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write(
            "#chr\tstart\tend\tgene_id\tS3\tS1\tS6\tS2\tS5\tS4\tBED_ONLY\n"
            "chr1\t0\t1\tprotein_ENSG1.2\t1\t3\t2\t9\t0\t4\t9\n"
            "chr1\t1\t2\tENSG1.3\tNA\t5\t3\t1\t0\t8\t8\n"
            "chr1\t2\t3\tENSG2\tNA\t2\t1\t3\t7\t4\t6\n"
            "chr1\t3\t4\tENSG3\t8\t8\t8\t8\t8\t8\t8\n"
        )
    covariates = tmp_path / "covariates.tsv"
    covariates.write_text("sample_id\tage\nS6\t0\nS4\t2\nS2\t1\nS1\t0\nS3\t-1\n")
    result = run_export(tmp_path, bed=bed, covariates=covariates, selected=0)
    assert result.returncode == 0, result.stderr
    header, genes, values = read_matrix(tmp_path)
    assert header == ["gene_id", "S3", "S1", "S6", "S2", "S4"]
    assert genes == ["ENSG1", "ENSG2", "ENSG3"]
    # Independent least-squares reference, in BED sample order.
    design = np.column_stack([np.ones(5), [-1, 0, 0, 1, 2]])
    for index, expression in enumerate(([1, 3, 2, 1, 4], [np.nan, 2, 1, 3, 4])):
        expression = np.array(expression)
        usable = np.isfinite(expression)
        residual = expression[usable] - design[usable] @ np.linalg.lstsq(
            design[usable], expression[usable], rcond=None
        )[0]
        expected = (residual - residual.mean()) / residual.std(ddof=0)
        np.testing.assert_allclose(values[index, usable], expected, atol=1e-12)
    assert np.isnan(values[1, 0])
    assert np.isnan(values[2]).all()
    with gzip.open(tmp_path / "gene_qc.tsv.gz", "rt") as handle:
        qc = list(csv.DictReader(handle, delimiter="\t"))
    assert qc[2]["exclusion_reason"] == "invalid_or_zero_residual_sd"
    assert [row["status"] for row in qc] == ["included", "included", "excluded"]
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["additional_covariate_names"] == ["age"]
    assert summary["duplicate_feature_count"] == 1


@pytest.mark.parametrize("selected", [-1, True, 1.5, "1", None, 2])
def test_export_rejects_invalid_or_unselectable_pc_count(tmp_path, selected):
    result = run_export(tmp_path, selected=selected)
    assert result.returncode != 0
    assert "selected_pc_count" in result.stderr or "PC counts" in result.stderr
    assert not (tmp_path / "matrix.tsv.gz").exists()


def test_export_reports_unlocalized_cloud_input(tmp_path):
    result = run_export(tmp_path, bed="gs://bucket/phenotypes.bed.gz")
    assert result.returncode != 0
    assert "localization error" in result.stderr.lower()


@pytest.mark.parametrize(
    ("sample_count", "pc_values", "reason"),
    [
        (6, [-1, 1, -1, 1, -1, 1], "invalid_or_zero_residual_sd"),
        (6, [1, 1, 1, 1, 1, 1], "rank_deficiency"),
        (3, [-1, 1, -1], "insufficient_dof"),
    ],
)
def test_export_keeps_excluded_gene_as_na_row(tmp_path, sample_count, pc_values, reason):
    bed = tmp_path / "constant.bed"
    samples = [f"S{index}" for index in range(sample_count)]
    bed.write_text(
        "#chr\tstart\tend\tgene_id\t" + "\t".join(samples) + "\n"
        + "chr1\t0\t1\tENSG1\t" + "\t".join(["8"] * sample_count) + "\n"
    )
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\tPC2\n" + "".join(
        f"{sample}\t{value}\t{index}\n"
        for index, (sample, value) in enumerate(zip(samples, pc_values))
    ))
    result = run_export(tmp_path, bed=bed, pcs=pcs)
    assert result.returncode == 0, result.stderr
    assert read_matrix(tmp_path)[1] == ["ENSG1"]
    assert np.isnan(read_matrix(tmp_path)[2]).all()
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["exclusion_counts"][reason] == 1
    assert summary["included_gene_count"] == 0
