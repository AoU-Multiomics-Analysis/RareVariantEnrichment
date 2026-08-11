import json
import gzip
from pathlib import Path
import subprocess
import sys

import pytest

from rare_variant_enrichment.cli import build_parser
from rare_variant_enrichment.phenotypes import classify_outlier, prepare_phenotypes


def test_classify_outlier_supports_all_tail_modes():
    assert classify_outlier(-3.0, 2.5, "absolute")
    assert classify_outlier(3.0, 2.5, "positive")
    assert classify_outlier(-3.0, 2.5, "negative")
    assert not classify_outlier(-3.0, 2.5, "positive")


def test_prepare_phenotypes_uses_tss_end_and_shared_samples(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\n"
        "chr1\t99\t100\tGENE1\t3.0\t0.0\tNA\n"
        "chr2\t199\t200\tGENE2\t-4.0\t1.0\t2.5\n"
    )
    vcf_samples = tmp_path / "vcf_samples.txt"
    vcf_samples.write_text("S1\nS3\nS4\n")

    feature_output = tmp_path / "features.tsv"
    sample_output = tmp_path / "shared_samples.txt"
    qc_output = tmp_path / "phenotype_qc.json"
    prepare_phenotypes(
        bed, vcf_samples, ["chr1", "chr2"], [2.0, 3.0], "absolute",
        feature_output, sample_output, qc_output,
    )

    assert feature_output.read_text().splitlines()[1] == "chr1\t100\tGENE1"
    assert sample_output.read_text().splitlines() == ["S1", "S3"]
    qc = json.loads(qc_output.read_text())
    assert qc["shared_sample_count"] == 2
    assert qc["bed_only_sample_count"] == 1
    assert qc["vcf_only_sample_count"] == 1
    assert qc["bed_sample_count"] == 3
    assert qc["vcf_sample_count"] == 3
    assert qc["selected_chromosomes"] == ["chr1", "chr2"]
    assert "bed_only_samples" not in qc
    assert "vcf_only_samples" not in qc
    assert qc["non_missing_observations"] == 3


def test_prepare_phenotypes_rejects_non_unit_tss_interval(tmp_path: Path):
    bed = tmp_path / "bad.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t90\t100\tGENE1\t2.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    with pytest.raises(ValueError, match="one base"):
        prepare_phenotypes(bed, samples, ["chr1"], [2.0], "absolute",
                           tmp_path / "f.tsv", tmp_path / "s.txt", tmp_path / "q.json")


@pytest.mark.parametrize("duplicate", ["feature", "sample"])
def test_prepare_phenotypes_rejects_duplicate_ids(tmp_path: Path, duplicate: str):
    if duplicate == "feature":
        bed_text = (
            "#chr\tstart\tend\tgene_id\tS1\n"
            "chr1\t0\t1\tGENE1\t1.0\n"
            "chr1\t1\t2\tGENE1\t2.0\n"
        )
    else:
        bed_text = "#chr\tstart\tend\tgene_id\tS1\tS1\nchr1\t0\t1\tGENE1\t1.0\t2.0\n"
    bed = tmp_path / "duplicate.bed"
    bed.write_text(bed_text)
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match=f"Duplicate {duplicate}"):
        prepare_phenotypes(bed, samples, ["chr1"], [2.0], "absolute",
                           tmp_path / "f.tsv", tmp_path / "s.txt", tmp_path / "q.json")


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), "invalid"])
def test_prepare_phenotypes_rejects_non_finite_or_non_numeric_thresholds(
    tmp_path: Path, threshold: object
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t0\t1\tGENE1\t1.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="threshold"):
        prepare_phenotypes(bed, samples, ["chr1"], [threshold], "absolute",  # type: ignore[list-item]
                           tmp_path / "f.tsv", tmp_path / "s.txt", tmp_path / "q.json")


def test_prepare_phenotypes_rejects_negative_z_threshold(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t0\t1\tGENE1\t1.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="z-score thresholds must be non-negative"):
        prepare_phenotypes(
            bed,
            samples,
            ["chr1"],
            [-1.0],
            "absolute",
            tmp_path / "f.tsv",
            tmp_path / "s.txt",
            tmp_path / "q.json",
        )


def test_prepare_phenotypes_streams_gzip_and_counts_shared_observations(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write(
            "#chr\tstart\tend\tgene_id\tS1\tS2\n"
            "chr1\t0\t1\tGENE1\t2.0\tNA\n"
            "chr2\t1\t2\tGENE2\t.\tNaN\n"
            "chr3\t2\t3\tGENE3\t9.0\t9.0\n"
        )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    features = tmp_path / "features.tsv"
    qc_output = tmp_path / "qc.json"

    prepare_phenotypes(bed, samples, ["chr1", "chr2"], [2.0], "absolute",
                       features, tmp_path / "shared.txt", qc_output)

    assert features.read_text().splitlines() == [
        "chrom\ttss\tfeature_id", "chr1\t1\tGENE1", "chr2\t2\tGENE2"
    ]
    qc = json.loads(qc_output.read_text())
    assert qc["non_missing_observations"] == 1
    assert qc["outlier_observations"] == {"2.0": 1}
    assert qc["input_feature_count"] == 3
    assert qc["selected_feature_count"] == 2
    assert qc["unselected_feature_count"] == 1
    assert qc["missing_z_observations"] == 1


def test_prepare_phenotypes_rejects_missing_requested_chromosome(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t0\t1\tGENE1\t1.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="chr2"):
        prepare_phenotypes(bed, samples, ["chr1", "chr2"], [2.0], "absolute",
                           tmp_path / "f.tsv", tmp_path / "s.txt", tmp_path / "q.json")


def test_prepare_phenotypes_cli_parses_lists_and_writes_outputs(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t0\t1\tGENE1\t3.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    features = tmp_path / "features.tsv"
    shared = tmp_path / "shared.txt"
    qc = tmp_path / "qc.json"
    command = [
        sys.executable, "-m", "rare_variant_enrichment.cli", "prepare-phenotypes",
        "--phenotype-bed", str(bed), "--vcf-samples", str(samples),
        "--chromosomes", "chr1", "--z-thresholds", "2.5", "--tail", "absolute",
        "--feature-output", str(features), "--sample-output", str(shared),
        "--qc-output", str(qc),
    ]

    result = subprocess.run(command, text=True, capture_output=True, check=False)

    assert result.returncode == 0, result.stderr
    assert features.read_text().splitlines() == ["chrom\ttss\tfeature_id", "chr1\t1\tGENE1"]
    parsed = build_parser().parse_args(command[3:])
    assert parsed.chromosomes == ["chr1"]
    assert parsed.z_thresholds == [2.5]
