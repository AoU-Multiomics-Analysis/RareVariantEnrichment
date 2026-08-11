import csv
import json
from pathlib import Path
import shutil

import pytest

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.phenotypes import prepare_phenotypes
from rare_variant_enrichment.statistics import calculate_enrichment
from rare_variant_enrichment.variants import classify_chromosome


requires_htslib = pytest.mark.skipif(
    shutil.which("tabix") is None or shutil.which("bgzip") is None,
    reason="htslib executables are required",
)


@requires_htslib
def test_miniature_pipeline_emits_hand_checked_threshold_combinations(
    prepared_fixture, tmp_path: Path
):
    feature_tsv = tmp_path / "features.tsv"
    shared = tmp_path / "shared.txt"
    phenotype_qc = tmp_path / "phenotype_qc.json"
    prepare_phenotypes(
        prepared_fixture.bed,
        prepared_fixture.samples,
        ["chr1"],
        [2.0, 3.0],
        "absolute",
        feature_tsv,
        shared,
        phenotype_qc,
    )

    assert shared.read_text().splitlines() == ["S1", "S2", "S3"]
    phenotype_summary = json.loads(phenotype_qc.read_text())
    assert phenotype_summary["bed_only_sample_count"] == 1
    assert phenotype_summary["non_missing_observations"] == 11
    assert phenotype_summary["outlier_observations"] == {"2.0": 3, "3.0": 2}

    chromosome_carriers = tmp_path / "chr1.carriers.tsv"
    regions = tmp_path / "chr1.regions.bed"
    chromosome_qc = tmp_path / "chr1.qc.json"
    classify_chromosome(
        prepared_fixture.vcf_gz,
        feature_tsv,
        shared,
        "chr1",
        [1, 2],
        [1, 2],
        100,
        chromosome_carriers,
        regions,
        chromosome_qc,
    )

    chromosome_summary = json.loads(chromosome_qc.read_text())
    assert chromosome_summary["missing_genotypes"] == 1
    assert chromosome_summary["tabix_query_count"] == 1
    assert "S1\tGENE4\tAC=1\t0" in chromosome_carriers.read_text().splitlines()
    assert prepared_fixture.vcf_tbi.is_file()

    all_carriers = tmp_path / "carriers.tsv"
    all_qc = tmp_path / "chromosome_qc.tsv"
    gather_outputs([chromosome_carriers], [chromosome_qc], all_carriers, all_qc)

    enrichment = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"
    calculate_enrichment(
        prepared_fixture.bed,
        shared,
        all_carriers,
        feature_tsv,
        [1, 2],
        [1, 2],
        [2.0, 3.0],
        [10, 100],
        "absolute",
        enrichment,
        summary,
    )

    rows = list(csv.DictReader(enrichment.open(), delimiter="\t"))
    assert len(rows) == 2 * 4 * 2
    by_key = {
        (row["z_threshold"], row["ac_class"], row["distance_bp"]): row
        for row in rows
    }
    assert set(row["ac_class"] for row in rows) == {"AC=1", "AC=2", "AC<=1", "AC<=2"}

    exact_singletons = by_key[("2.0", "AC=1", "10")]
    assert (
        exact_singletons["total_observations"],
        exact_singletons["n11"],
        exact_singletons["n10"],
        exact_singletons["n01"],
        exact_singletons["n00"],
    ) == ("11", "3", "0", "0", "8")
    assert float(exact_singletons["fisher_p_value"]) == pytest.approx(1 / 165)
    assert float(exact_singletons["fisher_fdr_bh"]) == pytest.approx(0.019393939393939394)

    cumulative_doubletons = by_key[("2.0", "AC<=2", "100")]
    assert (
        cumulative_doubletons["n11"],
        cumulative_doubletons["n10"],
        cumulative_doubletons["n01"],
        cumulative_doubletons["n00"],
    ) == ("3", "0", "4", "4")
    assert float(cumulative_doubletons["fisher_fdr_bh"]) == pytest.approx(
        0.3151515151515151
    )

    stricter_singletons = by_key[("3.0", "AC=1", "10")]
    assert (
        stricter_singletons["n11"],
        stricter_singletons["n10"],
        stricter_singletons["n01"],
        stricter_singletons["n00"],
    ) == ("2", "0", "1", "8")
    assert float(stricter_singletons["fisher_fdr_bh"]) == pytest.approx(
        0.08727272727272727
    )
    assert all(
        0.0 <= float(row["fisher_p_value"]) <= float(row["fisher_fdr_bh"]) <= 1.0
        for row in rows
    )

    run_summary = json.loads(summary.read_text())
    assert run_summary["missing_z_observations"] == 1
    assert run_summary["emitted_rows"] == 16


def test_pipeline_rejects_duplicate_normalized_z_thresholds_during_preparation(
    tmp_path: Path,
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\n"
        "chr1\t99\t100\tGENE1\t3.5\n"
    )
    samples = tmp_path / "vcf_samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="z-score thresholds must be unique"):
        prepare_phenotypes(
            bed,
            samples,
            ["chr1"],
            [2, 2.0],
            "absolute",
            tmp_path / "features.tsv",
            tmp_path / "shared.txt",
            tmp_path / "phenotype_qc.json",
        )


@pytest.mark.parametrize("stage", ["prepare", "calculate"])
def test_pipeline_rejects_boolean_z_thresholds_consistently(stage: str, tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\n"
        "chr1\t99\t100\tGENE1\t3.5\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="z-score thresholds must be finite numeric values"):
        if stage == "prepare":
            prepare_phenotypes(
                bed,
                samples,
                ["chr1"],
                [True],
                "absolute",
                tmp_path / "features.tsv",
                tmp_path / "shared.txt",
                tmp_path / "phenotype_qc.json",
            )
        else:
            carriers = tmp_path / "carriers.tsv"
            carriers.write_text(
                "sample_id\tfeature_id\tac_class\tminimum_distance_bp\n"
            )
            features = tmp_path / "features.tsv"
            features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
            calculate_enrichment(
                bed,
                samples,
                carriers,
                features,
                [1],
                [],
                [True],
                [100],
                "absolute",
                tmp_path / "enrichment.tsv",
                tmp_path / "summary.json",
            )
