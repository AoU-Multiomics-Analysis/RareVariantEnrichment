import csv
import json
from pathlib import Path
import subprocess
import sys

import pytest

from rare_variant_enrichment.statistics import (
    benjamini_hochberg,
    calculate_enrichment,
    fisher_exact_two_sided,
)


def test_fisher_exact_matches_hand_checked_table():
    assert fisher_exact_two_sided(1, 9, 11, 3) == pytest.approx(0.002759456, rel=1e-6)


def test_fisher_exact_preserves_precision_at_billion_scale_small_margin():
    p_value = fisher_exact_two_sided(1, 19_999_999, 9, 7_979_999_991)
    assert p_value == pytest.approx(0.024720616835182417, rel=1e-12, abs=0.0)


def test_fisher_exact_handles_all_zero_cells_and_rejects_negative_cells():
    assert fisher_exact_two_sided(0, 0, 0, 0) == 1.0
    with pytest.raises(ValueError, match="non-negative integers"):
        fisher_exact_two_sided(-1, 1, 1, 1)


def test_bh_is_monotone_in_rank_order():
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03])
    assert adjusted == pytest.approx([0.03, 0.04, 0.04])


def test_calculate_enrichment_counts_distance_specific_carriers(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\n"
        "chr1\t99\t100\tGENE1\t3\t0\t-3\t0\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\nS3\nS4\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t10\n"
        "S2\tGENE1\tAC=1\tbaseline\tall_rare_variants\t100\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    calculate_enrichment(
        bed,
        samples,
        carriers,
        features,
        [1],
        [],
        [2.0],
        [50, 100],
        "absolute",
        output,
        summary,
    )

    reader = csv.DictReader(output.open(), delimiter="\t")
    rows = list(reader)
    assert reader.fieldnames == [
        "z_threshold",
        "tail",
        "annotation_family",
        "annotation_class",
        "ac_class",
        "ac_kind",
        "ac_value",
        "distance_bp",
        "total_observations",
        "outlier_observations",
        "nonoutlier_observations",
        "n11",
        "n10",
        "n01",
        "n00",
        "outlier_carrier_rate",
        "nonoutlier_carrier_rate",
        "carrier_rate_ratio",
        "odds_ratio",
        "odds_ratio_corrected_0_5",
        "fisher_p_value",
        "fisher_fdr_bh",
    ]
    assert rows[0]["distance_bp"] == "50"
    assert (rows[0]["n11"], rows[0]["n10"], rows[0]["n01"], rows[0]["n00"]) == (
        "1",
        "1",
        "0",
        "2",
    )
    assert (rows[1]["n11"], rows[1]["n10"], rows[1]["n01"], rows[1]["n00"]) == (
        "1",
        "1",
        "1",
        "1",
    )
    assert "screening" in json.loads(summary.read_text())["statistical_limitation"].lower()


def test_calculate_enrichment_stratifies_annotations_with_one_global_bh_correction(
    tmp_path: Path,
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\tS5\tS6\n"
        "chr1\t99\t100\tGENE1\t3\t3\t3\t0\t0\t0\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\nS3\nS4\nS5\nS6\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t0\n"
        "S4\tGENE1\tAC=1\tbaseline\tall_rare_variants\t0\n"
        "S1\tGENE1\tAC=1\tconsequence\tstop_gained\t0\n"
        "S2\tGENE1\tAC=1\tconsequence\tstop_gained\t0\n"
        "S4\tGENE1\tAC=1\tconsequence\tmissense_variant\t0\n"
        "S5\tGENE1\tAC=1\tconsequence\tmissense_variant\t0\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
    output = tmp_path / "enrichment.tsv"

    calculate_enrichment(
        bed,
        samples,
        carriers,
        features,
        [1],
        [],
        [2.0],
        [0],
        "absolute",
        output,
        tmp_path / "summary.json",
        consequence_classes=["stop_gained", "missense_variant"],
    )

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    by_annotation = {
        (row["annotation_family"], row["annotation_class"]): row for row in rows
    }
    assert list(by_annotation) == [
        ("baseline", "all_rare_variants"),
        ("consequence", "stop_gained"),
        ("consequence", "missense_variant"),
    ]
    assert {
        annotation: (row["n11"], row["n10"], row["n01"], row["n00"])
        for annotation, row in by_annotation.items()
    } == {
        ("baseline", "all_rare_variants"): ("1", "2", "1", "2"),
        ("consequence", "stop_gained"): ("2", "1", "0", "3"),
        ("consequence", "missense_variant"): ("0", "3", "2", "1"),
    }
    assert {
        annotation: float(row["fisher_fdr_bh"]) for annotation, row in by_annotation.items()
    } == pytest.approx(
        {
            ("baseline", "all_rare_variants"): 1.0,
            ("consequence", "stop_gained"): 0.6,
            ("consequence", "missense_variant"): 0.6,
        }
    )


def test_calculate_enrichment_emits_configured_zero_carrier_annotation_rows_and_provenance(
    tmp_path: Path,
):
    """Removing configured annotation strata or provenance must fail this contract."""
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\n"
        "chr1\t99\t100\tENSG1.1\t3\t0\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\n")
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tENSG1.1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\t"
        "annotation_class\tminimum_distance_bp\n"
        "S1\tENSG1.1\tAC=1\tbaseline\tall_rare_variants\t0\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary_path = tmp_path / "enrichment.json"

    calculate_enrichment(
        bed,
        samples,
        carriers,
        features,
        [1],
        [],
        [2.0],
        [100],
        "absolute",
        output,
        summary_path,
        consequence_classes=["missense_variant", "stop_gained"],
        loftee_enabled=True,
        vat_index_provenance="generated",
        maximum_gvs_maf=0.01,
        annotation_chunk_size_bp=10_000_000,
    )

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert [(row["annotation_family"], row["annotation_class"], row["n11"]) for row in rows] == [
        ("baseline", "all_rare_variants", "1"),
        ("consequence", "stop_gained", "0"),
        ("consequence", "missense_variant", "0"),
        ("loftee", "HC", "0"),
        ("loftee", "LC", "0"),
    ]
    fdr_by_p_value = sorted(
        (float(row["fisher_p_value"]), float(row["fisher_fdr_bh"])) for row in rows
    )
    assert [fdr for _, fdr in fdr_by_p_value] == sorted(fdr for _, fdr in fdr_by_p_value)
    summary = json.loads(summary_path.read_text())
    assert summary["analysis_parameters"] == {
        "annotation_chunk_size_bp": 10_000_000,
        "annotation_classes": [
            {"family": "baseline", "label": "all_rare_variants"},
            {"family": "consequence", "label": "stop_gained"},
            {"family": "consequence", "label": "missense_variant"},
            {"family": "loftee", "label": "HC"},
            {"family": "loftee", "label": "LC"},
        ],
        "consequence_classes": ["stop_gained", "missense_variant"],
        "cumulative_allele_count_maxima": [],
        "distance_thresholds_bp": [100],
        "exact_allele_counts": [1],
        "loftee_enabled": True,
        "maximum_gvs_maf": 0.01,
        "outlier_tail": "absolute",
        "z_thresholds": [2.0],
    }
    assert summary["provenance"] | {"software_versions": None} == {
        "annotation_chunk_size_bp": 10_000_000,
        "container_image": None,
        "consequence_classes": ["stop_gained", "missense_variant"],
        "loftee_enabled": True,
        "maximum_gvs_maf": 0.01,
        "max_retries": 0,
        "selected_chromosomes": ["chr1"],
        "severity_order_version": "Ensembl release 116",
        "software_versions": None,
        "vat_index": "generated",
        "vcf_index": "unknown",
    }


@pytest.mark.parametrize(
    ("keyword", "value", "message"),
    [
        ("consequence_classes", ["not_an_ensembl_consequence"], "Unknown consequence classes"),
        ("maximum_gvs_maf", float("nan"), "maximum_gvs_maf must be a finite number from 0 to 0.5"),
        ("maximum_gvs_maf", 0.500_001, "maximum_gvs_maf must be a finite number from 0 to 0.5"),
        ("annotation_chunk_size_bp", 0, "annotation_chunk_size_bp must be a positive integer"),
        ("vat_index_provenance", "invalid", "vat_index_provenance must be generated, supplied, or unknown"),
    ],
)
def test_calculate_enrichment_validates_annotation_configuration(
    tmp_path: Path, keyword: str, value: object, message: str
):
    """Removing the annotation-configuration validators must fail this contract."""
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    with pytest.raises(ValueError, match=message):
        calculate_enrichment(
            bed,
            samples,
            carriers,
            features,
            [1],
            [],
            [2.0],
            [100],
            "absolute",
            tmp_path / "enrichment.tsv",
            tmp_path / "summary.json",
            **{keyword: value},  # type: ignore[arg-type]
        )


def test_calculate_enrichment_rejects_unconfigured_carrier_annotation(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tconsequence\tstop_gained\t0\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    with pytest.raises(ValueError, match="Carrier annotation class is not configured"):
        calculate_enrichment(
            bed,
            samples,
            carriers,
            features,
            [1],
            [],
            [2.0],
            [0],
            "absolute",
            tmp_path / "enrichment.tsv",
            tmp_path / "summary.json",
        )


def test_calculate_enrichment_uses_exact_prepared_feature_set(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\n"
        "chr1\t99\t100\tGENE1\t3\t0\n"
        "chr2\t199\t200\tGENE2\t3\t0\n"
    )
    selected_features = tmp_path / "features.tsv"
    selected_features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"

    calculate_enrichment(
        bed,
        samples,
        carriers,
        selected_features,
        [1],
        [],
        [2.0],
        [100],
        "absolute",
        output,
        summary,
    )

    row = next(csv.DictReader(output.open(), delimiter="\t"))
    assert (row["total_observations"], row["outlier_observations"]) == ("2", "1")
    assert json.loads(summary.read_text())["feature_count"] == 1


def test_calculate_enrichment_embeds_counts_and_reproducibility_provenance(
    tmp_path: Path,
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    selected_features = tmp_path / "features.tsv"
    selected_features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
    )
    phenotype_qc = tmp_path / "phenotype_qc.json"
    phenotype_qc.write_text(
        json.dumps(
            {
                "bed_only_sample_count": 2,
                "shared_sample_count": 1,
                "vcf_only_sample_count": 3,
                "selected_chromosomes": ["chr1"],
            }
        )
    )
    chromosome_qc = tmp_path / "chromosome_qc.tsv"
    chromosome_qc.write_text(
        "chromosome\talt_alleles\tboundary_variant_feature_pairs\tinfo_ac_alt_alleles\n"
        "chr1\t4\t2\t3\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary_path = tmp_path / "summary.json"

    calculate_enrichment(
        bed,
        samples,
        carriers,
        selected_features,
        [1],
        [],
        [2.0],
        [100],
        "absolute",
        output,
        summary_path,
        phenotype_qc_path=phenotype_qc,
        chromosome_qc_path=chromosome_qc,
        selected_chromosomes=["chr1"],
        container_image="example.invalid/rare-variant@sha256:abc",
        workflow_version="0.2.0",
        max_retries=2,
        index_provenance="supplied",
    )

    summary = json.loads(summary_path.read_text())
    assert summary["phenotype_qc"]["bed_only_sample_count"] == 2
    assert summary["chromosome_qc"]["totals"] == {
        "alt_alleles": 4,
        "boundary_variant_feature_pairs": 2,
        "info_ac_alt_alleles": 3,
    }
    provenance = summary["provenance"]
    assert provenance["selected_chromosomes"] == ["chr1"]
    assert provenance["container_image"].endswith("@sha256:abc")
    assert provenance["max_retries"] == 2
    assert provenance["vcf_index"] == "supplied"
    assert provenance["software_versions"]["workflow"] == "0.2.0"
    assert provenance["software_versions"]["rare_variant_enrichment"] == "0.2.0"
    assert "S1" not in summary_path.read_text()


def test_calculate_enrichment_emits_literal_zero_carrier_and_missing_value_table(
    tmp_path: Path,
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\n"
        "chr1\t99\t100\tGENE1\t3\t0\t-3\tNA\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\nS3\nS4\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t50\n"
        "S2\tGENE1\tAC=1\tbaseline\tall_rare_variants\t100\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary_path = tmp_path / "summary.json"
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    calculate_enrichment(
        bed,
        samples,
        carriers,
        features,
        [1, 2],
        [1],
        [2.0],
        [50, 100],
        "absolute",
        output,
        summary_path,
    )

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    literal_table = [
        (
            row["ac_class"],
            row["distance_bp"],
            row["total_observations"],
            row["outlier_observations"],
            row["nonoutlier_observations"],
            row["n11"],
            row["n10"],
            row["n01"],
            row["n00"],
            row["carrier_rate_ratio"],
            row["odds_ratio"],
            row["odds_ratio_corrected_0_5"],
        )
        for row in rows
    ]
    assert literal_table == [
        ("AC=1", "50", "3", "2", "1", "1", "1", "0", "1", "NA", "NA", "3.0"),
        ("AC=1", "100", "3", "2", "1", "1", "1", "1", "0", "0.5", "0.0", "0.3333333333333333"),
        ("AC=2", "50", "3", "2", "1", "0", "2", "0", "1", "NA", "NA", "0.6"),
        ("AC=2", "100", "3", "2", "1", "0", "2", "0", "1", "NA", "NA", "0.6"),
        ("AC<=1", "50", "3", "2", "1", "0", "2", "0", "1", "NA", "NA", "0.6"),
        ("AC<=1", "100", "3", "2", "1", "0", "2", "0", "1", "NA", "NA", "0.6"),
    ]
    summary = json.loads(summary_path.read_text())
    assert summary["missing_z_observations"] == 1
    assert summary["emitted_rows"] == 6


def test_calculate_enrichment_uses_na_when_a_rate_denominator_is_zero(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
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
        [2.0],
        [0],
        "absolute",
        tmp_path / "enrichment.tsv",
        tmp_path / "summary.json",
    )

    row = next(csv.DictReader((tmp_path / "enrichment.tsv").open(), delimiter="\t"))
    assert (row["outlier_carrier_rate"], row["nonoutlier_carrier_rate"]) == ("0.0", "NA")
    assert (row["carrier_rate_ratio"], row["odds_ratio"]) == ("NA", "NA")


@pytest.mark.parametrize(
    ("exact_ac", "cumulative_ac_max", "z_thresholds", "distance_thresholds", "message"),
    [
        ([1, 1], [], [2.0], [50], "Exact AC values must be unique"),
        ([], [1, 1], [2.0], [50], "Cumulative AC maxima must be unique"),
        ([1], [], [2.0, 2], [50], "z-score thresholds must be unique"),
        ([1], [], [2.0], [50, 50], "Distance thresholds must be unique"),
    ],
)
def test_calculate_enrichment_rejects_duplicate_configuration_values(
    tmp_path: Path,
    exact_ac: list[int],
    cumulative_ac_max: list[int],
    z_thresholds: list[float],
    distance_thresholds: list[int],
    message: str,
):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    with pytest.raises(ValueError, match=message):
        calculate_enrichment(
            bed,
            samples,
            carriers,
            features,
            exact_ac,
            cumulative_ac_max,
            z_thresholds,
            distance_thresholds,
            "absolute",
            tmp_path / "enrichment.tsv",
            tmp_path / "summary.json",
        )


def test_calculate_cli_dispatches_all_arguments(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t99\t100\tGENE1\t3\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t0\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rare_variant_enrichment.cli",
            "calculate",
            "--phenotype-bed",
            str(bed),
            "--shared-samples",
            str(samples),
            "--carriers",
            str(carriers),
            "--features",
            str(features),
            "--exact-ac",
            "1",
            "--cumulative-ac-max",
            "1",
            "--z-thresholds",
            "2.0",
            "--distance-thresholds",
            "0,100",
            "--tail",
            "absolute",
            "--output-tsv",
            str(output),
            "--output-json",
            str(summary),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert [(row["ac_class"], row["distance_bp"], row["n11"]) for row in rows] == [
        ("AC=1", "0", "1"),
        ("AC=1", "100", "1"),
        ("AC<=1", "0", "0"),
        ("AC<=1", "100", "0"),
    ]
    assert json.loads(summary.read_text())["emitted_rows"] == 4
