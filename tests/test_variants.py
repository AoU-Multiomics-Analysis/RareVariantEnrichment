import pytest

from rare_variant_enrichment.variants import (
    FeatureTss,
    build_ac_classes,
    merge_tss_windows,
    nearby_features,
    parse_variant_alleles,
)


def test_build_ac_classes_is_stable_and_deduplicated():
    classes = build_ac_classes([1, 2, 2, 3], [1, 3])

    assert [item.label for item in classes] == ["AC=1", "AC=2", "AC=3", "AC<=1", "AC<=3"]


def test_ac_classes_match_their_exact_or_cumulative_cutoff():
    exact, cumulative = build_ac_classes([2], [2])

    assert exact.contains(2)
    assert not exact.contains(1)
    assert cumulative.contains(1)
    assert cumulative.contains(2)
    assert not cumulative.contains(3)


def test_multiallelic_record_assigns_ac_and_carriers_per_alt():
    fields = "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=1,2\tGT\t0/1\t0/2\t0/2".split("\t")

    alleles = parse_variant_alleles(fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"})

    assert [(allele.alt, allele.ac, allele.carriers) for allele in alleles] == [
        ("C", 1, ("S1",)),
        ("G", 2, ("S2", "S3")),
    ]


def test_genotypes_supply_global_ac_but_only_shared_carriers():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1\t./.\t0/1".split("\t")

    allele = parse_variant_alleles(fields, ["S1", "S2", "S3", "S4"], {"S1", "S2", "S3"})[0]

    assert allele.ac == 4
    assert allele.carriers == ("S1", "S2")


def test_genotypes_count_phased_and_unphased_alt_indices_per_alt():
    fields = "chr2\t200\t.\tT\tC,G\t.\tPASS\t.\tGT\t1|2\t2/2\t0/1".split("\t")

    alleles = parse_variant_alleles(fields, ["S1", "S2", "S3"], {"S1", "S3"})

    assert [(allele.ac, allele.carriers) for allele in alleles] == [
        (2, ("S1", "S3")),
        (3, ("S1",)),
    ]


def test_parse_variant_alleles_rejects_mismatched_info_ac_cardinality():
    fields = "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=1\tGT\t0/1".split("\t")

    with pytest.raises(ValueError, match="INFO/AC.*ALT"):
        parse_variant_alleles(fields, ["S1"], {"S1"})


def test_windows_merge_and_distance_boundary_is_inclusive():
    features = [FeatureTss("chr1", 100, "G1"), FeatureTss("chr1", 120, "G2")]

    assert merge_tss_windows(features, 10) == [("chr1", 89, 130)]
    assert [feature.feature_id for feature in nearby_features(features, 110, 10)] == ["G1", "G2"]


def test_tss_windows_preserve_nonoverlapping_chromosomes_and_bed_boundaries():
    features = [
        FeatureTss("chr1", 10, "G1"),
        FeatureTss("chr1", 31, "G2"),
        FeatureTss("chr2", 1, "G3"),
    ]

    assert merge_tss_windows(features, 10) == [
        ("chr1", 0, 20),
        ("chr1", 20, 41),
        ("chr2", 0, 11),
    ]
    assert nearby_features(features, 21, 10) == [FeatureTss("chr1", 31, "G2")]
