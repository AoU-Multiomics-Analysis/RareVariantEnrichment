import pytest

from rare_variant_enrichment.variants import (
    FeatureTssIndex,
    FeatureTss,
    QueryChunk,
    build_ac_classes,
    build_query_chunks,
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


def test_build_ac_classes_rejects_no_configured_family():
    with pytest.raises(ValueError, match="At least one AC class is required"):
        build_ac_classes([], [])


def test_multiallelic_record_assigns_ac_and_carriers_per_alt():
    fields = "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=1,2\tGT\t0/1\t0/2\t0/2".split("\t")

    alleles = parse_variant_alleles(fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"})

    assert [(allele.alt, allele.ac, allele.carriers) for allele in alleles] == [
        ("C", 1, ("S1",)),
        ("G", 2, ("S2", "S3")),
    ]


def test_info_ac_dot_falls_back_to_genotypes_and_records_source_qc():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\tAC=.\tGT\t0/1\t0/0".split("\t")
    qc: dict[str, int] = {}

    alleles = parse_variant_alleles(
        fields, ["S1", "S2"], {"S1", "S2"}, qc=qc
    )

    assert [(allele.ac, allele.carriers) for allele in alleles] == [(1, ("S1",))]
    assert qc["info_ac_alt_alleles"] == 0
    assert qc["genotype_ac_fallback_alt_alleles"] == 1
    assert qc["unavailable_ac_alt_alleles"] == 0


def test_multiallelic_info_ac_falls_back_per_alt():
    fields = "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=.,2\tGT\t0/1\t0/2\t0/2".split("\t")
    qc: dict[str, int] = {}

    alleles = parse_variant_alleles(
        fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"}, qc=qc
    )

    assert [(allele.alt, allele.ac, allele.carriers) for allele in alleles] == [
        ("C", 1, ("S1",)),
        ("G", 2, ("S2", "S3")),
    ]
    assert qc["info_ac_alt_alleles"] == 1
    assert qc["genotype_ac_fallback_alt_alleles"] == 1


def test_negative_info_ac_is_rejected():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\tAC=-1\tGT\t0/1".split("\t")

    with pytest.raises(ValueError, match="INFO/AC values must be non-negative"):
        parse_variant_alleles(fields, ["S1"], {"S1"})


def test_info_ac_is_authoritative_and_complete_genotype_disagreement_is_counted():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\tAC=2\tGT\t0/1\t0/0".split("\t")
    qc: dict[str, int] = {}

    allele = parse_variant_alleles(
        fields, ["S1", "S2"], {"S1", "S2"}, qc=qc
    )[0]

    assert allele.ac == 2
    assert allele.carriers == ("S1",)
    assert qc["info_genotype_ac_compared_alt_alleles"] == 1
    assert qc["info_genotype_ac_mismatch_alt_alleles"] == 1


def test_partial_calls_count_known_alt_alleles_but_fully_missing_calls_do_not():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t1/.\t./.\t0/0".split("\t")
    qc: dict[str, int] = {}

    allele = parse_variant_alleles(
        fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"}, qc=qc
    )[0]

    assert allele.ac == 1
    assert allele.carriers == ("S1",)
    assert qc["partial_genotype_calls"] == 1
    assert qc["fully_missing_genotype_calls"] == 1
    assert qc["missing_genotypes"] == 2


def test_fully_missing_calls_leave_missing_info_ac_unavailable():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t./.\t.".split("\t")
    qc: dict[str, int] = {}

    alleles = parse_variant_alleles(
        fields, ["S1", "S2"], {"S1", "S2"}, qc=qc
    )

    assert alleles == []
    assert qc["fully_missing_genotype_calls"] == 2
    assert qc["unavailable_ac_alt_alleles"] == 1


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


def test_feature_tss_index_returns_only_bisected_interval_candidates():
    features = [
        FeatureTss("chr1", 10_000 + offset * 100, f"G{offset}")
        for offset in reversed(range(200))
    ]
    index = FeatureTssIndex(features)

    matches = index.nearby(15_050, 50)

    assert [feature.feature_id for feature in matches] == ["G50", "G51"]


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


def test_query_chunks_split_merged_windows_without_overlap():
    features = [
        FeatureTss("chr1", 10, "ENSG1.1"),
        FeatureTss("chr1", 25, "ENSG2.2"),
    ]

    chunks = build_query_chunks(features, max_distance=10, chunk_size_bp=10)

    assert chunks == [
        QueryChunk("chr1", 1, 10),
        QueryChunk("chr1", 11, 20),
        QueryChunk("chr1", 21, 30),
        QueryChunk("chr1", 31, 35),
    ]
    assert [chunk.tabix_region for chunk in chunks] == [
        "chr1:1-10",
        "chr1:11-20",
        "chr1:21-30",
        "chr1:31-35",
    ]
    assert all(left.end + 1 == right.start for left, right in zip(chunks, chunks[1:]))


@pytest.mark.parametrize("chunk_size_bp", [0, -1, True])
def test_query_chunks_reject_invalid_chunk_sizes(chunk_size_bp: int):
    with pytest.raises(ValueError, match=r"^annotation_chunk_size_bp must be a positive integer$"):
        build_query_chunks([FeatureTss("chr1", 10, "ENSG1")], 10, chunk_size_bp)


def test_query_chunks_reject_features_from_multiple_chromosomes():
    features = [FeatureTss("chr1", 10, "ENSG1"), FeatureTss("chr2", 20, "ENSG2")]

    with pytest.raises(ValueError, match=r"^Query chunks require features from one chromosome$"):
        build_query_chunks(features, max_distance=10, chunk_size_bp=10)
