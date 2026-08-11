from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from rare_variant_enrichment.annotations import (
    AnnotationClass,
    GeneAnnotation,
    VatSchema,
    VariantKey,
    build_annotation_classes,
    collapse_loftee,
    most_severe_consequence,
    normalize_gene_id,
    parse_consequence_terms,
    parse_gvs_max_af,
)


MTTOVCF_HEADER = [
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "LoF", "LoF_filter", "LoF_flags", "LoF_info", "gvs_max_af",
    "gvs_max_subpop",
]


def test_vat_schema_accepts_exact_mttovcf_transcript_header():
    schema = VatSchema.from_header(MTTOVCF_HEADER)

    assert (schema.chromosome, schema.position, schema.ref, schema.alt) == (0, 1, 2, 3)
    assert (schema.gene_id, schema.consequence, schema.lof, schema.gvs_max_af) == (
        5, 9, 11, 15,
    )


def test_vat_schema_allows_absent_lof_without_changing_required_indices():
    header = [column for column in MTTOVCF_HEADER if column != "LoF"]

    schema = VatSchema.from_header(header)

    assert schema.lof is None
    assert (schema.chromosome, schema.position, schema.ref, schema.alt) == (0, 1, 2, 3)
    assert (schema.gene_id, schema.consequence, schema.gvs_max_af) == (5, 9, 14)


def test_vat_schema_reports_all_missing_required_columns_in_schema_order():
    with pytest.raises(
        ValueError,
        match=r"^Missing required VAT columns: chrom, ref, consequence$",
    ):
        VatSchema.from_header(["pos", "alt", "gene_id", "gvs_max_af"])


def test_vat_schema_reports_all_duplicate_required_columns_in_schema_order():
    header = [*MTTOVCF_HEADER, "alt", "chrom", "chrom"]

    with pytest.raises(
        ValueError,
        match=r"^Duplicate required VAT columns: chrom, alt$",
    ):
        VatSchema.from_header(header)


def test_vat_schema_round_trips_and_rejects_duplicate_json_keys(tmp_path: Path):
    schema = VatSchema.from_header(MTTOVCF_HEADER)
    path = tmp_path / "schema.json"

    schema.write_json(path)

    assert VatSchema.read_json(path) == schema
    path.write_text('{"header": [], "header": []}', encoding="utf-8")
    with pytest.raises(ValueError, match=r"^VAT schema JSON contains duplicate key: header$"):
        VatSchema.read_json(path)


def test_schema_and_annotation_value_types_are_immutable():
    schema = VatSchema.from_header(MTTOVCF_HEADER)
    variant = VariantKey("chr1", 10, "A", "C")
    annotation = GeneAnnotation("missense_variant", "HC")

    with pytest.raises(FrozenInstanceError):
        schema.position = 2  # type: ignore[misc]
    assert variant == VariantKey("chr1", 10, "A", "C")
    assert annotation == GeneAnnotation("missense_variant", "HC")


def test_gene_ids_strip_only_terminal_numeric_ensembl_versions():
    assert normalize_gene_id("ENSG00000123456.17") == "ENSG00000123456"
    assert normalize_gene_id("ENSG00000123456") == "ENSG00000123456"
    assert normalize_gene_id("GENE.A") == "GENE.A"
    assert normalize_gene_id("GENE.17") == "GENE.17"


def test_gene_id_normalization_preserves_surrounding_whitespace():
    padded = " ENSG00000123456.17 "

    assert normalize_gene_id(padded) == padded
    assert normalize_gene_id(padded) != "ENSG00000123456"


def test_empty_normalized_gene_id_reports_the_input():
    with pytest.raises(ValueError, match=r"^Gene ID is empty after normalization: ''$"):
        normalize_gene_id("")


def test_delimited_consequences_collapse_by_ensembl_severity():
    terms = (
        *parse_consequence_terms("intron_variant,splice_region_variant"),
        *parse_consequence_terms("missense_variant&stop_gained"),
    )

    selected, unknown = most_severe_consequence(terms)

    assert selected == "stop_gained"
    assert unknown == ()


def test_consequence_parser_deduplicates_and_preserves_order():
    assert parse_consequence_terms(" missense_variant & stop_gained,missense_variant ") == (
        "missense_variant", "stop_gained",
    )
    assert parse_consequence_terms("NA") == ()
    assert parse_consequence_terms(".") == ()


def test_most_severe_consequence_returns_unknown_terms_for_qc():
    selected, unknown = most_severe_consequence(
        ["future_variant", "missense_variant", "future_variant", ""]
    )

    assert selected == "missense_variant"
    assert unknown == ("future_variant",)


def test_loftee_and_frequency_normalization_are_deterministic():
    assert collapse_loftee(["LC", "hc", "."]) == "HC"
    assert collapse_loftee([".", "low_confidence"]) is None
    converted = parse_gvs_max_af("0.999")
    assert (converted.status, converted.maf, converted.converted) == ("valid", 0.001, True)
    assert parse_gvs_max_af("NA").status == "missing"
    assert parse_gvs_max_af("not-a-number").status == "non_numeric"
    assert parse_gvs_max_af("1.1").status == "out_of_range"
    assert parse_gvs_max_af("0.4").converted is False


def test_boolean_and_non_finite_af_inputs_report_specific_categories():
    with pytest.raises(ValueError, match=r"^gvs_max_af must not be a Boolean: True$"):
        parse_gvs_max_af(True)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match=r"^gvs_max_af must be finite: 'NaN'$"):
        parse_gvs_max_af("NaN")
    with pytest.raises(ValueError, match=r"^gvs_max_af must be finite: inf$"):
        parse_gvs_max_af(float("inf"))  # type: ignore[arg-type]


def test_annotation_classes_have_stable_family_order():
    classes = build_annotation_classes(["frameshift_variant", "stop_gained"], True)

    assert classes == [
        AnnotationClass("baseline", "all_rare_variants"),
        AnnotationClass("consequence", "frameshift_variant"),
        AnnotationClass("consequence", "stop_gained"),
        AnnotationClass("loftee", "HC"),
        AnnotationClass("loftee", "LC"),
    ]


def test_annotation_classes_report_sorted_duplicate_consequence_names():
    with pytest.raises(
        ValueError,
        match=r"^Duplicate consequence classes: frameshift_variant, stop_gained$",
    ):
        build_annotation_classes(
            ["stop_gained", "frameshift_variant", "stop_gained", "frameshift_variant"],
            False,
        )


def test_annotation_classes_can_contain_baseline_only():
    assert build_annotation_classes([], False) == [
        AnnotationClass("baseline", "all_rare_variants")
    ]
