from pathlib import Path

import pytest

from rare_variant_enrichment.annotations import VariantKey
from rare_variant_enrichment.carrier_annotations import (
    TranscriptCarrierRow,
    TranscriptCarrierSchema,
    collapse_transcript_rows,
    initial_variant_classes,
    parse_optional_unit_interval,
)


HEADER = (
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "revel", "LoF", "LoF_filter", "LoF_flags", "LoF_info",
    "gvs_max_af", "gvs_max_subpop",
)


def _row(
    gene_id: str,
    symbol: str | None,
    consequence: str,
    loftee: str | None,
    revel: float | None,
    gvs_max_af: float | None,
    *,
    key: VariantKey = VariantKey("chr1", 100, "A", "C"),
) -> TranscriptCarrierRow:
    return TranscriptCarrierRow(
        key, gene_id, symbol, (consequence,), loftee, revel, gvs_max_af
    )


def test_transcript_carrier_schema_round_trips(tmp_path: Path):
    schema = TranscriptCarrierSchema.from_header(HEADER)
    path = tmp_path / "schema.json"
    schema.write_json(path)

    assert TranscriptCarrierSchema.read_json(path) == schema
    assert schema.gene_id == HEADER.index("gene_id")
    assert schema.revel == HEADER.index("revel")


@pytest.mark.parametrize("value", ["", ".", "NA", "null"])
def test_optional_unit_interval_accepts_missing(value: str):
    assert parse_optional_unit_interval(value, "REVEL") is None


@pytest.mark.parametrize("value, expected", [("0", 0.0), ("0.81", 0.81), ("1", 1.0)])
def test_optional_unit_interval_accepts_bounds(value: str, expected: float):
    assert parse_optional_unit_interval(value, "REVEL") == expected


@pytest.mark.parametrize("value", ["NaN", "inf", "-0.1", "1.1", "bad"])
def test_optional_unit_interval_rejects_invalid_values(value: str):
    with pytest.raises(ValueError, match="REVEL"):
        parse_optional_unit_interval(value, "REVEL")


def test_schema_rejects_missing_and_duplicate_required_columns():
    with pytest.raises(ValueError, match="Missing.*revel"):
        TranscriptCarrierSchema.from_header(tuple(value for value in HEADER if value != "revel"))
    with pytest.raises(ValueError, match="Duplicate.*gene_id"):
        TranscriptCarrierSchema.from_header((*HEADER, "gene_id"))


def test_collapse_selects_most_severe_hc_and_maximum_scores():
    collapsed = collapse_transcript_rows([
        _row("ENSG1.1", "GENE1", "missense_variant", "LC", 0.42, 0.001),
        _row("ENSG1.2", "GENE1", "splice_acceptor_variant", "HC", 0.81, 0.004),
        _row("ENSG1.3", "GENE1", "synonymous_variant", None, None, None),
    ])

    assert collapsed.gene_id == "ENSG1"
    assert collapsed.most_severe_consequence == "splice_acceptor_variant"
    assert collapsed.all_consequences == (
        "splice_acceptor_variant", "missense_variant", "synonymous_variant"
    )
    assert collapsed.loftee == "HC"
    assert collapsed.revel == 0.81
    assert collapsed.gvs_max_af == 0.004
    assert initial_variant_classes(collapsed) == (
        "lof_hc", "lof_hc_or_lc", "splice_core"
    )


@pytest.mark.parametrize(
    "consequence, expected",
    [
        ("missense_variant", ("missense",)),
        ("splice_donor_variant", ("splice_core",)),
        ("splice_region_variant", ("splice_region",)),
        ("synonymous_variant", ()),
    ],
)
def test_initial_classes_use_the_collapsed_most_severe_consequence(
    consequence: str, expected: tuple[str, ...]
):
    collapsed = collapse_transcript_rows([
        _row("ENSG1", "GENE1", consequence, None, None, None)
    ])
    assert initial_variant_classes(collapsed) == expected


def test_collapse_rejects_mixed_alleles_genes_and_symbols():
    base = _row("ENSG1", "GENE1", "missense_variant", None, None, None)
    with pytest.raises(ValueError, match="same exact allele"):
        collapse_transcript_rows([
            base,
            _row(
                "ENSG1", "GENE1", "missense_variant", None, None, None,
                key=VariantKey("chr1", 101, "A", "C"),
            ),
        ])
    with pytest.raises(ValueError, match="same normalized gene"):
        collapse_transcript_rows([
            base, _row("ENSG2", "GENE2", "missense_variant", None, None, None)
        ])
    with pytest.raises(ValueError, match="gene symbols"):
        collapse_transcript_rows([
            base, _row("ENSG1", "OTHER", "missense_variant", None, None, None)
        ])
