from pathlib import Path

import pytest

from rare_variant_enrichment.annotation_storage import VatChunkStore
from rare_variant_enrichment.annotations import GeneAnnotation, VatSchema, VariantKey


MTTOVCF_HEADER = [
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "LoF", "LoF_filter", "LoF_flags", "LoF_info", "gvs_max_af",
    "gvs_max_subpop",
]


def row(
    gene_id: str,
    consequence: str,
    gvs_max_af: str,
    lof: str,
    *,
    position: int = 100,
    transcript: str = "ENST000001",
) -> list[str]:
    return [
        "chr1", str(position), "A", "C", "rsTest", gene_id, "GENE1",
        transcript, "true", consequence, "p.Ala1Val", lof, ".", ".", ".",
        gvs_max_af, "afr",
    ]


def store(tmp_path: Path, *, maximum_gvs_maf: float = 0.01) -> VatChunkStore:
    return VatChunkStore(
        tmp_path,
        VatSchema.from_header(MTTOVCF_HEADER),
        maximum_gvs_maf,
        ["stop_gained", "missense_variant", "frameshift_variant"],
    )


def test_chunk_store_collapses_transcripts_and_uses_hc_over_lc(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG000001.1", "missense_variant", "0.001", "LC"))
        chunk_store.ingest(
            row(
                "ENSG000001.2",
                "stop_gained",
                "0.001",
                "HC",
                transcript="ENST000002",
            )
        )

        qc = chunk_store.finalize()
        key = VariantKey("chr1", 100, "A", "C")

        assert chunk_store.qualifying_maf(key) == 0.001
        assert chunk_store.gene_annotation(key, "ENSG000001") == GeneAnnotation(
            "stop_gained", "HC"
        )
        assert qc["vat_rows"] == 2
        assert qc["unique_vat_alleles"] == 1
        assert qc["unique_vat_allele_gene_pairs"] == 1


def test_chunk_store_counts_only_exact_duplicate_transcript_rows(tmp_path: Path):
    duplicate = row("ENSG000001.1", "missense_variant", "0.001", "LC")
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(duplicate)
        chunk_store.ingest(duplicate)
        chunk_store.ingest(
            row(
                "ENSG000001.1",
                "missense_variant",
                "0.001",
                "LC",
                transcript="ENST000002",
            )
        )

        qc = chunk_store.finalize()

        assert qc["vat_rows"] == 3
        assert qc["duplicate_vat_rows"] == 1
        serialized_qc = repr(qc)
        assert "ENST000001" not in serialized_qc
        assert "ENSG000001" not in serialized_qc
        assert "rsTest" not in serialized_qc


def test_chunk_store_excludes_inconsistent_and_common_frequency(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG1", "frameshift_variant", "0.001", "HC", position=100))
        chunk_store.ingest(
            row(
                "ENSG1",
                "frameshift_variant",
                "0.002",
                "HC",
                position=100,
                transcript="ENST000002",
            )
        )
        chunk_store.ingest(row("ENSG2", "frameshift_variant", "0.02", "HC", position=200))

        qc = chunk_store.finalize()

        assert chunk_store.qualifying_maf(VariantKey("chr1", 100, "A", "C")) is None
        assert chunk_store.qualifying_maf(VariantKey("chr1", 200, "A", "C")) is None
        assert qc["inconsistent_frequency_alleles"] == 1
        assert qc["above_maf_threshold_alleles"] == 1


def test_chunk_store_excludes_complementary_raw_frequency_values(tmp_path: Path):
    """Comparing normalized MAF would incorrectly keep complementary raw AFs."""
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG1", "missense_variant", "0.001", "HC"))
        chunk_store.ingest(
            row(
                "ENSG1",
                "missense_variant",
                "0.999",
                "HC",
                transcript="ENST000002",
            )
        )

        qc = chunk_store.finalize()

        assert chunk_store.qualifying_maf(VariantKey("chr1", 100, "A", "C")) is None
        assert qc["inconsistent_frequency_alleles"] == 1


def test_chunk_store_reports_terminal_frequency_exclusions_per_allele(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG1", "missense_variant", ".", "HC", position=101))
        chunk_store.ingest(row("ENSG1", "missense_variant", "not-a-number", "HC", position=102))
        chunk_store.ingest(row("ENSG1", "missense_variant", "1.1", "HC", position=103))

        qc = chunk_store.finalize()

        assert all(
            chunk_store.qualifying_maf(VariantKey("chr1", position, "A", "C")) is None
            for position in (101, 102, 103)
        )
        assert qc["missing_frequency_alleles"] == 1
        assert qc["non_numeric_frequency_alleles"] == 1
        assert qc["out_of_range_frequency_alleles"] == 1


def test_chunk_store_converts_high_af_and_retains_inclusive_threshold(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG1", "missense_variant", "0.99", "HC", position=101))
        chunk_store.ingest(row("ENSG1", "missense_variant", "0.01", "HC", position=102))

        qc = chunk_store.finalize()

        assert chunk_store.qualifying_maf(VariantKey("chr1", 101, "A", "C")) == 0.01
        assert chunk_store.qualifying_maf(VariantKey("chr1", 102, "A", "C")) == 0.01
        assert qc["observed_raw_gvs_max_af"] == 0.99
        assert qc["converted_gvs_max_af_values"] == 1


def test_chunk_store_returns_no_annotation_for_a_different_gene(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row("ENSG000001.4", "missense_variant", "0.001", "HC"))
        chunk_store.finalize()

        assert chunk_store.gene_annotation(
            VariantKey("chr1", 100, "A", "C"), "ENSG000002"
        ) == GeneAnnotation(None, None)


def test_chunk_store_records_unknown_terms_without_exposing_identifiers_in_qc(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(
            row("ENSG1", "future_variant&missense_variant&future_variant", "0.001", "HC")
        )

        qc = chunk_store.finalize()

        assert chunk_store.gene_annotation(
            VariantKey("chr1", 100, "A", "C"), "ENSG1"
        ) == GeneAnnotation("missense_variant", "HC")
        assert qc["unknown_consequence_terms"] == 1
        assert qc["unknown_consequence_rows"] == 1
        assert all("future_variant" not in key for key in qc)
        assert all("ENSG1" not in key for key in qc)


def test_chunk_store_keeps_intergenic_frequency_without_a_gene_annotation(tmp_path: Path):
    with store(tmp_path) as chunk_store:
        chunk_store.ingest(row(".", "intergenic_variant", "0.001", "."))

        qc = chunk_store.finalize()

        assert chunk_store.qualifying_maf(VariantKey("chr1", 100, "A", "C")) == 0.001
        assert chunk_store.gene_annotation(
            VariantKey("chr1", 100, "A", "C"), "ENSG000001"
        ) == GeneAnnotation(None, None)
        assert qc["unique_vat_allele_gene_pairs"] == 0


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        (row("ENSG1", "missense_variant", "0.001", "HC")[:4], "VAT row"),
        (
            row("ENSG1", "missense_variant", "0.001", "HC", position=100)[:1]
            + ["not-a-position"]
            + row("ENSG1", "missense_variant", "0.001", "HC")[2:],
            "VAT position",
        ),
        (
            [""] + row("ENSG1", "missense_variant", "0.001", "HC")[1:],
            "VAT chromosome",
        ),
    ],
)
def test_chunk_store_rejects_malformed_required_coordinates(
    tmp_path: Path, fields: list[str], message: str
):
    with store(tmp_path) as chunk_store:
        with pytest.raises(ValueError, match=message):
            chunk_store.ingest(fields)


def test_chunk_store_rejects_lookups_before_finalize_and_ingestion_after_finalize(tmp_path: Path):
    key = VariantKey("chr1", 100, "A", "C")
    with store(tmp_path) as chunk_store:
        with pytest.raises(RuntimeError, match="finalize"):
            chunk_store.qualifying_maf(key)
        with pytest.raises(RuntimeError, match="finalize"):
            chunk_store.gene_annotation(key, "ENSG1")

        chunk_store.ingest(row("ENSG1", "missense_variant", "0.001", "HC"))
        chunk_store.finalize()

        with pytest.raises(RuntimeError, match="finalize"):
            chunk_store.ingest(row("ENSG1", "missense_variant", "0.001", "HC"))


def test_chunk_store_deletes_its_database_when_ingestion_fails(tmp_path: Path):
    chunk_store = store(tmp_path)
    database_path = chunk_store.path

    with pytest.raises(ValueError, match="VAT position"):
        with chunk_store:
            malformed = row("ENSG1", "missense_variant", "0.001", "HC")
            malformed[1] = "not-a-position"
            chunk_store.ingest(malformed)

    assert not database_path.exists()
