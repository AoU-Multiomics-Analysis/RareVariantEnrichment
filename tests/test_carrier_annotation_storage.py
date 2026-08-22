from pathlib import Path

import pytest

from rare_variant_enrichment.annotations import VariantKey
from rare_variant_enrichment.carrier_annotation_storage import CarrierAnnotationChunkStore
from rare_variant_enrichment.carrier_annotations import TranscriptCarrierSchema


HEADER = (
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "revel", "LoF", "LoF_filter", "LoF_flags", "LoF_info",
    "gvs_max_af", "gvs_max_subpop",
)


def _fields(gene: str, consequence: str, revel: str, loftee: str = ".") -> list[str]:
    values = {name: "." for name in HEADER}
    values.update({
        "chrom": "chr1", "pos": "100", "ref": "A", "alt": "C",
        "gene_id": gene, "gene_symbol": "GENE1", "consequence": consequence,
        "revel": revel, "LoF": loftee, "gvs_max_af": "0.004",
    })
    return [values[name] for name in HEADER]


def test_chunk_store_collapses_each_exact_allele_gene_pair(tmp_path: Path):
    schema = TranscriptCarrierSchema.from_header(HEADER)
    with CarrierAnnotationChunkStore(tmp_path, schema) as store:
        first = _fields("ENSG1.1", "missense_variant", "0.42", "LC")
        second = _fields("ENSG1.2", "splice_acceptor_variant", "0.81", "HC")
        store.ingest(first)
        store.ingest(second)
        store.ingest(first)
        qc = store.finalize()
        annotations = store.annotations_for_allele(VariantKey("chr1", 100, "A", "C"))

    assert len(annotations) == 1
    assert annotations[0].gene_id == "ENSG1"
    assert annotations[0].most_severe_consequence == "splice_acceptor_variant"
    assert annotations[0].revel == 0.81
    assert qc == {
        "transcript_rows": 3,
        "duplicate_transcript_rows": 1,
        "unique_annotation_alleles": 1,
        "unique_annotation_allele_gene_pairs": 1,
    }


def test_chunk_store_requires_finalize_before_query(tmp_path: Path):
    schema = TranscriptCarrierSchema.from_header(HEADER)
    with CarrierAnnotationChunkStore(tmp_path, schema) as store:
        with pytest.raises(RuntimeError, match="finalized"):
            store.annotations_for_allele(VariantKey("chr1", 100, "A", "C"))


def test_chunk_store_removes_its_temporary_database(tmp_path: Path):
    schema = TranscriptCarrierSchema.from_header(HEADER)
    store = CarrierAnnotationChunkStore(tmp_path, schema)
    database = store.path
    assert database.exists()
    store.close()
    assert not database.exists()
