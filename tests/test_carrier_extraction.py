import csv
import gzip
import json
from pathlib import Path

import pytest

from rare_variant_enrichment.carrier_extraction import (
    build_chromosome_chunks,
    extract_chromosome_carriers,
    prepare_carrier_inputs,
)
from rare_variant_enrichment.variants import QueryChunk


def test_build_chromosome_chunks_covers_contig_without_overlap():
    assert build_chromosome_chunks("chr1", 1000, 400) == (
        QueryChunk("chr1", 1, 400),
        QueryChunk("chr1", 401, 800),
        QueryChunk("chr1", 801, 1000),
    )


@pytest.mark.parametrize("chunk_size", [True, 0, -1])
def test_build_chromosome_chunks_rejects_invalid_size(chunk_size: int):
    with pytest.raises(ValueError, match="chunk_size_bp must be a positive integer"):
        build_chromosome_chunks("chr1", 1000, chunk_size)


def test_prepare_records_schema_contigs_and_index_provenance(
    prepared_fixture, tmp_path: Path
):
    schema = tmp_path / "schema.json"
    qc = tmp_path / "prepare.qc.json"

    prepare_carrier_inputs(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        ["chr1"],
        "supplied",
        "generated",
        schema,
        qc,
    )

    payload = json.loads(qc.read_text())
    assert payload["selected_chromosomes"] == ["chr1"]
    assert payload["vcf_contig_lengths"] == {"chr1": 1000}
    assert payload["sample_count"] == 3
    assert payload["revel_available"] is True
    assert payload["vcf_index_provenance"] == "supplied"
    assert payload["transcript_index_provenance"] == "generated"
    assert schema.exists()
    assert "S1" not in qc.read_text()


def test_extract_emits_exact_gene_matches_dosage_and_unclassified_rows(
    prepared_fixture, tmp_path: Path
):
    audit = tmp_path / "chr1.audit.tsv.gz"
    qc = tmp_path / "chr1.qc.json"

    extract_chromosome_carriers(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        prepared_fixture.carrier_schema,
        "chr1",
        250,
        audit,
        qc,
    )

    with gzip.open(audit, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert "distance" not in reader.fieldnames
    at_100 = [row for row in rows if row["variant_id"] == "chr1:100:A:C"]
    assert {row["gene_id"] for row in at_100} == {"ENSG000001", "ENSG000099"}
    gene1 = next(row for row in at_100 if row["gene_id"] == "ENSG000001")
    assert gene1["most_severe_consequence"] == "stop_gained"
    assert gene1["loftee"] == "HC"
    assert gene1["revel"] == "0.81"
    assert gene1["variant_classes"] == "lof_hc,lof_hc_or_lc"

    missense = next(row for row in rows if row["variant_id"] == "chr1:150:G:T")
    assert missense["sample_id"] == "S2"
    assert missense["sample_alt_allele_count"] == "2"
    assert missense["variant_classes"] == "missense"
    assert any(
        row["most_severe_consequence"] == "synonymous_variant"
        and row["variant_classes"] == ""
        for row in rows
    )

    payload = json.loads(qc.read_text())
    assert payload["quality_or_frequency_filters_applied"] is False
    assert payload["carrier_audit_rows"] == len(rows)
    assert "sample_ids" not in payload
