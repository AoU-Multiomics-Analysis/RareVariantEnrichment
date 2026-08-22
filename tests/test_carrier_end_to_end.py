import csv
import gzip
import json
from pathlib import Path

from rare_variant_enrichment.carrier_aggregation import gather_variant_carriers
from rare_variant_enrichment.carrier_extraction import (
    extract_chromosome_carriers,
    prepare_carrier_inputs,
)


def test_gene_matched_carrier_extraction_end_to_end(prepared_fixture, tmp_path: Path):
    schema = tmp_path / "schema.json"
    preparation_qc = tmp_path / "prepare.qc.json"
    chromosome_audit = tmp_path / "chr1.audit.tsv.gz"
    chromosome_qc = tmp_path / "chr1.qc.json"
    final_audit = tmp_path / "variant_carrier_audit.tsv.gz"
    carriers = tmp_path / "variant_carriers.tsv.gz"
    gathered_qc = tmp_path / "variant_carriers.qc.json"

    prepare_carrier_inputs(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        ["chr1"],
        "supplied",
        "supplied",
        schema,
        preparation_qc,
    )
    extract_chromosome_carriers(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        schema,
        "chr1",
        250,
        chromosome_audit,
        chromosome_qc,
    )
    gather_variant_carriers(
        [chromosome_audit],
        [chromosome_qc],
        preparation_qc,
        final_audit,
        carriers,
        gathered_qc,
    )

    with gzip.open(final_audit, "rt", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(audit_rows) == 5
    gene1 = next(row for row in audit_rows if row["gene_id"] == "ENSG000001")
    assert gene1["most_severe_consequence"] == "stop_gained"
    assert gene1["revel"] == "0.81"
    assert {row["gene_id"] for row in audit_rows if row["variant_id"] == "chr1:100:A:C"} == {
        "ENSG000001", "ENSG000099"
    }
    missense = next(row for row in audit_rows if row["variant_id"] == "chr1:150:G:T")
    assert missense["sample_alt_allele_count"] == "2"
    assert missense["variant_classes"] == "missense"
    assert any(row["variant_classes"] == "" for row in audit_rows)

    with gzip.open(carriers, "rt", newline="") as handle:
        carrier_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(
        row["sample_id"] == "S1"
        and row["gene_id"] == "ENSG000001"
        and row["variant_class"] == "lof_hc"
        and row["n_variants"] == "1"
        and row["variant_ids"] == "chr1:100:A:C"
        for row in carrier_rows
    )

    chromosome_payload = json.loads(chromosome_qc.read_text())
    gathered_payload = json.loads(gathered_qc.read_text())
    assert chromosome_payload["carrier_audit_rows"] == len(audit_rows)
    assert gathered_payload["audit_row_count"] == len(audit_rows)
    assert gathered_payload["carrier_row_count"] == len(carrier_rows)
