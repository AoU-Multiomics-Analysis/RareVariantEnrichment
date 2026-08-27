import csv
import gzip
import hashlib
import json
from pathlib import Path

from rare_variant_enrichment.carrier_aggregation import gather_variant_carriers
from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER


def _write_audit(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _row(variant_id: str, classes: str, *, sample: str = "S1") -> dict[str, str]:
    chrom, pos, ref, alt = variant_id.split(":")
    return {
        "sample_id": sample,
        "gene_id": "ENSG1",
        "gene_symbol": "GENE1",
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_id": variant_id,
        "variant_ac": "1",
        "variant_af": "0.001",
        "sample_alt_allele_count": "1",
        "most_severe_consequence": "stop_gained",
        "all_consequences": "stop_gained",
        "unknown_consequences": "",
        "loftee": "HC",
        "revel": "0.81",
        "gvs_max_af": "0.001",
        "variant_classes": classes,
    }


def test_gather_deduplicates_audit_and_aggregates_classes(tmp_path: Path):
    chr2_audit = tmp_path / "chr2.tsv.gz"
    chr1_audit = tmp_path / "chr1.tsv.gz"
    duplicate = _row("chr1:100:A:C", "lof_hc,lof_hc_or_lc")
    _write_audit(chr2_audit, [_row("chr2:200:G:T", "lof_hc", sample="S1")])
    _write_audit(chr1_audit, [duplicate, duplicate])
    chr2_qc = tmp_path / "chr2.json"
    chr1_qc = tmp_path / "chr1.json"
    chr2_qc.write_text(json.dumps({"chromosome": "chr2", "carrier_audit_rows": 1}))
    chr1_qc.write_text(json.dumps({"chromosome": "chr1", "carrier_audit_rows": 2}))
    preparation_qc = tmp_path / "prepare.json"
    preparation_qc.write_text(json.dumps({
        "vcf_index_provenance": "supplied",
        "transcript_index_provenance": "generated",
    }))

    audit = tmp_path / "variant_carrier_audit.tsv.gz"
    carriers = tmp_path / "variant_carriers.tsv.gz"
    qc = tmp_path / "variant_carriers.qc.json"
    gather_variant_carriers(
        [chr2_audit, chr1_audit],
        [chr2_qc, chr1_qc],
        preparation_qc,
        audit,
        carriers,
        qc,
    )

    with gzip.open(audit, "rt", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["variant_id"] for row in audit_rows] == ["chr1:100:A:C", "chr2:200:G:T"]

    with gzip.open(carriers, "rt", newline="") as handle:
        carrier_rows = list(csv.DictReader(handle, delimiter="\t"))
    by_class = {row["variant_class"]: row for row in carrier_rows}
    assert by_class["lof_hc"]["n_variants"] == "2"
    assert by_class["lof_hc"]["variant_ids"] == "chr1:100:A:C,chr2:200:G:T"
    assert by_class["lof_hc_or_lc"]["n_variants"] == "1"

    payload = json.loads(qc.read_text())
    assert payload["duplicate_audit_rows"] == 1
    assert payload["audit_row_count"] == 2
    assert payload["carrier_row_count"] == 2
    assert payload["audit_artifact"] == {
        "logical_name": "variant_carrier_audit.tsv.gz",
        "header": list(AUDIT_HEADER),
        "row_count": 2,
        "size_bytes": audit.stat().st_size,
        "sha256": hashlib.sha256(audit.read_bytes()).hexdigest(),
    }
    assert payload["quality_or_frequency_filters_applied"] is False
    assert payload["vcf_index_provenance"] == "supplied"
    assert payload["transcript_index_provenance"] == "generated"
