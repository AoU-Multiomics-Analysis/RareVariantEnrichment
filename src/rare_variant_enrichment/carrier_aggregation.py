"""Gather chromosome carrier audits and build the class-level carrier table."""

from __future__ import annotations

import csv
import gzip
from io import TextIOWrapper
import json
import math
from pathlib import Path
import sqlite3
import tempfile
from typing import Any, Sequence

from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER
from rare_variant_enrichment.io import open_text, write_json


CARRIER_HEADER = (
    "sample_id",
    "gene_id",
    "gene_symbol",
    "variant_class",
    "n_variants",
    "variant_ids",
)
_VALID_CLASSES = {
    "lof_hc", "lof_hc_or_lc", "missense", "splice_core", "splice_region"
}


def gather_variant_carriers(
    audit_inputs: Sequence[Path],
    qc_inputs: Sequence[Path],
    preparation_qc_path: Path,
    audit_output: Path,
    carrier_output: Path,
    qc_output: Path,
) -> None:
    if len(audit_inputs) != len(qc_inputs):
        raise ValueError("Audit and QC input counts must match")
    preparation_qc = _read_json_object(preparation_qc_path)
    chromosome_qc = [_read_json_object(path) for path in qc_inputs]
    chromosomes = [payload.get("chromosome") for payload in chromosome_qc]
    if len(set(chromosomes)) != len(chromosomes):
        raise ValueError("Chromosome QC inputs must have unique chromosomes")

    temporary = tempfile.NamedTemporaryFile(
        prefix="carrier-gather-", suffix=".sqlite3", dir=audit_output.parent, delete=False
    )
    temporary.close()
    database_path = Path(temporary.name)
    duplicate_rows = 0
    raw_rows = 0
    symbols: dict[str, str] = {}
    try:
        connection = sqlite3.connect(database_path)
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        connection.executescript(
            """
            CREATE TABLE audit (
                sample_id TEXT NOT NULL,
                gene_id TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                variant_id TEXT NOT NULL,
                row_data TEXT NOT NULL,
                PRIMARY KEY (sample_id, gene_id, chrom, pos, ref, alt)
            ) WITHOUT ROWID;

            CREATE TABLE class_variant (
                sample_id TEXT NOT NULL,
                gene_id TEXT NOT NULL,
                gene_symbol TEXT NOT NULL,
                variant_class TEXT NOT NULL,
                variant_id TEXT NOT NULL,
                PRIMARY KEY (sample_id, gene_id, variant_class, variant_id)
            ) WITHOUT ROWID;
            """
        )
        for audit_path in audit_inputs:
            with open_text(audit_path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if tuple(reader.fieldnames or ()) != AUDIT_HEADER:
                    raise ValueError(f"Carrier audit has an invalid header: {audit_path}")
                for row in reader:
                    raw_rows += 1
                    _validate_audit_row(row)
                    symbol = row["gene_symbol"]
                    previous_symbol = symbols.get(row["gene_id"], "")
                    if symbol and previous_symbol and symbol != previous_symbol:
                        raise ValueError(f"Conflicting symbols for gene {row['gene_id']}")
                    if symbol:
                        symbols[row["gene_id"]] = symbol
                    row_data = json.dumps(
                        [row[name] for name in AUDIT_HEADER],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    )
                    key = (
                        row["sample_id"], row["gene_id"], row["chrom"], int(row["pos"]),
                        row["ref"], row["alt"],
                    )
                    existing = connection.execute(
                        """
                        SELECT row_data FROM audit
                        WHERE sample_id = ? AND gene_id = ? AND chrom = ?
                          AND pos = ? AND ref = ? AND alt = ?
                        """,
                        key,
                    ).fetchone()
                    if existing is not None:
                        if existing[0] != row_data:
                            raise ValueError("Conflicting duplicate carrier audit row")
                        duplicate_rows += 1
                        continue
                    connection.execute(
                        """
                        INSERT INTO audit (
                            sample_id, gene_id, gene_symbol, chrom, pos, ref, alt,
                            variant_id, row_data
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (*key[:2], symbol, *key[2:], row["variant_id"], row_data),
                    )
                    for variant_class in _parse_classes(row["variant_classes"]):
                        connection.execute(
                            """
                            INSERT INTO class_variant (
                                sample_id, gene_id, gene_symbol, variant_class, variant_id
                            ) VALUES (?, ?, ?, ?, ?)
                            """,
                            (
                                row["sample_id"], row["gene_id"], symbol,
                                variant_class, row["variant_id"],
                            ),
                        )
        connection.commit()
        expected_raw_rows = sum(
            int(payload.get("carrier_audit_rows", 0)) for payload in chromosome_qc
        )
        if expected_raw_rows != raw_rows:
            raise ValueError("Chromosome QC audit row counts do not match audit inputs")

        audit_count = _write_audit(connection, audit_output)
        carrier_count = _write_carriers(connection, carrier_output)
        unique_sample_count = _scalar(connection, "SELECT count(DISTINCT sample_id) FROM audit")
        unique_gene_count = _scalar(connection, "SELECT count(DISTINCT gene_id) FROM audit")
        unique_allele_count = _scalar(
            connection,
            "SELECT count(*) FROM (SELECT chrom, pos, ref, alt FROM audit GROUP BY chrom, pos, ref, alt)",
        )
        unique_pair_count = _scalar(
            connection,
            "SELECT count(*) FROM (SELECT gene_id, chrom, pos, ref, alt FROM audit GROUP BY gene_id, chrom, pos, ref, alt)",
        )
        class_counts = {
            row[0]: int(row[1])
            for row in connection.execute(
                "SELECT variant_class, count(*) FROM class_variant GROUP BY variant_class ORDER BY variant_class"
            )
        }
        write_json(
            qc_output,
            {
                "preparation_qc": preparation_qc,
                "chromosome_qc": chromosome_qc,
                "audit_inputs": [str(path) for path in audit_inputs],
                "qc_inputs": [str(path) for path in qc_inputs],
                "duplicate_audit_rows": duplicate_rows,
                "audit_row_count": audit_count,
                "carrier_row_count": carrier_count,
                "unique_sample_count": unique_sample_count,
                "unique_gene_count": unique_gene_count,
                "unique_allele_count": unique_allele_count,
                "unique_allele_gene_pair_count": unique_pair_count,
                "class_variant_counts": class_counts,
                "vcf_index_provenance": preparation_qc.get("vcf_index_provenance"),
                "transcript_index_provenance": preparation_qc.get(
                    "transcript_index_provenance"
                ),
                "quality_or_frequency_filters_applied": False,
            },
        )
    finally:
        if "connection" in locals():
            connection.close()
        database_path.unlink(missing_ok=True)


def _write_audit(connection: sqlite3.Connection, path: Path) -> int:
    count = 0
    with _open_deterministic_gzip_text(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(AUDIT_HEADER)
        rows = connection.execute(
            "SELECT row_data FROM audit ORDER BY chrom, pos, ref, alt, gene_id, sample_id"
        )
        for (row_data,) in rows:
            writer.writerow(json.loads(row_data))
            count += 1
    return count


def _write_carriers(connection: sqlite3.Connection, path: Path) -> int:
    count = 0
    with _open_deterministic_gzip_text(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(CARRIER_HEADER)
        rows = connection.execute(
            """
            SELECT sample_id, gene_id, max(gene_symbol), variant_class,
                   count(*), group_concat(variant_id, ',')
            FROM (
                SELECT * FROM class_variant
                ORDER BY sample_id, gene_id, variant_class, variant_id
            )
            GROUP BY sample_id, gene_id, variant_class
            ORDER BY sample_id, gene_id, variant_class
            """
        )
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def _validate_audit_row(row: dict[str, str]) -> None:
    for field in ("sample_id", "gene_id", "chrom", "pos", "ref", "alt", "variant_id"):
        if not row[field]:
            raise ValueError(f"Carrier audit field is empty: {field}")
    try:
        position = int(row["pos"])
        ac = int(row["variant_ac"])
        dosage = int(row["sample_alt_allele_count"])
    except ValueError as error:
        raise ValueError("Carrier audit integer field is invalid") from error
    if position < 1 or ac < 0 or dosage < 1:
        raise ValueError("Carrier audit integer field is outside its valid range")
    expected_id = f"{row['chrom']}:{position}:{row['ref']}:{row['alt']}"
    if row["variant_id"] != expected_id:
        raise ValueError("Carrier audit variant_id does not match its allele fields")
    for field in ("variant_af", "revel", "gvs_max_af"):
        if not row[field]:
            continue
        try:
            value = float(row[field])
        except ValueError as error:
            raise ValueError(f"Carrier audit {field} must be numeric or missing") from error
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(f"Carrier audit {field} must be from 0 through 1")
    _parse_classes(row["variant_classes"])


def _parse_classes(value: str) -> tuple[str, ...]:
    classes = tuple(item for item in value.split(",") if item)
    if len(classes) != len(set(classes)):
        raise ValueError("Carrier audit has duplicate variant classes")
    invalid = sorted(set(classes) - _VALID_CLASSES)
    if invalid:
        raise ValueError("Carrier audit has invalid variant classes: " + ", ".join(invalid))
    return classes


def _read_json_object(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except json.JSONDecodeError as error:
        raise ValueError(f"QC JSON is invalid: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"QC JSON must contain an object: {path}")
    return payload


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"QC JSON contains duplicate key: {key}")
        payload[key] = value
    return payload


def _scalar(connection: sqlite3.Connection, query: str) -> int:
    row = connection.execute(query).fetchone()
    assert row is not None
    return int(row[0])


def _open_deterministic_gzip_text(path: Path) -> TextIOWrapper:
    raw = path.open("wb")
    compressed = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
    return TextIOWrapper(compressed, encoding="utf-8", newline="")
