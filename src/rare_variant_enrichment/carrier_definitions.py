"""Validate and apply named carrier-definition rules."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import hashlib
from io import TextIOWrapper
import json
import logging
import math
from pathlib import Path
import re
import sqlite3
import tempfile
from typing import Any

from rare_variant_enrichment import __version__
from rare_variant_enrichment.artifacts import file_artifact
from rare_variant_enrichment.carrier_aggregation import _parse_classes, _validate_audit_row
from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER
from rare_variant_enrichment.io import open_text, write_json


SUPPORTED_BASE_CLASSES = frozenset(
    {"lof_hc", "lof_hc_or_lc", "missense", "splice_core", "splice_region"}
)
_NAME_PATTERN = re.compile(r"[A-Za-z][A-Za-z0-9_.-]*\Z")
_GENE_ID_PATTERN = re.compile(r"(ENSG[0-9]+)(?:\.[0-9]+)?\Z")
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")
_ASCII_WHITESPACE = " \t\r\n\f\v"
CARRIER_DEFINITION_HEADER = (
    "sample_id",
    "gene_id",
    "gene_symbol",
    "carrier_definition",
    "n_variants",
    "variant_ids",
)
PROGRESS_INTERVAL_AUDIT_ROWS = 1_000_000
LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class CarrierDefinition:
    """One ordered carrier rule from a schema-versioned configuration."""

    name: str
    variant_classes: tuple[str, ...]
    minimum_revel: float | None = None

    def matches(self, classes: frozenset[str], revel: float | None) -> bool:
        class_match = any(value in classes for value in self.variant_classes)
        threshold_match = self.minimum_revel is None or (
            revel is not None and revel >= self.minimum_revel
        )
        return class_match and threshold_match


@dataclass(frozen=True)
class CarrierDefinitionConfig:
    """Validated carrier-definition configuration."""

    schema_version: int
    definitions: tuple[CarrierDefinition, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.definitions)


def read_carrier_definition_config(path: Path) -> CarrierDefinitionConfig:
    """Read and strictly validate a schema-version 1 definition file."""
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except json.JSONDecodeError as error:
        raise ValueError(f"Carrier-definition JSON is invalid: {path}") from error

    if not isinstance(payload, dict):
        raise ValueError("Carrier-definition JSON must contain an object")
    unknown_top_level = set(payload) - {"schema_version", "definitions"}
    if unknown_top_level:
        raise ValueError(
            "Carrier-definition JSON has an unknown top-level key: "
            + ", ".join(sorted(unknown_top_level))
        )
    missing_top_level = {"schema_version", "definitions"} - set(payload)
    if missing_top_level:
        raise ValueError(
            "Carrier-definition JSON is missing a top-level key: "
            + ", ".join(sorted(missing_top_level))
        )

    schema_version = payload["schema_version"]
    if isinstance(schema_version, bool) or not isinstance(schema_version, int):
        raise ValueError("schema_version must be the integer 1")
    if schema_version != 1:
        raise ValueError(f"Unsupported carrier-definition schema_version: {schema_version}")

    raw_definitions = payload["definitions"]
    if not isinstance(raw_definitions, list):
        raise ValueError("definitions must be an array")
    if not raw_definitions:
        raise ValueError("definitions must contain at least one carrier definition")

    definitions: list[CarrierDefinition] = []
    names: set[str] = set()
    for index, raw_definition in enumerate(raw_definitions):
        definition = _parse_definition(raw_definition, index)
        if definition.name in names:
            raise ValueError(f"Carrier-definition JSON has duplicate definition name: {definition.name}")
        names.add(definition.name)
        definitions.append(definition)

    return CarrierDefinitionConfig(
        schema_version=schema_version,
        definitions=tuple(definitions),
    )


def build_carrier_definitions(
    audit_path: Path,
    extraction_qc_path: Path,
    config_path: Path,
    output_path: Path,
    qc_path: Path,
    *,
    container_image: str,
) -> None:
    """Build named sample-gene carrier sets from the canonical extraction audit."""
    if not container_image.strip():
        raise ValueError("container_image must not be empty")

    extraction_qc = _read_json_object(extraction_qc_path, "Extraction QC")
    audit_artifact = _validate_audit_binding(audit_path, extraction_qc)
    config = read_carrier_definition_config(config_path)
    LOGGER.info(
        "Building %d carrier definitions from %d audit rows",
        len(config.definitions),
        audit_artifact["row_count"],
    )

    temporary = tempfile.NamedTemporaryFile(
        prefix="carrier-definitions-",
        suffix=".sqlite3",
        dir=output_path.parent,
        delete=False,
    )
    temporary.close()
    database_path = Path(temporary.name)
    duplicate_rows = 0
    present_revel_rows = 0
    missing_revel_rows = 0
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

            CREATE TABLE gene_symbol (
                gene_id TEXT PRIMARY KEY,
                symbol TEXT NOT NULL
            ) WITHOUT ROWID;

            CREATE TABLE definition_variant (
                sample_id TEXT NOT NULL,
                gene_id TEXT NOT NULL,
                definition_index INTEGER NOT NULL,
                variant_id TEXT NOT NULL,
                PRIMARY KEY (sample_id, gene_id, definition_index, variant_id)
            ) WITHOUT ROWID;
            """
        )

        with open_text(audit_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != AUDIT_HEADER:
                raise ValueError("Carrier audit header does not match extraction QC")
            for line_number, row in enumerate(reader, start=2):
                processed_rows = line_number - 1
                if processed_rows % PROGRESS_INTERVAL_AUDIT_ROWS == 0:
                    LOGGER.info("Processed %d carrier audit rows", processed_rows)
                if None in row:
                    raise ValueError(
                        f"Carrier audit line {line_number} has too many columns"
                    )
                if any(value is None for value in row.values()):
                    raise ValueError(
                        f"Carrier audit line {line_number} has too few columns"
                    )
                _validate_audit_row(row)
                sample_id = row["sample_id"].strip(_ASCII_WHITESPACE)
                if not sample_id:
                    raise ValueError(
                        f"Carrier audit line {line_number} has an empty sample ID"
                    )
                gene_id = _normalize_gene_id(row["gene_id"], line_number)
                gene_symbol = row["gene_symbol"].strip(_ASCII_WHITESPACE)
                if gene_symbol == ".":
                    gene_symbol = ""
                if gene_symbol:
                    previous_symbol = connection.execute(
                        "SELECT symbol FROM gene_symbol WHERE gene_id = ?",
                        (gene_id,),
                    ).fetchone()
                    if previous_symbol is not None and previous_symbol[0] != gene_symbol:
                        raise ValueError(f"Conflicting symbols for gene {gene_id}")
                    connection.execute(
                        "INSERT OR IGNORE INTO gene_symbol (gene_id, symbol) VALUES (?, ?)",
                        (gene_id, gene_symbol),
                    )

                canonical_values = []
                for field in AUDIT_HEADER:
                    if field == "sample_id":
                        canonical_values.append(sample_id)
                    elif field == "gene_id":
                        canonical_values.append(gene_id)
                    elif field == "gene_symbol":
                        canonical_values.append(gene_symbol)
                    else:
                        canonical_values.append(row[field])
                row_data = json.dumps(
                    canonical_values,
                    ensure_ascii=False,
                    separators=(",", ":"),
                )
                key = (
                    sample_id,
                    gene_id,
                    row["chrom"],
                    int(row["pos"]),
                    row["ref"],
                    row["alt"],
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
                    (*key[:2], gene_symbol, *key[2:], row["variant_id"], row_data),
                )
                classes = frozenset(_parse_classes(row["variant_classes"]))
                revel = float(row["revel"]) if row["revel"] else None
                if revel is None:
                    missing_revel_rows += 1
                else:
                    present_revel_rows += 1
                for definition_index, definition in enumerate(config.definitions):
                    if definition.matches(classes, revel):
                        connection.execute(
                            """
                            INSERT OR IGNORE INTO definition_variant (
                                sample_id, gene_id, definition_index, variant_id
                            ) VALUES (?, ?, ?, ?)
                            """,
                            (sample_id, gene_id, definition_index, row["variant_id"]),
                        )

        connection.commit()
        deduplicated_rows = _scalar(connection, "SELECT count(*) FROM audit")
        if deduplicated_rows + duplicate_rows != audit_artifact["row_count"]:
            raise ValueError("Carrier audit row counts do not reconcile")

        output_row_count = _write_materialized_table(connection, config, output_path)
        output_artifact = file_artifact(
            output_path,
            "carrier_definitions.tsv.gz",
            CARRIER_DEFINITION_HEADER,
            output_row_count,
        )
        definition_counts = _definition_counts(connection, config)
        write_json(
            qc_path,
            {
                "schema": "aou.carrier-definitions-manifest.v1",
                "schema_version": config.schema_version,
                "definition_order": list(config.names),
                "definitions": [_definition_payload(item) for item in config.definitions],
                "input_artifacts": {
                    "audit": audit_artifact,
                    "extraction_qc": _plain_file_artifact(
                        extraction_qc_path, "variant_carriers.qc.json"
                    ),
                    "definition_config": _plain_file_artifact(
                        config_path, "carrier_definitions.json"
                    ),
                },
                "extraction_provenance": {
                    "vcf_index_provenance": extraction_qc.get(
                        "vcf_index_provenance"
                    ),
                    "transcript_index_provenance": extraction_qc.get(
                        "transcript_index_provenance"
                    ),
                    "quality_or_frequency_filters_applied": extraction_qc.get(
                        "quality_or_frequency_filters_applied"
                    ),
                },
                "audit_counts": {
                    "input_rows": audit_artifact["row_count"],
                    "deduplicated_rows": deduplicated_rows,
                    "duplicate_rows": duplicate_rows,
                    "present_revel_rows": present_revel_rows,
                    "missing_revel_rows": missing_revel_rows,
                },
                "definition_counts": definition_counts,
                "output_artifact": output_artifact,
                "provenance": {
                    "package_versions": {"rare_variant_enrichment": __version__},
                    "container_image": container_image,
                },
            },
            sort_keys=False,
        )
        LOGGER.info(
            "Built %d carrier-definition rows from %d deduplicated audit rows",
            output_row_count,
            deduplicated_rows,
        )
    finally:
        if "connection" in locals():
            connection.close()
        database_path.unlink(missing_ok=True)


def _parse_definition(raw_definition: Any, index: int) -> CarrierDefinition:
    if not isinstance(raw_definition, dict):
        raise ValueError(f"Definition {index} must be an object")
    allowed_keys = {"name", "variant_classes", "minimum_revel"}
    unknown_keys = set(raw_definition) - allowed_keys
    if unknown_keys:
        raise ValueError(
            f"Definition {index} has an unknown definition key: "
            + ", ".join(sorted(unknown_keys))
        )
    missing_keys = {"name", "variant_classes"} - set(raw_definition)
    if missing_keys:
        raise ValueError(
            f"Definition {index} is missing a required key: "
            + ", ".join(sorted(missing_keys))
        )

    name = raw_definition["name"]
    if not isinstance(name, str) or _NAME_PATTERN.fullmatch(name) is None:
        raise ValueError(f"Definition {index} has an invalid name")

    raw_classes = raw_definition["variant_classes"]
    if not isinstance(raw_classes, list) or not raw_classes:
        raise ValueError(f"Definition {name} must have at least one variant class")
    if not all(isinstance(value, str) for value in raw_classes):
        raise ValueError(f"Definition {name} has a non-string variant class")
    if len(raw_classes) != len(set(raw_classes)):
        raise ValueError(f"Definition {name} has a duplicate variant class")
    invalid_classes = sorted(set(raw_classes) - SUPPORTED_BASE_CLASSES)
    if invalid_classes:
        raise ValueError(
            f"Definition {name} has an unsupported variant class: "
            + ", ".join(invalid_classes)
        )

    minimum_revel = None
    if "minimum_revel" in raw_definition:
        minimum_revel = raw_definition["minimum_revel"]
        if isinstance(minimum_revel, bool) or not isinstance(minimum_revel, (int, float)):
            raise ValueError(f"Definition {name} minimum_revel must be numeric")
        minimum_revel = float(minimum_revel)
        if not math.isfinite(minimum_revel) or not 0.0 <= minimum_revel <= 1.0:
            raise ValueError(f"Definition {name} minimum_revel must be from 0 through 1")

    return CarrierDefinition(
        name=name,
        variant_classes=tuple(raw_classes),
        minimum_revel=minimum_revel,
    )


def _validate_audit_binding(
    audit_path: Path, extraction_qc: dict[str, Any]
) -> dict[str, object]:
    raw_artifact = extraction_qc.get("audit_artifact")
    if not isinstance(raw_artifact, dict):
        raise ValueError("Extraction QC is missing audit_artifact")
    required_keys = {"logical_name", "header", "row_count", "size_bytes", "sha256"}
    if set(raw_artifact) != required_keys:
        raise ValueError("Extraction QC audit_artifact has invalid keys")
    if raw_artifact["logical_name"] != "variant_carrier_audit.tsv.gz":
        raise ValueError("Extraction QC audit_artifact has an invalid logical name")
    if raw_artifact["header"] != list(AUDIT_HEADER):
        raise ValueError("Extraction QC audit_artifact has an invalid header")
    for key in ("row_count", "size_bytes"):
        value = raw_artifact[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"Extraction QC audit_artifact has an invalid {key}")
    expected_digest = raw_artifact["sha256"]
    if not isinstance(expected_digest, str) or _SHA256_PATTERN.fullmatch(expected_digest) is None:
        raise ValueError("Extraction QC audit_artifact has an invalid SHA-256")

    actual_digest = _sha256(audit_path)
    if actual_digest != expected_digest:
        raise ValueError("Carrier audit SHA-256 does not match extraction QC")
    if audit_path.stat().st_size != raw_artifact["size_bytes"]:
        raise ValueError("Carrier audit byte size does not match extraction QC")
    header, row_count = _scan_audit(audit_path)
    if header != AUDIT_HEADER:
        raise ValueError("Carrier audit header does not match extraction QC")
    if row_count != raw_artifact["row_count"]:
        raise ValueError("Carrier audit row count does not match extraction QC")
    return {
        "logical_name": raw_artifact["logical_name"],
        "header": list(header),
        "row_count": row_count,
        "size_bytes": audit_path.stat().st_size,
        "sha256": actual_digest,
    }


def _scan_audit(path: Path) -> tuple[tuple[str, ...], int]:
    try:
        with open_text(path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = tuple(next(reader))
            except StopIteration as error:
                raise ValueError("Carrier audit is empty") from error
            row_count = sum(1 for _ in reader)
    except (OSError, UnicodeError) as error:
        raise ValueError("Carrier audit cannot be read") from error
    return header, row_count


def _normalize_gene_id(value: str, line_number: int) -> str:
    normalized = value.strip(_ASCII_WHITESPACE)
    match = _GENE_ID_PATTERN.fullmatch(normalized)
    if match is None:
        raise ValueError(
            f"Carrier audit line {line_number} has an invalid Ensembl gene ID"
        )
    return match.group(1)


def _write_materialized_table(
    connection: sqlite3.Connection,
    config: CarrierDefinitionConfig,
    path: Path,
) -> int:
    output_rows = 0
    with _open_deterministic_gzip_text(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(CARRIER_DEFINITION_HEADER)
        rows = connection.execute(
            """
            SELECT sample_id, gene_id, gene_symbol, definition_index,
                   count(*), group_concat(variant_id, ',')
            FROM (
                SELECT dv.sample_id, dv.gene_id,
                       coalesce(gs.symbol, '') AS gene_symbol,
                       dv.definition_index, dv.variant_id
                FROM definition_variant AS dv
                LEFT JOIN gene_symbol AS gs ON gs.gene_id = dv.gene_id
                ORDER BY dv.sample_id, dv.gene_id, dv.definition_index, dv.variant_id
            )
            GROUP BY sample_id, gene_id, gene_symbol, definition_index
            ORDER BY sample_id, gene_id, definition_index
            """
        )
        for sample_id, gene_id, gene_symbol, definition_index, count, variant_ids in rows:
            writer.writerow(
                (
                    sample_id,
                    gene_id,
                    gene_symbol,
                    config.definitions[int(definition_index)].name,
                    count,
                    variant_ids,
                )
            )
            output_rows += 1
    return output_rows


def _definition_counts(
    connection: sqlite3.Connection, config: CarrierDefinitionConfig
) -> dict[str, object]:
    counts: dict[str, object] = {}
    for definition_index, definition in enumerate(config.definitions):
        values = connection.execute(
            """
            SELECT count(*), count(DISTINCT variant_id),
                   count(DISTINCT sample_id || char(0) || gene_id)
            FROM definition_variant
            WHERE definition_index = ?
            """,
            (definition_index,),
        ).fetchone()
        assert values is not None
        pair_count = int(values[2])
        counts[definition.name] = {
            "matched_audit_rows": int(values[0]),
            "distinct_variants": int(values[1]),
            "distinct_sample_gene_pairs": pair_count,
            "output_rows": pair_count,
        }
    return counts


def _definition_payload(definition: CarrierDefinition) -> dict[str, object]:
    payload: dict[str, object] = {
        "name": definition.name,
        "variant_classes": list(definition.variant_classes),
    }
    if definition.minimum_revel is not None:
        payload["minimum_revel"] = definition.minimum_revel
    return payload


def _plain_file_artifact(path: Path, logical_name: str) -> dict[str, object]:
    return {
        "logical_name": logical_name,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path, context: str) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except json.JSONDecodeError as error:
        raise ValueError(f"{context} JSON is invalid: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"{context} JSON must contain an object")
    return payload


def _scalar(connection: sqlite3.Connection, query: str) -> int:
    row = connection.execute(query).fetchone()
    assert row is not None
    return int(row[0])


def _open_deterministic_gzip_text(path: Path) -> TextIOWrapper:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0)
    return TextIOWrapper(compressed, encoding="utf-8", newline="")


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"Carrier-definition JSON contains duplicate JSON key: {key}")
        payload[key] = value
    return payload
