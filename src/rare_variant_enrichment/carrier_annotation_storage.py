"""Bounded SQLite storage for one transcript-annotation coordinate chunk."""

from __future__ import annotations

import json
import sqlite3
import tempfile
from pathlib import Path
from types import TracebackType
from typing import Sequence

from rare_variant_enrichment.annotations import VariantKey
from rare_variant_enrichment.carrier_annotations import (
    CollapsedCarrierAnnotation,
    TranscriptCarrierSchema,
    collapse_transcript_rows,
    parse_transcript_carrier_row,
)


_MISSING_GENE_IDS = {"", ".", "na", "nan", "null"}


class CarrierAnnotationChunkStore:
    """Store and collapse transcript rows for one bounded tabix query."""

    def __init__(self, directory: Path, schema: TranscriptCarrierSchema):
        temporary = tempfile.NamedTemporaryFile(
            prefix="carrier-annotation-chunk-",
            suffix=".sqlite3",
            dir=directory,
            delete=False,
        )
        temporary.close()
        self.path = Path(temporary.name)
        self.schema = schema
        self._closed = False
        self._finalized = False
        self._transcript_rows = 0
        self._duplicate_rows = 0
        self._missing_gene_rows = 0
        try:
            self.connection = sqlite3.connect(self.path)
            self.connection.execute("PRAGMA journal_mode=OFF")
            self.connection.execute("PRAGMA synchronous=OFF")
            self.connection.execute("PRAGMA temp_store=FILE")
            self.connection.execute("PRAGMA cache_size=-4096")
            self.connection.execute(
                """
                CREATE TABLE transcript_row (
                    chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    gene_id TEXT NOT NULL,
                    row_data TEXT NOT NULL PRIMARY KEY
                ) WITHOUT ROWID
                """
            )
        except BaseException:
            self.path.unlink(missing_ok=True)
            raise

    def __enter__(self) -> "CarrierAnnotationChunkStore":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            self.connection.close()
        finally:
            self.path.unlink(missing_ok=True)

    def ingest(self, fields: Sequence[str]) -> None:
        self._require_open()
        if self._finalized:
            raise RuntimeError("Carrier annotation store is already finalized")
        if len(fields) != len(self.schema.header):
            raise ValueError("Transcript row column count does not match its header")
        self._transcript_rows += 1
        if fields[self.schema.gene_id].strip().casefold() in _MISSING_GENE_IDS:
            self._missing_gene_rows += 1
            return
        parsed = parse_transcript_carrier_row(fields, self.schema)
        row_data = json.dumps(list(fields), ensure_ascii=False, separators=(",", ":"))
        inserted = self.connection.execute(
            """
            INSERT OR IGNORE INTO transcript_row (
                chromosome, position, ref, alt, gene_id, row_data
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                parsed.key.chromosome,
                parsed.key.position,
                parsed.key.ref,
                parsed.key.alt,
                parsed.gene_id,
                row_data,
            ),
        )
        if inserted.rowcount == 0:
            self._duplicate_rows += 1

    def finalize(self) -> dict[str, int]:
        self._require_open()
        self.connection.commit()
        self._finalized = True
        allele_count = self.connection.execute(
            """
            SELECT count(*) FROM (
                SELECT chromosome, position, ref, alt FROM transcript_row
                GROUP BY chromosome, position, ref, alt
            )
            """
        ).fetchone()
        pair_count = self.connection.execute(
            """
            SELECT count(*) FROM (
                SELECT chromosome, position, ref, alt, gene_id FROM transcript_row
                GROUP BY chromosome, position, ref, alt, gene_id
            )
            """
        ).fetchone()
        assert allele_count is not None and pair_count is not None
        return {
            "transcript_rows": self._transcript_rows,
            "duplicate_transcript_rows": self._duplicate_rows,
            "missing_gene_transcript_rows": self._missing_gene_rows,
            "unique_annotation_alleles": int(allele_count[0]),
            "unique_annotation_allele_gene_pairs": int(pair_count[0]),
        }

    def annotations_for_allele(
        self, key: VariantKey
    ) -> tuple[CollapsedCarrierAnnotation, ...]:
        self._require_open()
        if not self._finalized:
            raise RuntimeError("Carrier annotation store must be finalized before query")
        rows = self.connection.execute(
            """
            SELECT gene_id, row_data
            FROM transcript_row
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
            ORDER BY gene_id, row_data
            """,
            (key.chromosome, key.position, key.ref, key.alt),
        )
        grouped: dict[str, list] = {}
        for gene_id, row_data in rows:
            parsed = parse_transcript_carrier_row(json.loads(row_data), self.schema)
            grouped.setdefault(gene_id, []).append(parsed)
        return tuple(collapse_transcript_rows(grouped[gene]) for gene in sorted(grouped))

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("Carrier annotation store is closed")
