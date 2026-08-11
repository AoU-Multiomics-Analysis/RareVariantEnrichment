"""Bounded, disk-backed transcript annotation aggregation for one VAT chunk."""

from __future__ import annotations

import json
import math
import sqlite3
import tempfile
from pathlib import Path
from types import TracebackType
from typing import Sequence

from rare_variant_enrichment.annotations import (
    ENSEMBL_CONSEQUENCE_ORDER,
    GeneAnnotation,
    VatSchema,
    VariantKey,
    collapse_loftee,
    most_severe_consequence,
    normalize_gene_id,
    parse_consequence_terms,
    parse_gvs_max_af,
)


_MISSING_GENE_IDS = {"", ".", "na", "nan", "null"}
_FREQUENCY_TOLERANCE = 1e-12


class VatChunkStore:
    """Aggregate one coordinate chunk without retaining its transcript rows in memory."""

    def __init__(
        self,
        directory: Path,
        schema: VatSchema,
        maximum_gvs_maf: float,
        configured_consequences: Sequence[str],
    ):
        if (
            isinstance(maximum_gvs_maf, bool)
            or not isinstance(maximum_gvs_maf, (int, float))
            or not math.isfinite(maximum_gvs_maf)
            or not 0.0 <= maximum_gvs_maf <= 0.5
        ):
            raise ValueError("maximum_gvs_maf must be a finite number from 0 to 0.5")

        temporary = tempfile.NamedTemporaryFile(
            prefix="vat-chunk-",
            suffix=".sqlite3",
            dir=directory,
            delete=False,
        )
        temporary.close()
        self.path = Path(temporary.name)
        self.schema = schema
        self.maximum_gvs_maf = float(maximum_gvs_maf)
        self.configured_consequences = frozenset(configured_consequences)
        self._closed = False
        self._finalized = False
        self._qc: dict[str, int | float | str] | None = None
        self._vat_rows = 0
        self._duplicate_vat_rows = 0
        self._observed_raw_gvs_max_af: float | None = None
        self._converted_gvs_max_af_values = 0
        self._consequence_terms_parsed = 0
        self._recognized_consequence_terms = 0
        self._loftee_hc_values = 0
        self._loftee_lc_values = 0
        self._loftee_missing_values = 0
        self._loftee_unrecognized_values = 0

        try:
            self.connection = sqlite3.connect(self.path)
            self.connection.execute("PRAGMA journal_mode=OFF")
            self.connection.execute("PRAGMA synchronous=OFF")
            self.connection.execute("PRAGMA temp_store=FILE")
            self.connection.execute("PRAGMA cache_size=-4096")
            self.connection.executescript(
                """
                CREATE TABLE variant_frequency (
                    chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    raw_af REAL,
                    maf REAL,
                    status TEXT NOT NULL,
                    converted INTEGER NOT NULL,
                    PRIMARY KEY (chromosome, position, ref, alt)
                ) WITHOUT ROWID;

                CREATE TABLE gene_annotation (
                    chromosome TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    ref TEXT NOT NULL,
                    alt TEXT NOT NULL,
                    gene_id TEXT NOT NULL,
                    consequence_rank INTEGER,
                    consequence TEXT,
                    loftee_rank INTEGER,
                    loftee TEXT,
                    PRIMARY KEY (chromosome, position, ref, alt, gene_id)
                ) WITHOUT ROWID;

                CREATE TABLE unknown_consequence (
                    term TEXT PRIMARY KEY,
                    row_count INTEGER NOT NULL
                ) WITHOUT ROWID;

                CREATE TABLE transcript_row (
                    row_data TEXT PRIMARY KEY
                ) WITHOUT ROWID;
                """
            )
        except BaseException:
            self.path.unlink(missing_ok=True)
            raise

    def __enter__(self) -> "VatChunkStore":
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
        self._require_ingestable()
        key = self._variant_key(fields)
        serialized_row = json.dumps(list(fields), ensure_ascii=False, separators=(",", ":"))
        inserted = self.connection.execute(
            "INSERT OR IGNORE INTO transcript_row (row_data) VALUES (?)",
            (serialized_row,),
        )
        if inserted.rowcount == 0:
            self._duplicate_vat_rows += 1
        frequency_status, raw_af, maf, converted = self._frequency(fields[self.schema.gvs_max_af])
        self._vat_rows += 1
        if raw_af is not None:
            self._observed_raw_gvs_max_af = max(
                raw_af,
                self._observed_raw_gvs_max_af if self._observed_raw_gvs_max_af is not None else raw_af,
            )
        if converted:
            self._converted_gvs_max_af_values += 1
        self._upsert_frequency(key, frequency_status, raw_af, maf, converted)

        terms = parse_consequence_terms(fields[self.schema.consequence])
        selected, unknown_terms = most_severe_consequence(terms)
        self._consequence_terms_parsed += len(terms)
        self._recognized_consequence_terms += len(terms) - len(unknown_terms)
        for term in unknown_terms:
            self.connection.execute(
                """
                INSERT INTO unknown_consequence (term, row_count) VALUES (?, 1)
                ON CONFLICT (term) DO UPDATE SET row_count = row_count + 1
                """,
                (term,),
            )

        raw_gene_id = fields[self.schema.gene_id]
        if raw_gene_id.casefold() in _MISSING_GENE_IDS:
            return
        gene_id = normalize_gene_id(raw_gene_id)
        loftee, loftee_rank = self._loftee(fields)
        consequence_rank = (
            ENSEMBL_CONSEQUENCE_ORDER.index(selected) if selected is not None else None
        )
        self.connection.execute(
            """
            INSERT INTO gene_annotation (
                chromosome, position, ref, alt, gene_id,
                consequence_rank, consequence, loftee_rank, loftee
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT (chromosome, position, ref, alt, gene_id) DO UPDATE SET
                consequence_rank = CASE
                    WHEN excluded.consequence_rank IS NULL THEN gene_annotation.consequence_rank
                    WHEN gene_annotation.consequence_rank IS NULL
                         OR excluded.consequence_rank < gene_annotation.consequence_rank
                    THEN excluded.consequence_rank
                    ELSE gene_annotation.consequence_rank
                END,
                consequence = CASE
                    WHEN excluded.consequence_rank IS NULL THEN gene_annotation.consequence
                    WHEN gene_annotation.consequence_rank IS NULL
                         OR excluded.consequence_rank < gene_annotation.consequence_rank
                    THEN excluded.consequence
                    ELSE gene_annotation.consequence
                END,
                loftee_rank = CASE
                    WHEN excluded.loftee_rank IS NULL THEN gene_annotation.loftee_rank
                    WHEN gene_annotation.loftee_rank IS NULL
                         OR excluded.loftee_rank < gene_annotation.loftee_rank
                    THEN excluded.loftee_rank
                    ELSE gene_annotation.loftee_rank
                END,
                loftee = CASE
                    WHEN excluded.loftee_rank IS NULL THEN gene_annotation.loftee
                    WHEN gene_annotation.loftee_rank IS NULL
                         OR excluded.loftee_rank < gene_annotation.loftee_rank
                    THEN excluded.loftee
                    ELSE gene_annotation.loftee
                END
            """,
            (*self._key_fields(key), gene_id, consequence_rank, selected, loftee_rank, loftee),
        )

    def finalize(self) -> dict[str, int | float | str]:
        self._require_open()
        if self._finalized:
            assert self._qc is not None
            return dict(self._qc)

        self.connection.commit()
        frequency_counts = self.connection.execute(
            """
            SELECT
                count(*),
                sum(status = 'missing'),
                sum(status = 'non_numeric'),
                sum(status = 'out_of_range'),
                sum(status = 'inconsistent'),
                sum(status = 'valid' AND maf > ?)
            FROM variant_frequency
            """,
            (self.maximum_gvs_maf,),
        ).fetchone()
        unknown_counts = self.connection.execute(
            "SELECT count(*), coalesce(sum(row_count), 0) FROM unknown_consequence"
        ).fetchone()
        pair_count = self.connection.execute("SELECT count(*) FROM gene_annotation").fetchone()
        configured_count = self._configured_consequence_count()
        assert frequency_counts is not None and unknown_counts is not None and pair_count is not None

        self._qc = {
            "vat_rows": self._vat_rows,
            "duplicate_vat_rows": self._duplicate_vat_rows,
            "unique_vat_alleles": int(frequency_counts[0]),
            "unique_vat_allele_gene_pairs": int(pair_count[0]),
            "observed_raw_gvs_max_af": (
                self._observed_raw_gvs_max_af
                if self._observed_raw_gvs_max_af is not None
                else "none"
            ),
            "converted_gvs_max_af_values": self._converted_gvs_max_af_values,
            "missing_frequency_alleles": int(frequency_counts[1] or 0),
            "non_numeric_frequency_alleles": int(frequency_counts[2] or 0),
            "out_of_range_frequency_alleles": int(frequency_counts[3] or 0),
            "inconsistent_frequency_alleles": int(frequency_counts[4] or 0),
            "above_maf_threshold_alleles": int(frequency_counts[5] or 0),
            "consequence_terms_parsed": self._consequence_terms_parsed,
            "recognized_consequence_terms": self._recognized_consequence_terms,
            "unknown_consequence_terms": int(unknown_counts[0]),
            "unknown_consequence_rows": int(unknown_counts[1]),
            "configured_consequence_annotations": configured_count,
            "loftee_hc_values": self._loftee_hc_values,
            "loftee_lc_values": self._loftee_lc_values,
            "loftee_missing_values": self._loftee_missing_values,
            "loftee_unrecognized_values": self._loftee_unrecognized_values,
        }
        self._finalized = True
        return dict(self._qc)

    def qualifying_maf(self, key: VariantKey) -> float | None:
        self._require_finalized()
        row = self.connection.execute(
            """
            SELECT maf FROM variant_frequency
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
              AND status = 'valid' AND maf <= ?
            """,
            (*self._key_fields(key), self.maximum_gvs_maf),
        ).fetchone()
        return float(row[0]) if row is not None else None

    def has_allele(self, key: VariantKey) -> bool:
        self._require_finalized()
        row = self.connection.execute(
            """
            SELECT 1 FROM variant_frequency
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
            """,
            self._key_fields(key),
        ).fetchone()
        return row is not None

    def has_gene_annotation(self, key: VariantKey, normalized_gene_id: str) -> bool:
        self._require_finalized()
        row = self.connection.execute(
            """
            SELECT 1 FROM gene_annotation
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ? AND gene_id = ?
            """,
            (*self._key_fields(key), normalized_gene_id),
        ).fetchone()
        return row is not None

    def gene_annotation(self, key: VariantKey, normalized_gene_id: str) -> GeneAnnotation:
        self._require_finalized()
        row = self.connection.execute(
            """
            SELECT consequence, loftee FROM gene_annotation
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ? AND gene_id = ?
            """,
            (*self._key_fields(key), normalized_gene_id),
        ).fetchone()
        if row is None:
            return GeneAnnotation(None, None)
        consequence, loftee = row
        return GeneAnnotation(
            str(consequence) if consequence is not None else None,
            str(loftee) if loftee is not None else None,
        )

    def _variant_key(self, fields: Sequence[str]) -> VariantKey:
        required_index = max(
            self.schema.chromosome,
            self.schema.position,
            self.schema.ref,
            self.schema.alt,
            self.schema.gene_id,
            self.schema.consequence,
            self.schema.gvs_max_af,
            self.schema.lof if self.schema.lof is not None else 0,
        )
        if len(fields) <= required_index:
            raise ValueError("VAT row has fewer fields than its header")
        chromosome = fields[self.schema.chromosome]
        ref = fields[self.schema.ref]
        alt = fields[self.schema.alt]
        if not chromosome or chromosome == ".":
            raise ValueError("VAT chromosome is missing")
        if not ref or ref == ".":
            raise ValueError("VAT reference allele is missing")
        if not alt or alt == ".":
            raise ValueError("VAT alternate allele is missing")
        try:
            position = int(fields[self.schema.position])
        except (TypeError, ValueError) as error:
            raise ValueError(f"VAT position is invalid: {fields[self.schema.position]!r}") from error
        if position < 1:
            raise ValueError(f"VAT position is invalid: {fields[self.schema.position]!r}")
        return VariantKey(chromosome, position, ref, alt)

    def _frequency(self, value: str) -> tuple[str, float | None, float | None, bool]:
        try:
            parsed = parse_gvs_max_af(value)
        except ValueError:
            return "non_numeric", None, None, False
        if parsed.status != "valid":
            return parsed.status, None, None, False
        raw_af = float(value.strip())
        return "valid", raw_af, parsed.maf, parsed.converted

    def _upsert_frequency(
        self,
        key: VariantKey,
        status: str,
        raw_af: float | None,
        maf: float | None,
        converted: bool,
    ) -> None:
        existing = self.connection.execute(
            """
            SELECT raw_af, status, converted FROM variant_frequency
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
            """,
            self._key_fields(key),
        ).fetchone()
        if existing is None:
            self.connection.execute(
                """
                INSERT INTO variant_frequency (
                    chromosome, position, ref, alt, raw_af, maf, status, converted
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (*self._key_fields(key), raw_af, maf, status, int(converted)),
            )
            return

        previous_raw_af, previous_status, previous_converted = existing
        if previous_status == "inconsistent":
            return
        if status == "valid" and previous_raw_af is not None and raw_af is not None:
            if abs(float(previous_raw_af) - raw_af) > _FREQUENCY_TOLERANCE:
                self.connection.execute(
                    """
                    UPDATE variant_frequency
                    SET status = 'inconsistent', converted = ?
                    WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
                    """,
                    (int(bool(previous_converted) or converted), *self._key_fields(key)),
                )
                return
        elif status == "valid" and maf is not None:
            self.connection.execute(
                """
                UPDATE variant_frequency
                SET raw_af = ?, maf = ?, converted = ?
                WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
                """,
                (
                    raw_af,
                    maf,
                    int(bool(previous_converted) or converted),
                    *self._key_fields(key),
                ),
            )
            return
        elif status != "valid":
            self.connection.execute(
                """
                UPDATE variant_frequency
                SET status = ?, converted = ?
                WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
                """,
                (status, int(bool(previous_converted) or converted), *self._key_fields(key)),
            )
            return

        self.connection.execute(
            """
            UPDATE variant_frequency SET converted = ?
            WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
            """,
            (int(bool(previous_converted) or converted), *self._key_fields(key)),
        )

    def _loftee(self, fields: Sequence[str]) -> tuple[str | None, int | None]:
        if self.schema.lof is None:
            return None, None
        value = fields[self.schema.lof]
        selected = collapse_loftee([value])
        if selected == "HC":
            self._loftee_hc_values += 1
            return selected, 0
        if selected == "LC":
            self._loftee_lc_values += 1
            return selected, 1
        if value.strip().casefold() in _MISSING_GENE_IDS:
            self._loftee_missing_values += 1
        else:
            self._loftee_unrecognized_values += 1
        return None, None

    def _configured_consequence_count(self) -> int:
        if not self.configured_consequences:
            return 0
        placeholders = ", ".join("?" for _ in self.configured_consequences)
        row = self.connection.execute(
            "SELECT count(*) FROM gene_annotation WHERE consequence IN (" + placeholders + ")",
            tuple(sorted(self.configured_consequences)),
        ).fetchone()
        assert row is not None
        return int(row[0])

    @staticmethod
    def _key_fields(key: VariantKey) -> tuple[str, int, str, str]:
        return key.chromosome, key.position, key.ref, key.alt

    def _require_ingestable(self) -> None:
        self._require_open()
        if self._finalized:
            raise RuntimeError("VAT chunk store cannot ingest after finalize")

    def _require_finalized(self) -> None:
        self._require_open()
        if not self._finalized:
            raise RuntimeError("VAT chunk store must finalize before lookup")

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("VAT chunk store is closed")
