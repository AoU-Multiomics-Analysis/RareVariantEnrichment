"""Transcript annotation contracts for gene-matched carrier extraction."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Sequence

from rare_variant_enrichment.annotations import (
    ENSEMBL_CONSEQUENCE_ORDER,
    VariantKey,
    collapse_loftee,
    most_severe_consequence,
    normalize_gene_id,
    parse_consequence_terms,
)
from rare_variant_enrichment.io import write_json


_CONSEQUENCE_RANK = {
    consequence: rank for rank, consequence in enumerate(ENSEMBL_CONSEQUENCE_ORDER)
}
_SPLICE_CORE = {"splice_acceptor_variant", "splice_donor_variant"}
_SPLICE_REGION = {
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}


@dataclass(frozen=True)
class TranscriptCarrierSchema:
    header: tuple[str, ...]
    chromosome: int
    position: int
    ref: int
    alt: int
    gene_id: int
    gene_symbol: int
    consequence: int
    lof: int
    gvs_max_af: int
    revel: int

    @classmethod
    def from_header(cls, header: Sequence[str]) -> "TranscriptCarrierSchema":
        header_tuple = tuple(header)
        required = (
            "chrom", "pos", "ref", "alt", "gene_id", "gene_symbol",
            "consequence", "LoF", "gvs_max_af", "revel",
        )
        missing = [name for name in required if name not in header_tuple]
        duplicates = [name for name in required if header_tuple.count(name) > 1]
        if missing:
            raise ValueError("Missing required transcript columns: " + ", ".join(missing))
        if duplicates:
            raise ValueError("Duplicate required transcript columns: " + ", ".join(duplicates))
        index = {name: header_tuple.index(name) for name in required}
        return cls(
            header_tuple,
            index["chrom"],
            index["pos"],
            index["ref"],
            index["alt"],
            index["gene_id"],
            index["gene_symbol"],
            index["consequence"],
            index["LoF"],
            index["gvs_max_af"],
            index["revel"],
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "header": list(self.header),
            "chromosome": self.chromosome,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "gene_id": self.gene_id,
            "gene_symbol": self.gene_symbol,
            "consequence": self.consequence,
            "lof": self.lof,
            "gvs_max_af": self.gvs_max_af,
            "revel": self.revel,
        }

    def write_json(self, path: Path) -> None:
        write_json(path, self.as_dict())

    @classmethod
    def read_json(cls, path: Path) -> "TranscriptCarrierSchema":
        try:
            with path.open(encoding="utf-8") as handle:
                payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
        except json.JSONDecodeError as error:
            raise ValueError(f"Transcript schema JSON is invalid: {path}") from error
        if not isinstance(payload, dict) or not isinstance(payload.get("header"), list):
            raise ValueError("Transcript schema JSON must contain a header array")
        resolved = cls.from_header(payload["header"])
        if resolved.as_dict() != payload:
            raise ValueError("Transcript schema indices do not match its header")
        return resolved


@dataclass(frozen=True)
class TranscriptCarrierRow:
    key: VariantKey
    gene_id: str
    gene_symbol: str | None
    consequences: tuple[str, ...]
    loftee: str | None
    revel: float | None
    gvs_max_af: float | None


@dataclass(frozen=True)
class CollapsedCarrierAnnotation:
    key: VariantKey
    gene_id: str
    gene_symbol: str | None
    most_severe_consequence: str | None
    all_consequences: tuple[str, ...]
    unknown_consequences: tuple[str, ...]
    loftee: str | None
    revel: float | None
    gvs_max_af: float | None


def parse_optional_unit_interval(value: str, label: str) -> float | None:
    stripped = value.strip()
    if stripped.casefold() in {"", ".", "na", "null"}:
        return None
    try:
        parsed = float(stripped)
    except ValueError as error:
        raise ValueError(f"{label} must be numeric or missing") from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise ValueError(f"{label} must be from 0 through 1")
    return parsed


def parse_transcript_carrier_row(
    fields: Sequence[str], schema: TranscriptCarrierSchema
) -> TranscriptCarrierRow:
    if len(fields) != len(schema.header):
        raise ValueError("Transcript row column count does not match its header")
    try:
        position = int(fields[schema.position])
    except ValueError as error:
        raise ValueError("Transcript position must be an integer") from error
    gene_symbol = fields[schema.gene_symbol].strip()
    raw_loftee = fields[schema.lof].strip().upper()
    return TranscriptCarrierRow(
        key=VariantKey(
            fields[schema.chromosome], position, fields[schema.ref], fields[schema.alt]
        ),
        gene_id=normalize_gene_id(fields[schema.gene_id].strip()),
        gene_symbol=gene_symbol if gene_symbol not in {"", "."} else None,
        consequences=parse_consequence_terms(fields[schema.consequence]),
        loftee=raw_loftee if raw_loftee in {"HC", "LC"} else None,
        revel=parse_optional_unit_interval(fields[schema.revel], "REVEL"),
        gvs_max_af=parse_optional_unit_interval(fields[schema.gvs_max_af], "gvs_max_af"),
    )


def collapse_transcript_rows(
    rows: Sequence[TranscriptCarrierRow],
) -> CollapsedCarrierAnnotation:
    if not rows:
        raise ValueError("At least one transcript row is required")
    key = rows[0].key
    if any(row.key != key for row in rows[1:]):
        raise ValueError("Transcript rows must have the same exact allele")
    normalized_genes = {normalize_gene_id(row.gene_id) for row in rows}
    if len(normalized_genes) != 1:
        raise ValueError("Transcript rows must have the same normalized gene")
    symbols = {row.gene_symbol for row in rows if row.gene_symbol is not None}
    if len(symbols) > 1:
        raise ValueError("Transcript rows have conflicting gene symbols")

    terms = {term for row in rows for term in row.consequences}
    selected, unknown = most_severe_consequence(terms)
    ordered_terms = tuple(
        sorted(terms, key=lambda term: (_CONSEQUENCE_RANK.get(term, len(_CONSEQUENCE_RANK)), term))
    )
    revel_values = [row.revel for row in rows if row.revel is not None]
    frequency_values = [row.gvs_max_af for row in rows if row.gvs_max_af is not None]
    return CollapsedCarrierAnnotation(
        key=key,
        gene_id=next(iter(normalized_genes)),
        gene_symbol=next(iter(symbols)) if symbols else None,
        most_severe_consequence=selected,
        all_consequences=ordered_terms,
        unknown_consequences=tuple(sorted(unknown)),
        loftee=collapse_loftee(row.loftee or "" for row in rows),
        revel=max(revel_values, default=None),
        gvs_max_af=max(frequency_values, default=None),
    )


def initial_variant_classes(annotation: CollapsedCarrierAnnotation) -> tuple[str, ...]:
    classes: list[str] = []
    if annotation.loftee == "HC":
        classes.append("lof_hc")
    if annotation.loftee in {"HC", "LC"}:
        classes.append("lof_hc_or_lc")
    if annotation.most_severe_consequence == "missense_variant":
        classes.append("missense")
    if annotation.most_severe_consequence in _SPLICE_CORE:
        classes.append("splice_core")
    if annotation.most_severe_consequence in _SPLICE_REGION:
        classes.append("splice_region")
    return tuple(classes)


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"Transcript schema JSON contains duplicate key: {key}")
        payload[key] = value
    return payload
