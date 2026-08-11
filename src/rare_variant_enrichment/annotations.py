"""Immutable VAT schema and transcript-annotation semantic value types."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Literal, Sequence


ENSEMBL_CONSEQUENCE_ORDER = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
)

ENSEMBL_SEVERITY_ORDER_VERSION = "Ensembl release 116"

DEFAULT_CONSEQUENCE_CLASSES: tuple[str, ...] = (
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "synonymous_variant",
    "coding_sequence_variant",
)

_MISSING_CONSEQUENCE_TOKENS = {"", ".", "na", "nan", "null"}
_MISSING_FREQUENCY_TOKENS = {"", ".", "na", "null"}
_ENSEMBL_VERSION = re.compile(r"^(ENSG\d+)(?:\.\d+)$")
_CONSEQUENCE_SPLITTER = re.compile(r"[&,]")
_CONSEQUENCE_RANK = {
    consequence: rank for rank, consequence in enumerate(ENSEMBL_CONSEQUENCE_ORDER)
}


@dataclass(frozen=True)
class VatSchema:
    header: tuple[str, ...]
    chromosome: int
    position: int
    ref: int
    alt: int
    gene_id: int
    consequence: int
    gvs_max_af: int
    lof: int | None

    @classmethod
    def from_header(cls, header: Sequence[str]) -> VatSchema:
        header_tuple = tuple(header)
        required = ("chrom", "pos", "ref", "alt", "gene_id", "consequence", "gvs_max_af")
        missing = [name for name in required if name not in header_tuple]
        if missing:
            raise ValueError(f"Missing required VAT columns: {', '.join(missing)}")

        duplicates = [name for name in required if header_tuple.count(name) > 1]
        if duplicates:
            raise ValueError(f"Duplicate required VAT columns: {', '.join(duplicates)}")

        lof_count = header_tuple.count("LoF")
        if lof_count > 1:
            raise ValueError("Duplicate optional VAT column: LoF")

        indices = {name: header_tuple.index(name) for name in required}
        return cls(
            header=header_tuple,
            chromosome=indices["chrom"],
            position=indices["pos"],
            ref=indices["ref"],
            alt=indices["alt"],
            gene_id=indices["gene_id"],
            consequence=indices["consequence"],
            gvs_max_af=indices["gvs_max_af"],
            lof=header_tuple.index("LoF") if lof_count else None,
        )

    def write_json(self, path: Path) -> None:
        payload = {
            "header": list(self.header),
            "chromosome": self.chromosome,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "gene_id": self.gene_id,
            "consequence": self.consequence,
            "gvs_max_af": self.gvs_max_af,
            "lof": self.lof,
            "lof_enabled": self.lof is not None,
        }
        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True)
            handle.write("\n")

    @classmethod
    def read_json(cls, path: Path) -> VatSchema:
        try:
            with path.open(encoding="utf-8") as handle:
                payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
        except json.JSONDecodeError as error:
            raise ValueError(f"VAT schema JSON is invalid: {path}") from error

        if not isinstance(payload, dict):
            raise ValueError("VAT schema JSON must contain an object")
        expected = {
            "header", "chromosome", "position", "ref", "alt", "gene_id",
            "consequence", "gvs_max_af", "lof", "lof_enabled",
        }
        missing = sorted(expected - payload.keys())
        if missing:
            raise ValueError(f"VAT schema JSON is missing fields: {', '.join(missing)}")

        try:
            schema = cls(
                header=tuple(payload["header"]),
                chromosome=payload["chromosome"],
                position=payload["position"],
                ref=payload["ref"],
                alt=payload["alt"],
                gene_id=payload["gene_id"],
                consequence=payload["consequence"],
                gvs_max_af=payload["gvs_max_af"],
                lof=payload["lof"],
            )
        except (TypeError, ValueError) as error:
            raise ValueError("VAT schema JSON contains invalid schema values") from error

        resolved = cls.from_header(schema.header)
        if schema != resolved or payload["lof_enabled"] != (schema.lof is not None):
            raise ValueError("VAT schema JSON indices do not match its header")
        return schema


@dataclass(frozen=True)
class FrequencyValue:
    status: Literal["valid", "missing", "non_numeric", "out_of_range"]
    maf: float | None
    converted: bool


@dataclass(frozen=True, order=True)
class AnnotationClass:
    family: Literal["baseline", "consequence", "loftee"]
    label: str


@dataclass(frozen=True)
class VariantKey:
    chromosome: str
    position: int
    ref: str
    alt: str


@dataclass(frozen=True)
class GeneAnnotation:
    consequence: str | None
    loftee: str | None


def normalize_gene_id(value: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"Gene ID must be a string: {value!r}")
    normalized = value
    match = _ENSEMBL_VERSION.fullmatch(normalized)
    if match:
        normalized = match.group(1)
    if not normalized:
        raise ValueError(f"Gene ID is empty after normalization: {value!r}")
    return normalized


def parse_consequence_terms(value: str) -> tuple[str, ...]:
    if not isinstance(value, str):
        raise ValueError(f"Consequence value must be a string: {value!r}")
    terms: list[str] = []
    seen: set[str] = set()
    for raw_term in _CONSEQUENCE_SPLITTER.split(value):
        term = raw_term.strip()
        if term.casefold() in _MISSING_CONSEQUENCE_TOKENS or term in seen:
            continue
        terms.append(term)
        seen.add(term)
    return tuple(terms)


def most_severe_consequence(terms: Iterable[str]) -> tuple[str | None, tuple[str, ...]]:
    selected: str | None = None
    selected_rank = len(ENSEMBL_CONSEQUENCE_ORDER)
    unknown: list[str] = []
    unknown_seen: set[str] = set()
    for value in terms:
        for term in parse_consequence_terms(value):
            rank = _CONSEQUENCE_RANK.get(term)
            if rank is None:
                if term not in unknown_seen:
                    unknown.append(term)
                    unknown_seen.add(term)
            elif rank < selected_rank:
                selected = term
                selected_rank = rank
    return selected, tuple(unknown)


def collapse_loftee(values: Iterable[str]) -> str | None:
    normalized = {value.strip().upper() for value in values if isinstance(value, str)}
    if "HC" in normalized:
        return "HC"
    if "LC" in normalized:
        return "LC"
    return None


def parse_gvs_max_af(value: str) -> FrequencyValue:
    if isinstance(value, bool):
        raise ValueError(f"gvs_max_af must not be a Boolean: {value!r}")

    if isinstance(value, str):
        stripped = value.strip()
        if stripped.casefold() in _MISSING_FREQUENCY_TOKENS:
            return FrequencyValue("missing", None, False)
        display_value = value
    else:
        stripped = value
        display_value = value

    try:
        numeric = float(stripped)
    except (TypeError, ValueError):
        return FrequencyValue("non_numeric", None, False)
    if not math.isfinite(numeric):
        raise ValueError(f"gvs_max_af must be finite: {display_value!r}")
    if not 0.0 <= numeric <= 1.0:
        return FrequencyValue("out_of_range", None, False)
    decimal_value = Decimal(stripped) if isinstance(stripped, str) else Decimal(str(stripped))
    maf = float(min(decimal_value, Decimal("1") - decimal_value))
    return FrequencyValue("valid", maf, numeric > 0.5)


def build_annotation_classes(
    consequence_classes: Sequence[str], loftee_enabled: bool
) -> list[AnnotationClass]:
    configured = list(consequence_classes)
    duplicates = sorted({value for value in configured if configured.count(value) > 1})
    if duplicates:
        raise ValueError(f"Duplicate consequence classes: {', '.join(duplicates)}")
    unknown = sorted({value for value in configured if value not in _CONSEQUENCE_RANK})
    if unknown:
        raise ValueError(f"Unknown consequence classes: {', '.join(unknown)}")
    return [
        AnnotationClass("baseline", "all_rare_variants"),
        *[AnnotationClass("consequence", value) for value in configured],
        *(
            [AnnotationClass("loftee", "HC"), AnnotationClass("loftee", "LC")]
            if loftee_enabled
            else []
        ),
    ]


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"VAT schema JSON contains duplicate key: {key}")
        result[key] = value
    return result
