from dataclasses import dataclass
from typing import Literal, Sequence


@dataclass(frozen=True)
class AcClass:
    label: str
    kind: Literal["exact", "cumulative"]
    value: int

    def contains(self, ac: int) -> bool:
        return ac == self.value if self.kind == "exact" else ac <= self.value


@dataclass(frozen=True)
class VariantAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    ac: int
    carriers: tuple[str, ...]


@dataclass(frozen=True, order=True)
class FeatureTss:
    chrom: str
    tss: int
    feature_id: str


def build_ac_classes(
    exact_ac: Sequence[int], cumulative_ac_max: Sequence[int]
) -> list[AcClass]:
    exact = _unique_values(exact_ac)
    cumulative = _unique_values(cumulative_ac_max)
    return [
        *[AcClass(f"AC={value}", "exact", value) for value in exact],
        *[AcClass(f"AC<={value}", "cumulative", value) for value in cumulative],
    ]


def parse_variant_alleles(
    fields: Sequence[str], sample_ids: Sequence[str], shared_samples: set[str]
) -> list[VariantAllele]:
    if len(fields) < 9:
        raise ValueError("VCF record must contain fixed fields and FORMAT")

    chrom, pos_text, _, ref, alt_text, _, _, info, format_text = fields[:9]
    sample_fields = fields[9:]
    if len(sample_fields) != len(sample_ids):
        raise ValueError("VCF record sample count does not match sample IDs")

    try:
        pos = int(pos_text)
    except ValueError as error:
        raise ValueError(f"VCF position is not an integer: {pos_text}") from error

    alts = alt_text.split(",")
    genotype_index = _genotype_index(format_text)
    genotype_ac, carriers = _count_genotypes(
        sample_fields, sample_ids, shared_samples, genotype_index, len(alts)
    )
    info_ac = _info_ac_values(info, len(alts))
    ac_values = genotype_ac if info_ac is None else info_ac

    return [
        VariantAllele(chrom, pos, ref, alt, ac, tuple(alt_carriers))
        for alt, ac, alt_carriers in zip(alts, ac_values, carriers, strict=True)
    ]


def merge_tss_windows(
    features: Sequence[FeatureTss], max_distance: int
) -> list[tuple[str, int, int]]:
    _validate_max_distance(max_distance)
    windows = [
        (feature.chrom, max(0, feature.tss - max_distance - 1), feature.tss + max_distance)
        for feature in sorted(features)
    ]
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in windows:
        if merged and merged[-1][0] == chrom and start < merged[-1][2]:
            previous_chrom, previous_start, previous_end = merged[-1]
            merged[-1] = (previous_chrom, previous_start, max(previous_end, end))
        else:
            merged.append((chrom, start, end))
    return merged


def nearby_features(
    features: Sequence[FeatureTss], position: int, max_distance: int
) -> list[FeatureTss]:
    _validate_max_distance(max_distance)
    return [
        feature for feature in sorted(features) if abs(position - feature.tss) <= max_distance
    ]


def _unique_values(values: Sequence[int]) -> list[int]:
    unique: list[int] = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique


def _genotype_index(format_text: str) -> int:
    format_fields = format_text.split(":")
    try:
        return format_fields.index("GT")
    except ValueError as error:
        raise ValueError("VCF FORMAT field does not contain GT") from error


def _count_genotypes(
    sample_fields: Sequence[str],
    sample_ids: Sequence[str],
    shared_samples: set[str],
    genotype_index: int,
    alt_count: int,
) -> tuple[list[int], list[list[str]]]:
    ac_values = [0] * alt_count
    carriers: list[list[str]] = [[] for _ in range(alt_count)]
    for sample_id, sample_field in zip(sample_ids, sample_fields, strict=True):
        genotype_fields = sample_field.split(":")
        genotype = genotype_fields[genotype_index] if genotype_index < len(genotype_fields) else "."
        alt_indices = _alt_indices(genotype, alt_count)
        for allele_index in alt_indices:
            ac_values[allele_index - 1] += 1
        if sample_id in shared_samples:
            for allele_index in set(alt_indices):
                carriers[allele_index - 1].append(sample_id)
    return ac_values, carriers


def _alt_indices(genotype: str, alt_count: int) -> list[int]:
    if genotype in {"", "."}:
        return []
    alleles = genotype.replace("|", "/").split("/")
    alt_indices: list[int] = []
    for allele in alleles:
        if allele == ".":
            continue
        try:
            allele_index = int(allele)
        except ValueError as error:
            raise ValueError(f"Invalid GT allele index: {allele}") from error
        if allele_index < 0 or allele_index > alt_count:
            raise ValueError(f"GT allele index outside REF/ALT range: {allele_index}")
        if allele_index:
            alt_indices.append(allele_index)
    return alt_indices


def _info_ac_values(info: str, alt_count: int) -> list[int] | None:
    for item in info.split(";"):
        if not item.startswith("AC="):
            continue
        ac_text = item.removeprefix("AC=")
        values = ac_text.split(",")
        if len(values) != alt_count:
            raise ValueError(
                f"INFO/AC cardinality ({len(values)}) does not match ALT count ({alt_count})"
            )
        try:
            return [int(value) for value in values]
        except ValueError as error:
            raise ValueError("INFO/AC values must be integers") from error
    return None


def _validate_max_distance(max_distance: int) -> None:
    if max_distance < 0:
        raise ValueError("Maximum distance must be non-negative")
