from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
from typing import Literal, Sequence, TextIO

from rare_variant_enrichment.io import open_text, read_nonempty_lines, write_json


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
    if not exact and not cumulative:
        raise ValueError("At least one AC class is required")
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


def classify_chromosome(
    vcf_path: Path,
    features_path: Path,
    shared_samples_path: Path,
    chromosome: str,
    exact_ac: Sequence[int],
    cumulative_ac_max: Sequence[int],
    max_distance: int,
    carrier_output: Path,
    regions_output: Path,
    qc_output: Path,
) -> None:
    _, vcf_contigs = _read_vcf_header(vcf_path)
    if not vcf_contigs:
        raise ValueError("VCF header does not declare contigs")
    if chromosome not in vcf_contigs:
        raise ValueError(f"Requested chromosome is absent from VCF: {chromosome}")

    features = [feature for feature in _read_features(features_path) if feature.chrom == chromosome]
    shared_samples = set(read_nonempty_lines(shared_samples_path))
    ac_classes = build_ac_classes(exact_ac, cumulative_ac_max)
    regions = merge_tss_windows(features, max_distance)
    _write_regions(regions_output, regions)

    qc = {
        "alt_alleles": 0,
        "chromosome": chromosome,
        "emitted_keys": 0,
        "extracted_records": 0,
        "feature_count": len(features),
        "merged_region_count": len(regions),
        "missing_genotypes": 0,
        "tabix_query_count": 0,
        "variant_feature_pairs": 0,
    }
    minimum_distances: dict[tuple[str, str, str], int] = {}
    if not features:
        _write_carriers(carrier_output, minimum_distances)
        write_json(qc_output, qc)
        return

    sample_ids = _stream_tabix_records(
        vcf_path,
        regions_output,
        chromosome,
        features,
        max_distance,
        shared_samples,
        ac_classes,
        minimum_distances,
        qc,
    )
    if sample_ids is None:
        raise ValueError("Tabix output does not contain a VCF header")

    qc["emitted_keys"] = len(minimum_distances)
    _write_carriers(carrier_output, minimum_distances)
    write_json(qc_output, qc)


def _read_features(path: Path) -> list[FeatureTss]:
    with open_text(path) as handle:
        header = _read_first_nonempty_line(handle, "Feature TSV is empty")
        if header != ["chrom", "tss", "feature_id"]:
            raise ValueError("Feature TSV header must be chrom, tss, feature_id")

        features: list[FeatureTss] = []
        for line_number, raw_line in enumerate(handle, start=2):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"Feature TSV line {line_number} must have three columns")
            chrom, tss_text, feature_id = fields
            try:
                tss = int(tss_text)
            except ValueError as error:
                raise ValueError(f"Feature TSV line {line_number} has a non-integer TSS") from error
            if not chrom or not feature_id or tss < 1:
                raise ValueError(f"Feature TSV line {line_number} has invalid feature values")
            features.append(FeatureTss(chrom, tss, feature_id))
    return features


def _read_vcf_header(path: Path) -> tuple[list[str], set[str]]:
    contigs: set[str] = set()
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##contig=<"):
                match = re.search(r"(?:^|,)ID=([^,>]+)", line.removeprefix("##contig=<"))
                if match:
                    contigs.add(match.group(1))
            elif line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) < 9:
                    raise ValueError("VCF header must contain fixed columns and FORMAT")
                return fields[9:], contigs
    raise ValueError("VCF does not contain a #CHROM header")


def _stream_tabix_records(
    vcf_path: Path,
    regions_path: Path,
    chromosome: str,
    features: Sequence[FeatureTss],
    max_distance: int,
    shared_samples: set[str],
    ac_classes: Sequence[AcClass],
    minimum_distances: dict[tuple[str, str, str], int],
    qc: dict[str, int | str],
) -> list[str] | None:
    command = ["tabix", "-h", "-R", str(regions_path), str(vcf_path)]
    sample_ids: list[str] | None = None
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        assert process.stdout is not None
        for raw_line in process.stdout:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header = line.split("\t")
                    if len(header) < 9:
                        raise ValueError("Tabix VCF header must contain fixed columns and FORMAT")
                    sample_ids = header[9:]
                continue
            if sample_ids is None:
                raise ValueError("Tabix record appeared before the VCF header")

            fields = line.split("\t")
            qc["extracted_records"] += 1
            qc["missing_genotypes"] += _count_missing_genotypes(fields)
            alleles = parse_variant_alleles(fields, sample_ids, shared_samples)
            qc["alt_alleles"] += len(alleles)
            for allele in alleles:
                if allele.chrom != chromosome:
                    raise ValueError(
                        f"Tabix returned chromosome {allele.chrom} for requested chromosome {chromosome}"
                    )
                matched_features = nearby_features(features, allele.pos, max_distance)
                qc["variant_feature_pairs"] += len(matched_features)
                matching_classes = [ac_class for ac_class in ac_classes if ac_class.contains(allele.ac)]
                for feature in matched_features:
                    distance = abs(allele.pos - feature.tss)
                    for sample_id in allele.carriers:
                        for ac_class in matching_classes:
                            key = (sample_id, feature.feature_id, ac_class.label)
                            previous_distance = minimum_distances.get(key)
                            if previous_distance is None or distance < previous_distance:
                                minimum_distances[key] = distance
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, command)

    qc["tabix_query_count"] = 1
    return sample_ids


def _count_missing_genotypes(fields: Sequence[str]) -> int:
    if len(fields) < 9:
        raise ValueError("VCF record must contain fixed fields and FORMAT")
    genotype_index = _genotype_index(fields[8])
    missing = 0
    for sample_field in fields[9:]:
        genotype_fields = sample_field.split(":")
        genotype = genotype_fields[genotype_index] if genotype_index < len(genotype_fields) else "."
        if not genotype or "." in genotype.replace("|", "/").split("/"):
            missing += 1
    return missing


def _write_regions(path: Path, regions: Sequence[tuple[str, int, int]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\n")
        for chrom, start, end in regions:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def _write_carriers(path: Path, minimum_distances: dict[tuple[str, str, str], int]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
        for (sample_id, feature_id, ac_class), distance in sorted(minimum_distances.items()):
            handle.write(f"{sample_id}\t{feature_id}\t{ac_class}\t{distance}\n")


def _read_first_nonempty_line(handle: TextIO, empty_message: str) -> list[str]:
    for raw_line in handle:
        if raw_line.strip():
            return raw_line.rstrip("\r\n").split("\t")
    raise ValueError(empty_message)


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
