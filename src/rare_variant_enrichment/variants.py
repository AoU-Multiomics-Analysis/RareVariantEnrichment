from bisect import bisect_left, bisect_right
from dataclasses import dataclass
import math
from pathlib import Path
import re
import subprocess
from typing import Callable, Iterator, Literal, MutableMapping, Sequence, TextIO

from rare_variant_enrichment.annotation_storage import VatChunkStore
from rare_variant_enrichment.annotations import (
    GeneAnnotation,
    VatSchema,
    VariantKey,
    build_annotation_classes,
    normalize_gene_id,
)
from rare_variant_enrichment.io import open_text, read_nonempty_lines, write_json
from rare_variant_enrichment.storage import MinimumDistanceStore


_VAT_QC_SUM_KEYS = (
    "vat_rows",
    "duplicate_vat_rows",
    "unique_vat_alleles",
    "unique_vat_allele_gene_pairs",
    "converted_gvs_max_af_values",
    "missing_frequency_alleles",
    "non_numeric_frequency_alleles",
    "out_of_range_frequency_alleles",
    "inconsistent_frequency_alleles",
    "above_maf_threshold_alleles",
    "consequence_terms_parsed",
    "recognized_consequence_terms",
    "unknown_consequence_terms",
    "unknown_consequence_rows",
    "configured_consequence_annotations",
    "loftee_hc_values",
    "loftee_lc_values",
    "loftee_missing_values",
    "loftee_unrecognized_values",
)


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


@dataclass(frozen=True, order=True)
class QueryChunk:
    chromosome: str
    start: int
    end: int

    @property
    def tabix_region(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"


class FeatureTssIndex:
    """One-time positional index for repeated inclusive TSS-window lookups."""

    def __init__(self, features: Sequence[FeatureTss]):
        self.features = tuple(
            sorted(features, key=lambda feature: (feature.tss, feature.chrom, feature.feature_id))
        )
        self.positions = tuple(feature.tss for feature in self.features)

    def nearby(self, position: int, max_distance: int) -> list[FeatureTss]:
        _validate_max_distance(max_distance)
        start = bisect_left(self.positions, position - max_distance)
        end = bisect_right(self.positions, position + max_distance)
        return list(self.features[start:end])


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
    fields: Sequence[str],
    sample_ids: Sequence[str],
    shared_samples: set[str],
    *,
    qc: MutableMapping[str, int] | None = None,
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
    parse_qc = qc if qc is not None else {}
    _initialize_variant_qc(parse_qc)
    parse_qc["alt_alleles"] += len(alts)
    genotype_ac, carriers, called_allele_count, genotypes_complete = _count_genotypes(
        sample_fields,
        sample_ids,
        shared_samples,
        genotype_index,
        len(alts),
        parse_qc,
    )
    info_ac = _info_ac_values(info, len(alts))
    alleles: list[VariantAllele] = []
    for alt, parsed_info_ac, genotype_count, alt_carriers in zip(
        alts, info_ac, genotype_ac, carriers, strict=True
    ):
        if parsed_info_ac is None:
            if called_allele_count == 0:
                parse_qc["unavailable_ac_alt_alleles"] += 1
                continue
            ac = genotype_count
            parse_qc["genotype_ac_fallback_alt_alleles"] += 1
        else:
            ac = parsed_info_ac
            parse_qc["info_ac_alt_alleles"] += 1
            if genotypes_complete:
                parse_qc["info_genotype_ac_compared_alt_alleles"] += 1
                if ac != genotype_count:
                    parse_qc["info_genotype_ac_mismatch_alt_alleles"] += 1
            else:
                parse_qc["info_genotype_ac_unchecked_alt_alleles"] += 1
        alleles.append(VariantAllele(chrom, pos, ref, alt, ac, tuple(alt_carriers)))
    return alleles


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


def build_query_chunks(
    features: Sequence[FeatureTss], max_distance: int, chunk_size_bp: int
) -> list[QueryChunk]:
    _validate_max_distance(max_distance)
    if isinstance(chunk_size_bp, bool) or not isinstance(chunk_size_bp, int) or chunk_size_bp < 1:
        raise ValueError("annotation_chunk_size_bp must be a positive integer")
    chromosomes = {feature.chrom for feature in features}
    if len(chromosomes) > 1:
        raise ValueError("Query chunks require features from one chromosome")
    windows = sorted(
        (feature.chrom, max(1, feature.tss - max_distance), feature.tss + max_distance)
        for feature in features
    )
    return [
        QueryChunk(chromosome, chunk_start, min(end, chunk_start + chunk_size_bp - 1))
        for chromosome, start, end in _merge_inclusive_windows(windows)
        for chunk_start in range(start, end + 1, chunk_size_bp)
    ]


def _merge_inclusive_windows(
    windows: Sequence[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chromosome, start, end in windows:
        if merged and merged[-1][0] == chromosome and start <= merged[-1][2]:
            previous_chromosome, previous_start, previous_end = merged[-1]
            merged[-1] = (previous_chromosome, previous_start, max(previous_end, end))
        else:
            merged.append((chromosome, start, end))
    return merged


def nearby_features(
    features: Sequence[FeatureTss], position: int, max_distance: int
) -> list[FeatureTss]:
    return FeatureTssIndex(features).nearby(position, max_distance)


def classify_chromosome(
    vcf_path: Path,
    vat_path: Path,
    vat_schema_path: Path,
    features_path: Path,
    shared_samples_path: Path,
    chromosome: str,
    exact_ac: Sequence[int],
    cumulative_ac_max: Sequence[int],
    consequence_classes: Sequence[str],
    maximum_gvs_maf: float,
    max_distance: int,
    annotation_chunk_size_bp: int,
    carrier_output: Path,
    regions_output: Path,
    qc_output: Path,
) -> None:
    sample_ids, vcf_contigs = _read_vcf_header(vcf_path)
    if not vcf_contigs:
        raise ValueError("VCF header does not declare contigs")
    if chromosome not in vcf_contigs:
        raise ValueError(f"Requested chromosome is absent from VCF: {chromosome}")

    _validate_maximum_gvs_maf(maximum_gvs_maf)
    vat_schema = VatSchema.read_json(vat_schema_path)
    build_annotation_classes(consequence_classes, vat_schema.lof is not None)
    features = [feature for feature in _read_features(features_path) if feature.chrom == chromosome]
    shared_samples = set(read_nonempty_lines(shared_samples_path))
    ac_classes = build_ac_classes(exact_ac, cumulative_ac_max)
    chunks = build_query_chunks(features, max_distance, annotation_chunk_size_bp)
    _write_query_chunks(regions_output, chunks)

    qc: dict[str, int | float | str] = {
        "alt_alleles": 0,
        "annotation_chunk_count": len(chunks),
        "above_maf_threshold_alleles": 0,
        "baseline_emitted_keys": 0,
        "boundary_variant_feature_pairs": 0,
        "chromosome": chromosome,
        "classified_alt_alleles": 0,
        "configured_consequence_annotations": 0,
        "consequence_emitted_keys": 0,
        "consequence_terms_parsed": 0,
        "converted_gvs_max_af_values": 0,
        "duplicate_vat_rows": 0,
        "emitted_keys": 0,
        "extracted_records": 0,
        "feature_count": len(features),
        "gene_matched_variant_feature_pairs": 0,
        "gene_unmatched_variant_feature_pairs": 0,
        "inconsistent_frequency_alleles": 0,
        "loftee_emitted_keys": 0,
        "loftee_hc_values": 0,
        "loftee_lc_values": 0,
        "loftee_missing_values": 0,
        "loftee_unrecognized_values": 0,
        "merged_region_count": len(merge_tss_windows(features, max_distance)),
        "missing_frequency_alleles": 0,
        "missing_genotypes": 0,
        "non_numeric_frequency_alleles": 0,
        "observed_raw_gvs_max_af": "none",
        "out_of_range_frequency_alleles": 0,
        "recognized_consequence_terms": 0,
        "tabix_query_count": 0,
        "unique_vat_allele_gene_pairs": 0,
        "unique_vat_alleles": 0,
        "unknown_consequence_rows": 0,
        "unknown_consequence_terms": 0,
        "vat_joined_alt_alleles": 0,
        "vat_rows": 0,
        "vat_tabix_query_count": 0,
        "vat_unmatched_alt_alleles": 0,
        "variant_feature_pairs": 0,
        "vcf_tabix_query_count": 0,
    }
    _initialize_variant_qc(qc)  # type: ignore[arg-type]
    if not features:
        with MinimumDistanceStore(carrier_output.parent) as minimum_distances:
            minimum_distances.write_tsv(carrier_output, "sample")
        write_json(qc_output, qc)
        return

    feature_index = FeatureTssIndex(features)
    configured_consequence_set = frozenset(consequence_classes)
    with MinimumDistanceStore(carrier_output.parent) as minimum_distances:
        for chunk in chunks:
            with VatChunkStore(
                carrier_output.parent,
                vat_schema,
                maximum_gvs_maf,
                consequence_classes,
            ) as annotations:
                _stream_tabix_tsv(vat_path, chunk.tabix_region, annotations.ingest)
                qc["vat_tabix_query_count"] += 1
                _merge_annotation_qc(qc, annotations.finalize())

                for fields in _stream_tabix_vcf(vcf_path, chunk.tabix_region):
                    qc["extracted_records"] += 1
                    alleles = parse_variant_alleles(
                        fields, sample_ids, shared_samples, qc=qc  # type: ignore[arg-type]
                    )
                    qc["classified_alt_alleles"] += len(alleles)
                    for allele in alleles:
                        if allele.chrom != chromosome:
                            raise ValueError(
                                f"Tabix returned chromosome {allele.chrom} for requested "
                                f"chromosome {chromosome}"
                            )
                        key = VariantKey(allele.chrom, allele.pos, allele.ref, allele.alt)
                        if not annotations.has_allele(key):
                            qc["vat_unmatched_alt_alleles"] += 1
                            continue
                        qc["vat_joined_alt_alleles"] += 1
                        if annotations.qualifying_maf(key) is None:
                            continue

                        matched_features = feature_index.nearby(allele.pos, max_distance)
                        qc["variant_feature_pairs"] += len(matched_features)
                        for feature in matched_features:
                            distance = abs(allele.pos - feature.tss)
                            if distance == max_distance:
                                qc["boundary_variant_feature_pairs"] += 1
                            gene = normalize_gene_id(feature.feature_id)
                            if annotations.has_gene_annotation(key, gene):
                                qc["gene_matched_variant_feature_pairs"] += 1
                            else:
                                qc["gene_unmatched_variant_feature_pairs"] += 1
                            annotation = annotations.gene_annotation(key, gene)
                            _upsert_allele_carriers(
                                minimum_distances,
                                allele,
                                feature,
                                ac_classes,
                                annotation,
                                configured_consequence_set,
                            )
                qc["vcf_tabix_query_count"] += 1

        qc["emitted_keys"] = minimum_distances.count()
        family_counts = minimum_distances.count_by_annotation_family()
        qc["baseline_emitted_keys"] = family_counts.get("baseline", 0)
        qc["consequence_emitted_keys"] = family_counts.get("consequence", 0)
        qc["loftee_emitted_keys"] = family_counts.get("loftee", 0)
        minimum_distances.write_tsv(carrier_output, "sample")
    qc["tabix_query_count"] = qc["vcf_tabix_query_count"]
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


def _stream_tabix_tsv(
    path: Path, region: str, consume: Callable[[Sequence[str]], None]
) -> None:
    command = ["tabix", str(path), region]
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        assert process.stdout is not None
        for raw_line in process.stdout:
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            consume(line.split("\t"))
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, command)


def _stream_tabix_vcf(path: Path, region: str) -> Iterator[list[str]]:
    command = ["tabix", str(path), region]
    with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        assert process.stdout is not None
        for raw_line in process.stdout:
            line = raw_line.rstrip("\r\n")
            if line and not line.startswith("#"):
                yield line.split("\t")
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, command)


def _upsert_allele_carriers(
    minimum_distances: MinimumDistanceStore,
    allele: VariantAllele,
    feature: FeatureTss,
    ac_classes: Sequence[AcClass],
    annotation: GeneAnnotation,
    configured_consequences: frozenset[str],
) -> None:
    matching_classes = [ac_class for ac_class in ac_classes if ac_class.contains(allele.ac)]
    distance = abs(allele.pos - feature.tss)
    annotation_keys = [("baseline", "all_rare_variants")]
    if annotation.consequence in configured_consequences:
        annotation_keys.append(("consequence", annotation.consequence))
    if annotation.loftee in {"HC", "LC"}:
        annotation_keys.append(("loftee", annotation.loftee))
    for sample_id in allele.carriers:
        for ac_class in matching_classes:
            for family, annotation_class in annotation_keys:
                minimum_distances.upsert(
                    sample_id,
                    feature.feature_id,
                    ac_class.label,
                    family,
                    annotation_class,
                    distance,
                )


def _merge_annotation_qc(
    aggregate: dict[str, int | float | str], chunk: dict[str, int | float | str]
) -> None:
    for key in _VAT_QC_SUM_KEYS:
        aggregate[key] = int(aggregate[key]) + int(chunk[key])
    observed = chunk["observed_raw_gvs_max_af"]
    if observed != "none":
        previous = aggregate["observed_raw_gvs_max_af"]
        aggregate["observed_raw_gvs_max_af"] = (
            float(observed) if previous == "none" else max(float(previous), float(observed))
        )


def _write_query_chunks(path: Path, chunks: Sequence[QueryChunk]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\n")
        for chunk in chunks:
            handle.write(f"{chunk.chromosome}\t{chunk.start - 1}\t{chunk.end}\n")


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
    qc: MutableMapping[str, int],
) -> tuple[list[int], list[list[str]], int, bool]:
    ac_values = [0] * alt_count
    carriers: list[list[str]] = [[] for _ in range(alt_count)]
    called_allele_count = 0
    genotypes_complete = bool(sample_fields)
    for sample_id, sample_field in zip(sample_ids, sample_fields, strict=True):
        genotype_fields = sample_field.split(":")
        genotype = genotype_fields[genotype_index] if genotype_index < len(genotype_fields) else "."
        alt_indices, call_state, called_in_genotype = _parse_genotype(genotype, alt_count)
        called_allele_count += called_in_genotype
        if call_state == "fully_missing":
            qc["fully_missing_genotype_calls"] += 1
            qc["missing_genotypes"] += 1
            genotypes_complete = False
        elif call_state == "partial":
            qc["partial_genotype_calls"] += 1
            qc["missing_genotypes"] += 1
            genotypes_complete = False
        else:
            qc["complete_genotype_calls"] += 1
        for allele_index in alt_indices:
            ac_values[allele_index - 1] += 1
        if sample_id in shared_samples:
            for allele_index in set(alt_indices):
                carriers[allele_index - 1].append(sample_id)
    return ac_values, carriers, called_allele_count, genotypes_complete


def _parse_genotype(
    genotype: str, alt_count: int
) -> tuple[list[int], Literal["complete", "partial", "fully_missing"], int]:
    if genotype in {"", "."}:
        return [], "fully_missing", 0
    alleles = genotype.replace("|", "/").split("/")
    alt_indices: list[int] = []
    called_allele_count = 0
    for allele in alleles:
        if allele == ".":
            continue
        try:
            allele_index = int(allele)
        except ValueError as error:
            raise ValueError(f"Invalid GT allele index: {allele}") from error
        if allele_index < 0 or allele_index > alt_count:
            raise ValueError(f"GT allele index outside REF/ALT range: {allele_index}")
        called_allele_count += 1
        if allele_index:
            alt_indices.append(allele_index)
    if called_allele_count == 0:
        call_state: Literal["complete", "partial", "fully_missing"] = "fully_missing"
    elif called_allele_count < len(alleles):
        call_state = "partial"
    else:
        call_state = "complete"
    return alt_indices, call_state, called_allele_count


def _info_ac_values(info: str, alt_count: int) -> list[int | None]:
    for item in info.split(";"):
        if not item.startswith("AC="):
            continue
        ac_text = item.removeprefix("AC=")
        if ac_text == ".":
            return [None] * alt_count
        values = ac_text.split(",")
        if len(values) != alt_count:
            raise ValueError(
                f"INFO/AC cardinality ({len(values)}) does not match ALT count ({alt_count})"
            )
        parsed: list[int | None] = []
        for value in values:
            if value == ".":
                parsed.append(None)
                continue
            try:
                numeric_value = int(value)
            except ValueError as error:
                raise ValueError("INFO/AC values must be integers or missing") from error
            if numeric_value < 0:
                raise ValueError("INFO/AC values must be non-negative")
            parsed.append(numeric_value)
        return parsed
    return [None] * alt_count


def _initialize_variant_qc(qc: MutableMapping[str, int]) -> None:
    for key in (
        "alt_alleles",
        "complete_genotype_calls",
        "fully_missing_genotype_calls",
        "genotype_ac_fallback_alt_alleles",
        "info_ac_alt_alleles",
        "info_genotype_ac_compared_alt_alleles",
        "info_genotype_ac_mismatch_alt_alleles",
        "info_genotype_ac_unchecked_alt_alleles",
        "missing_genotypes",
        "partial_genotype_calls",
        "unavailable_ac_alt_alleles",
    ):
        qc.setdefault(key, 0)


def _validate_max_distance(max_distance: int) -> None:
    if max_distance < 0:
        raise ValueError("Maximum distance must be non-negative")


def _validate_maximum_gvs_maf(maximum_gvs_maf: float) -> None:
    if (
        isinstance(maximum_gvs_maf, bool)
        or not isinstance(maximum_gvs_maf, (int, float))
        or not math.isfinite(maximum_gvs_maf)
        or not 0.0 <= maximum_gvs_maf <= 0.5
    ):
        raise ValueError("maximum_gvs_maf must be a finite number from 0 to 0.5")
