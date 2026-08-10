import math
from pathlib import Path
from typing import Sequence, TextIO

from rare_variant_enrichment.io import open_text, read_nonempty_lines, write_json


MISSING_VALUES = {"", ".", "na", "nan"}


def classify_outlier(value: float, threshold: float, tail: str) -> bool:
    if tail == "absolute":
        return abs(value) >= threshold
    if tail == "positive":
        return value >= threshold
    if tail == "negative":
        return value <= -threshold
    raise ValueError(f"Unsupported tail mode: {tail}")


def prepare_phenotypes(
    phenotype_bed: Path,
    vcf_samples_path: Path,
    chromosomes: Sequence[str],
    z_thresholds: Sequence[float],
    tail: str,
    feature_output: Path,
    sample_output: Path,
    qc_output: Path,
) -> None:
    selected_chromosomes = _validate_chromosomes(chromosomes)
    thresholds = _validate_thresholds(z_thresholds)
    if tail not in {"absolute", "positive", "negative"}:
        raise ValueError(f"Unsupported tail mode: {tail}")

    vcf_samples = read_nonempty_lines(vcf_samples_path)
    _reject_duplicates(vcf_samples, "sample")
    vcf_sample_set = set(vcf_samples)

    with open_text(phenotype_bed) as bed_handle:
        header = _read_header(bed_handle)
        bed_samples = header[4:]
        _reject_duplicates(bed_samples, "sample")
        shared_samples = [sample for sample in bed_samples if sample in vcf_sample_set]
        if not shared_samples:
            raise ValueError("No shared BED and VCF samples")

        shared_columns = [index for index, sample in enumerate(bed_samples) if sample in vcf_sample_set]
        bed_sample_set = set(bed_samples)
        bed_only_samples = [sample for sample in bed_samples if sample not in vcf_sample_set]
        vcf_only_samples = [sample for sample in vcf_samples if sample not in bed_sample_set]

        sample_output.write_text("\n".join(shared_samples) + "\n")
        non_missing_observations = 0
        outlier_observations = {str(threshold): 0 for threshold in thresholds}
        seen_features: set[str] = set()
        seen_chromosomes: set[str] = set()
        feature_count = 0

        with feature_output.open("w", encoding="utf-8") as feature_handle:
            feature_handle.write("chrom\ttss\tfeature_id\n")
            for line_number, raw_line in enumerate(bed_handle, start=2):
                if not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) != len(header):
                    raise ValueError(f"Line {line_number} has {len(fields)} columns; expected {len(header)}")

                chrom, start_text, end_text, feature_id = fields[:4]
                start, end = _parse_tss_interval(start_text, end_text, line_number)
                if not chrom:
                    raise ValueError(f"Line {line_number} has an empty chromosome")
                if not feature_id:
                    raise ValueError(f"Line {line_number} has an empty feature ID")
                if feature_id in seen_features:
                    raise ValueError(f"Duplicate feature ID: {feature_id}")
                seen_features.add(feature_id)
                seen_chromosomes.add(chrom)

                values = [_parse_phenotype_value(value, line_number) for value in fields[4:]]
                if chrom not in selected_chromosomes:
                    continue

                feature_handle.write(f"{chrom}\t{end}\t{feature_id}\n")
                feature_count += 1
                for column_index in shared_columns:
                    value = values[column_index]
                    if value is None:
                        continue
                    non_missing_observations += 1
                    for threshold in thresholds:
                        if classify_outlier(value, threshold, tail):
                            outlier_observations[str(threshold)] += 1

    missing_chromosomes = [chrom for chrom in selected_chromosomes if chrom not in seen_chromosomes]
    if missing_chromosomes:
        raise ValueError(f"Requested chromosomes are absent from BED: {', '.join(missing_chromosomes)}")

    write_json(
        qc_output,
        {
            "bed_only_samples": bed_only_samples,
            "feature_count": feature_count,
            "non_missing_observations": non_missing_observations,
            "outlier_observations": outlier_observations,
            "shared_sample_count": len(shared_samples),
            "vcf_only_samples": vcf_only_samples,
        },
    )


def _read_header(handle: TextIO) -> list[str]:
    for raw_line in handle:
        if raw_line.strip():
            header = raw_line.rstrip("\r\n").split("\t")
            if len(header) < 5:
                raise ValueError("Phenotype BED header must contain four metadata columns and one sample")
            if header[:3] != ["#chr", "start", "end"] or not header[3]:
                raise ValueError("Phenotype BED header must begin with #chr, start, end, and a feature ID column")
            return header
    raise ValueError("Phenotype BED is empty")


def _parse_tss_interval(start_text: str, end_text: str, line_number: int) -> tuple[int, int]:
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as error:
        raise ValueError(f"Line {line_number} has non-integer BED coordinates") from error
    if start < 0 or end < 1:
        raise ValueError(f"Line {line_number} has invalid BED coordinates")
    if end - start != 1:
        raise ValueError(f"Line {line_number} TSS interval must be one base wide")
    return start, end


def _parse_phenotype_value(value_text: str, line_number: int) -> float | None:
    if value_text.strip().lower() in MISSING_VALUES:
        return None
    try:
        value = float(value_text)
    except ValueError as error:
        raise ValueError(f"Line {line_number} has a non-numeric phenotype value: {value_text}") from error
    if not math.isfinite(value):
        raise ValueError(f"Line {line_number} has a non-finite phenotype value: {value_text}")
    return value


def _validate_chromosomes(chromosomes: Sequence[str]) -> list[str]:
    if not chromosomes or any(not chromosome for chromosome in chromosomes):
        raise ValueError("At least one non-empty chromosome is required")
    values = list(chromosomes)
    _reject_duplicates(values, "chromosome")
    return values


def _validate_thresholds(z_thresholds: Sequence[float]) -> list[float]:
    if not z_thresholds:
        raise ValueError("At least one z-score threshold is required")
    thresholds: list[float] = []
    for threshold in z_thresholds:
        try:
            numeric_threshold = float(threshold)
        except (TypeError, ValueError) as error:
            raise ValueError("z-score thresholds must be finite numeric values") from error
        if not math.isfinite(numeric_threshold):
            raise ValueError("z-score thresholds must be finite numeric values")
        thresholds.append(numeric_threshold)
    return thresholds


def _reject_duplicates(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    for value in values:
        if value in seen:
            raise ValueError(f"Duplicate {label} ID: {value}")
        seen.add(value)
