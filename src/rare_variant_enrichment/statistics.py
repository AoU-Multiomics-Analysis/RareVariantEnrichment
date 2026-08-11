import csv
import json
import math
import platform
from pathlib import Path
import re
import sqlite3
from typing import Sequence

from rare_variant_enrichment import __version__
from rare_variant_enrichment.aggregation import _iter_carrier_file
from rare_variant_enrichment.annotations import AnnotationClass, build_annotation_classes
from rare_variant_enrichment.io import open_text, read_nonempty_lines, write_json
from rare_variant_enrichment.phenotypes import (
    _parse_phenotype_value,
    _parse_tss_interval,
    _read_header,
    _reject_duplicates,
    _validate_thresholds,
    classify_outlier,
)
from rare_variant_enrichment.variants import AcClass, build_ac_classes
from rare_variant_enrichment.storage import MinimumDistanceStore


ENRICHMENT_HEADER = (
    "z_threshold",
    "tail",
    "annotation_family",
    "annotation_class",
    "ac_class",
    "ac_kind",
    "ac_value",
    "distance_bp",
    "total_observations",
    "outlier_observations",
    "nonoutlier_observations",
    "n11",
    "n10",
    "n01",
    "n00",
    "outlier_carrier_rate",
    "nonoutlier_carrier_rate",
    "carrier_rate_ratio",
    "odds_ratio",
    "odds_ratio_corrected_0_5",
    "fisher_p_value",
    "fisher_fdr_bh",
)


def fisher_exact_two_sided(n11: int, n10: int, n01: int, n00: int) -> float:
    cells = (n11, n10, n01, n00)
    if any(isinstance(cell, bool) or not isinstance(cell, int) or cell < 0 for cell in cells):
        raise ValueError("Fisher table cells must be non-negative integers")

    row1 = n11 + n10
    row2 = n01 + n00
    col1 = n11 + n01
    total = row1 + row2
    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    if lower == upper:
        return 1.0

    candidate_mode = min(upper, max(lower, ((row1 + 1) * (col1 + 1)) // (total + 2)))
    mode = candidate_mode
    observed_log_relative = _hypergeometric_log_relative(
        mode, n11, row1, col1, total
    )
    inclusion_threshold = observed_log_relative + math.log1p(1e-12)
    if 0.0 <= inclusion_threshold:
        return 1.0

    normalizer = 1.0
    normalizer_compensation = 0.0
    included = 0.0
    included_compensation = 0.0

    log_relative = 0.0
    for x in range(mode, lower, -1):
        log_relative += _hypergeometric_log_ratio_down(x, row1, col1, total)
        relative_probability = math.exp(log_relative)
        normalizer, normalizer_compensation = _compensated_add(
            normalizer, normalizer_compensation, relative_probability
        )
        if log_relative <= inclusion_threshold:
            included, included_compensation = _compensated_add(
                included, included_compensation, relative_probability
            )
        if relative_probability == 0.0:
            break

    log_relative = 0.0
    for x in range(mode, upper):
        log_relative += _hypergeometric_log_ratio_up(x, row1, col1, total)
        relative_probability = math.exp(log_relative)
        normalizer, normalizer_compensation = _compensated_add(
            normalizer, normalizer_compensation, relative_probability
        )
        if log_relative <= inclusion_threshold:
            included, included_compensation = _compensated_add(
                included, included_compensation, relative_probability
            )
        if relative_probability == 0.0:
            break

    p_value = included / normalizer
    return min(1.0, max(0.0, p_value))


def _hypergeometric_log_relative(
    mode: int, target: int, row1: int, col1: int, total: int
) -> float:
    relative = 0.0
    if target < mode:
        for x in range(mode, target, -1):
            relative += _hypergeometric_log_ratio_down(x, row1, col1, total)
            if relative < -800.0:
                return float("-inf")
    else:
        for x in range(mode, target):
            relative += _hypergeometric_log_ratio_up(x, row1, col1, total)
            if relative < -800.0:
                return float("-inf")
    return relative


def _hypergeometric_log_ratio_down(x: int, row1: int, col1: int, total: int) -> float:
    return (
        math.log(x)
        + math.log(total - col1 - row1 + x)
        - math.log(col1 - x + 1)
        - math.log(row1 - x + 1)
    )


def _hypergeometric_log_ratio_up(x: int, row1: int, col1: int, total: int) -> float:
    return (
        math.log(col1 - x)
        + math.log(row1 - x)
        - math.log(x + 1)
        - math.log(total - col1 - row1 + x + 1)
    )


def _compensated_add(total: float, compensation: float, value: float) -> tuple[float, float]:
    adjusted = value - compensation
    updated = total + adjusted
    return updated, (updated - total) - adjusted


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, max(0.0, running))
    return adjusted


def calculate_enrichment(
    phenotype_bed: Path,
    shared_samples_path: Path,
    carriers_path: Path,
    selected_features_path: Path,
    exact_ac: Sequence[int],
    cumulative_ac_max: Sequence[int],
    z_thresholds: Sequence[float],
    distance_thresholds: Sequence[int],
    tail: str,
    output_tsv: Path,
    output_json: Path,
    *,
    consequence_classes: Sequence[str] = (),
    loftee_enabled: bool = False,
    phenotype_qc_path: Path | None = None,
    chromosome_qc_path: Path | None = None,
    selected_chromosomes: Sequence[str] | None = None,
    container_image: str | None = None,
    workflow_version: str = "unknown",
    max_retries: int = 0,
    index_provenance: str = "unknown",
    vat_index_provenance: str = "unknown",
    maximum_gvs_maf: float = 0.01,
    annotation_chunk_size_bp: int = 10_000_000,
) -> None:
    thresholds = _validate_statistics_thresholds(z_thresholds)
    distances = _validate_distances(distance_thresholds)
    ac_classes = _validate_ac_classes(exact_ac, cumulative_ac_max)
    annotation_classes = build_annotation_classes(consequence_classes, loftee_enabled)
    if tail not in {"absolute", "positive", "negative"}:
        raise ValueError(f"Unsupported tail mode: {tail}")
    if isinstance(max_retries, bool) or not isinstance(max_retries, int) or max_retries < 0:
        raise ValueError("max_retries must be a non-negative integer")
    if index_provenance not in {"generated", "supplied", "unknown"}:
        raise ValueError("index_provenance must be generated, supplied, or unknown")

    shared_samples = read_nonempty_lines(shared_samples_path)
    if not shared_samples:
        raise ValueError("At least one shared sample is required")
    _reject_duplicates(shared_samples, "sample")
    shared_sample_set = set(shared_samples)
    selected_features = _read_selected_features(selected_features_path)
    if not selected_features:
        raise ValueError("At least one selected feature is required")

    class_indexes = {ac_class.label: index for index, ac_class in enumerate(ac_classes)}
    annotation_indexes = {
        (annotation.family, annotation.label): index
        for index, annotation in enumerate(annotation_classes)
    }
    carrier_counts = [
        [[0 for _ in distances] for _ in ac_classes]
        for _ in annotation_classes
    ]
    carrier_outlier_counts = [
        [
            [[0 for _ in distances] for _ in ac_classes]
            for _ in annotation_classes
        ]
        for _ in thresholds
    ]
    outlier_counts = [0 for _ in thresholds]
    total_observations = 0
    missing_z_observations = 0
    seen_bed_features: set[str] = set()
    seen_selected_features: set[str] = set()
    unselected_feature_count = 0

    with MinimumDistanceStore(output_tsv.parent) as carrier_minima:
        for (
            sample_id,
            feature_id,
            ac_label,
            annotation_family,
            annotation_class,
            distance,
        ) in _iter_carrier_file(carriers_path):
            if sample_id not in shared_sample_set:
                raise ValueError(f"Carrier sample is absent from shared samples: {sample_id}")
            if feature_id not in selected_features:
                raise ValueError(f"Carrier feature is absent from selected features: {feature_id}")
            if ac_label not in class_indexes:
                raise ValueError(f"Carrier AC class is not configured: {ac_label}")
            annotation_key = (annotation_family, annotation_class)
            if annotation_key not in annotation_indexes:
                raise ValueError(
                    "Carrier annotation class is not configured: "
                    f"{annotation_family}/{annotation_class}"
                )
            carrier_minima.upsert(
                sample_id,
                feature_id,
                ac_label,
                annotation_family,
                annotation_class,
                distance,
            )

        carrier_key_count = carrier_minima.count()
        features_with_carriers = carrier_minima.distinct_feature_count()
        with open_text(phenotype_bed) as bed_handle:
            header = _read_header(bed_handle)
            bed_samples = header[4:]
            _reject_duplicates(bed_samples, "sample")
            bed_sample_indexes = {
                sample_id: index for index, sample_id in enumerate(bed_samples)
            }
            missing_samples = [
                sample for sample in shared_samples if sample not in bed_sample_indexes
            ]
            if missing_samples:
                raise ValueError(
                    f"Shared samples are absent from phenotype BED: {', '.join(missing_samples)}"
                )
            shared_columns = [bed_sample_indexes[sample] for sample in shared_samples]
            shared_sample_indexes = {
                sample_id: index for index, sample_id in enumerate(shared_samples)
            }

            for line_number, raw_line in enumerate(bed_handle, start=2):
                if not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) != len(header):
                    raise ValueError(
                        f"Line {line_number} has {len(fields)} columns; expected {len(header)}"
                    )
                chrom, start_text, end_text, feature_id = fields[:4]
                _, tss = _parse_tss_interval(start_text, end_text, line_number)
                if not chrom:
                    raise ValueError(f"Line {line_number} has an empty chromosome")
                if not feature_id:
                    raise ValueError(f"Line {line_number} has an empty feature ID")
                if feature_id in seen_bed_features:
                    raise ValueError(f"Duplicate feature ID: {feature_id}")
                seen_bed_features.add(feature_id)

                expected_location = selected_features.get(feature_id)
                if expected_location is None:
                    unselected_feature_count += 1
                    continue
                if expected_location != (chrom, tss):
                    raise ValueError(
                        f"Selected feature location does not match phenotype BED: {feature_id}"
                    )
                seen_selected_features.add(feature_id)

                shared_values = [
                    _parse_phenotype_value(fields[4 + column_index], line_number)
                    for column_index in shared_columns
                ]
                for value in shared_values:
                    if value is None:
                        missing_z_observations += 1
                        continue
                    total_observations += 1
                    for threshold_index, threshold in enumerate(thresholds):
                        if classify_outlier(value, threshold, tail):
                            outlier_counts[threshold_index] += 1

                for (
                    sample_id,
                    ac_label,
                    annotation_family,
                    annotation_class,
                    minimum_distance,
                ) in carrier_minima.iter_feature(feature_id):
                    value = shared_values[shared_sample_indexes[sample_id]]
                    if value is None:
                        continue
                    class_index = class_indexes[ac_label]
                    annotation_index = annotation_indexes[(annotation_family, annotation_class)]
                    outlier_flags = [
                        classify_outlier(value, threshold, tail) for threshold in thresholds
                    ]
                    for distance_index, distance_threshold in enumerate(distances):
                        if minimum_distance > distance_threshold:
                            continue
                        carrier_counts[annotation_index][class_index][distance_index] += 1
                        for threshold_index, is_outlier in enumerate(outlier_flags):
                            if is_outlier:
                                carrier_outlier_counts[threshold_index][annotation_index][
                                    class_index
                                ][distance_index] += 1

    missing_selected_features = sorted(set(selected_features) - seen_selected_features)
    if missing_selected_features:
        raise ValueError(
            "Selected features are absent from phenotype BED: "
            + ", ".join(missing_selected_features)
        )

    rows = _build_rows(
        thresholds,
        distances,
        ac_classes,
        annotation_classes,
        tail,
        total_observations,
        outlier_counts,
        carrier_counts,
        carrier_outlier_counts,
    )
    _write_rows(output_tsv, rows)
    feature_count = len(selected_features)
    summary: dict[str, object] = {
            "ac_classes": [ac_class.label for ac_class in ac_classes],
            "analysis_parameters": {
                "cumulative_allele_count_maxima": list(cumulative_ac_max),
                "distance_thresholds_bp": distances,
                "exact_allele_counts": list(exact_ac),
                "outlier_tail": tail,
                "z_thresholds": thresholds,
            },
            "carrier_key_count": carrier_key_count,
            "distance_thresholds_bp": distances,
            "emitted_rows": len(rows),
            "feature_count": feature_count,
            "features_with_carriers": features_with_carriers,
            "features_without_carriers": feature_count - features_with_carriers,
            "missing_z_observations": missing_z_observations,
            "shared_sample_count": len(shared_samples),
            "statistical_limitation": (
                "Pooled Fisher tests are screening statistics because samples and features recur "
                "as observations; use a dependence-aware model for confirmatory inference."
            ),
            "tail": tail,
            "total_observations": total_observations,
            "tss_coordinate_convention": "TSS equals BED end; distance boundaries are inclusive.",
            "z_thresholds": thresholds,
            "unselected_feature_count": unselected_feature_count,
            "provenance": {
                "container_image": container_image,
                "max_retries": max_retries,
                "selected_chromosomes": list(
                    selected_chromosomes
                    if selected_chromosomes is not None
                    else dict.fromkeys(chrom for chrom, _ in selected_features.values())
                ),
                "software_versions": {
                    "python": platform.python_version(),
                    "rare_variant_enrichment": __version__,
                    "sqlite": sqlite3.sqlite_version,
                    "workflow": workflow_version,
                },
                "vcf_index": index_provenance,
            },
        }
    if phenotype_qc_path is not None:
        summary["phenotype_qc"] = _read_phenotype_qc(phenotype_qc_path)
    if chromosome_qc_path is not None:
        summary["chromosome_qc"] = _read_chromosome_qc(chromosome_qc_path)
    write_json(output_json, summary)


def _read_selected_features(path: Path) -> dict[str, tuple[str, int]]:
    selected: dict[str, tuple[str, int]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("Selected feature TSV is empty") from error
        if header != ["chrom", "tss", "feature_id"]:
            raise ValueError("Selected feature TSV header must be chrom, tss, feature_id")
        for line_number, fields in enumerate(reader, start=2):
            if len(fields) != 3:
                raise ValueError(
                    f"Selected feature TSV line {line_number} must have three columns"
                )
            chrom, tss_text, feature_id = fields
            try:
                tss = int(tss_text)
            except ValueError as error:
                raise ValueError(
                    f"Selected feature TSV line {line_number} has a non-integer TSS"
                ) from error
            if not chrom or not feature_id or tss < 1:
                raise ValueError(
                    f"Selected feature TSV line {line_number} has invalid feature values"
                )
            if feature_id in selected:
                raise ValueError(f"Duplicate selected feature ID: {feature_id}")
            selected[feature_id] = (chrom, tss)
    return selected


def _read_phenotype_qc(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Phenotype QC JSON is invalid: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError("Phenotype QC JSON must contain an object")
    identifier_fields = {"bed_only_samples", "vcf_only_samples", "shared_samples"}
    if identifier_fields.intersection(payload):
        raise ValueError("Phenotype QC must publish overlap counts, not sample IDs")
    return payload


def _read_chromosome_qc(path: Path) -> dict[str, object]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "chromosome" not in reader.fieldnames:
            raise ValueError("Chromosome QC TSV must contain a chromosome column")
        per_chromosome: list[dict[str, object]] = []
        totals: dict[str, int | float] = {}
        seen_chromosomes: set[str] = set()
        for row in reader:
            chromosome = row.get("chromosome", "")
            if not chromosome or chromosome in seen_chromosomes:
                raise ValueError("Chromosome QC TSV must contain unique chromosomes")
            seen_chromosomes.add(chromosome)
            parsed_row: dict[str, object] = {"chromosome": chromosome}
            for key, text in row.items():
                if key == "chromosome":
                    continue
                value = _parse_qc_scalar(text or "")
                parsed_row[key] = value
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    totals[key] = totals.get(key, 0) + value
            per_chromosome.append(parsed_row)
    return {"per_chromosome": per_chromosome, "totals": totals}


def _parse_qc_scalar(text: str) -> object:
    if text == "":
        return None
    if re.fullmatch(r"-?[0-9]+", text):
        return int(text)
    try:
        numeric = float(text)
    except ValueError:
        return text
    return numeric if math.isfinite(numeric) else text


def _validate_statistics_thresholds(z_thresholds: Sequence[float]) -> list[float]:
    thresholds = _validate_thresholds(z_thresholds)
    if any(threshold < 0 for threshold in thresholds):
        raise ValueError("z-score thresholds must be non-negative")
    if len(thresholds) != len(set(thresholds)):
        raise ValueError("z-score thresholds must be unique")
    return thresholds


def _validate_distances(distance_thresholds: Sequence[int]) -> list[int]:
    if not distance_thresholds:
        raise ValueError("At least one distance threshold is required")
    distances = list(distance_thresholds)
    if any(isinstance(distance, bool) or not isinstance(distance, int) or distance < 0 for distance in distances):
        raise ValueError("Distance thresholds must be non-negative integers")
    if len(distances) != len(set(distances)):
        raise ValueError("Distance thresholds must be unique")
    return distances


def _validate_ac_classes(
    exact_ac: Sequence[int], cumulative_ac_max: Sequence[int]
) -> list[AcClass]:
    exact = _validate_ac_values(exact_ac, "Exact AC values")
    cumulative = _validate_ac_values(cumulative_ac_max, "Cumulative AC maxima")
    if not exact and not cumulative:
        raise ValueError("At least one AC class is required")
    return build_ac_classes(exact, cumulative)


def _validate_ac_values(values: Sequence[int], label: str) -> list[int]:
    validated = list(values)
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 1 for value in validated):
        raise ValueError(f"{label} must be positive integers")
    if len(validated) != len(set(validated)):
        raise ValueError(f"{label} must be unique")
    return validated


def _build_rows(
    thresholds: Sequence[float],
    distances: Sequence[int],
    ac_classes: Sequence[AcClass],
    annotation_classes: Sequence[AnnotationClass],
    tail: str,
    total_observations: int,
    outlier_counts: Sequence[int],
    carrier_counts: Sequence[Sequence[Sequence[int]]],
    carrier_outlier_counts: Sequence[Sequence[Sequence[Sequence[int]]]],
) -> list[list[str]]:
    rows: list[list[str]] = []
    p_values: list[float] = []
    for threshold_index, threshold in enumerate(thresholds):
        outlier_observations = outlier_counts[threshold_index]
        nonoutlier_observations = total_observations - outlier_observations
        for annotation_index, annotation in enumerate(annotation_classes):
            for class_index, ac_class in enumerate(ac_classes):
                for distance_index, distance in enumerate(distances):
                    n11 = carrier_outlier_counts[threshold_index][annotation_index][class_index][
                        distance_index
                    ]
                    n10 = outlier_observations - n11
                    n01 = carrier_counts[annotation_index][class_index][distance_index] - n11
                    n00 = nonoutlier_observations - n01
                    outlier_rate = _divide_or_none(n11, outlier_observations)
                    nonoutlier_rate = _divide_or_none(n01, nonoutlier_observations)
                    rate_ratio = (
                        None
                        if outlier_rate is None or nonoutlier_rate in {None, 0.0}
                        else outlier_rate / nonoutlier_rate
                    )
                    odds_ratio = _divide_or_none(n11 * n00, n10 * n01)
                    corrected_odds_ratio = (
                        (n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5))
                    )
                    fisher_p = fisher_exact_two_sided(n11, n10, n01, n00)
                    p_values.append(fisher_p)
                    rows.append(
                        [
                            str(threshold),
                            tail,
                            annotation.family,
                            annotation.label,
                            ac_class.label,
                            ac_class.kind,
                            str(ac_class.value),
                            str(distance),
                            str(total_observations),
                            str(outlier_observations),
                            str(nonoutlier_observations),
                            str(n11),
                            str(n10),
                            str(n01),
                            str(n00),
                            _format_statistic(outlier_rate),
                            _format_statistic(nonoutlier_rate),
                            _format_statistic(rate_ratio),
                            _format_statistic(odds_ratio),
                            _format_statistic(corrected_odds_ratio),
                            str(fisher_p),
                        ]
                    )

    for row, adjusted_p in zip(rows, benjamini_hochberg(p_values), strict=True):
        row.append(str(adjusted_p))
    return rows


def _divide_or_none(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else numerator / denominator


def _format_statistic(value: float | None) -> str:
    return "NA" if value is None else str(value)


def _write_rows(path: Path, rows: Sequence[Sequence[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(ENRICHMENT_HEADER)
        writer.writerows(rows)
