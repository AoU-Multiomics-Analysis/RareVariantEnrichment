import csv
import math
from pathlib import Path
from typing import Sequence

from rare_variant_enrichment.aggregation import _read_carrier_file
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


ENRICHMENT_HEADER = (
    "z_threshold",
    "tail",
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

    mode = min(upper, max(lower, ((row1 + 1) * (col1 + 1)) // (total + 2)))
    relative_probabilities = [0.0] * (upper - lower + 1)
    relative_probabilities[mode - lower] = 1.0

    for x in range(mode, upper):
        relative_probabilities[x + 1 - lower] = relative_probabilities[x - lower] * (
            (col1 - x) * (row1 - x)
            / ((x + 1) * (total - col1 - row1 + x + 1))
        )
    for x in range(mode, lower, -1):
        relative_probabilities[x - 1 - lower] = relative_probabilities[x - lower] * (
            x * (total - col1 - row1 + x)
            / ((col1 - x + 1) * (row1 - x + 1))
        )

    observed = relative_probabilities[n11 - lower]
    normalizer = math.fsum(relative_probabilities)
    p_value = math.fsum(
        probability
        for probability in relative_probabilities
        if probability <= observed * (1.0 + 1e-12)
    ) / normalizer
    return min(1.0, max(0.0, p_value))


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
    exact_ac: Sequence[int],
    cumulative_ac_max: Sequence[int],
    z_thresholds: Sequence[float],
    distance_thresholds: Sequence[int],
    tail: str,
    output_tsv: Path,
    output_json: Path,
) -> None:
    thresholds = _validate_statistics_thresholds(z_thresholds)
    distances = _validate_distances(distance_thresholds)
    ac_classes = _validate_ac_classes(exact_ac, cumulative_ac_max)
    if tail not in {"absolute", "positive", "negative"}:
        raise ValueError(f"Unsupported tail mode: {tail}")

    shared_samples = read_nonempty_lines(shared_samples_path)
    if not shared_samples:
        raise ValueError("At least one shared sample is required")
    _reject_duplicates(shared_samples, "sample")
    shared_sample_set = set(shared_samples)

    carrier_minima: dict[tuple[str, str, str], int] = {}
    _read_carrier_file(carriers_path, carrier_minima)
    class_indexes = {ac_class.label: index for index, ac_class in enumerate(ac_classes)}
    carriers_by_feature: dict[str, dict[str, dict[int, int]]] = {}
    for (sample_id, feature_id, ac_label), distance in carrier_minima.items():
        if sample_id not in shared_sample_set:
            raise ValueError(f"Carrier sample is absent from shared samples: {sample_id}")
        if ac_label not in class_indexes:
            raise ValueError(f"Carrier AC class is not configured: {ac_label}")
        feature_carriers = carriers_by_feature.setdefault(feature_id, {})
        sample_carriers = feature_carriers.setdefault(sample_id, {})
        sample_carriers[class_indexes[ac_label]] = distance

    carrier_counts = [
        [0 for _ in distances]
        for _ in ac_classes
    ]
    carrier_outlier_counts = [
        [[0 for _ in distances] for _ in ac_classes]
        for _ in thresholds
    ]
    outlier_counts = [0 for _ in thresholds]
    total_observations = 0
    missing_z_observations = 0
    seen_features: set[str] = set()

    with open_text(phenotype_bed) as bed_handle:
        header = _read_header(bed_handle)
        bed_samples = header[4:]
        _reject_duplicates(bed_samples, "sample")
        bed_sample_indexes = {sample_id: index for index, sample_id in enumerate(bed_samples)}
        missing_samples = [sample for sample in shared_samples if sample not in bed_sample_indexes]
        if missing_samples:
            raise ValueError(
                f"Shared samples are absent from phenotype BED: {', '.join(missing_samples)}"
            )
        shared_columns = [bed_sample_indexes[sample] for sample in shared_samples]

        for line_number, raw_line in enumerate(bed_handle, start=2):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"Line {line_number} has {len(fields)} columns; expected {len(header)}"
                )
            chrom, start_text, end_text, feature_id = fields[:4]
            _parse_tss_interval(start_text, end_text, line_number)
            if not chrom:
                raise ValueError(f"Line {line_number} has an empty chromosome")
            if not feature_id:
                raise ValueError(f"Line {line_number} has an empty feature ID")
            if feature_id in seen_features:
                raise ValueError(f"Duplicate feature ID: {feature_id}")
            seen_features.add(feature_id)

            values = [_parse_phenotype_value(value, line_number) for value in fields[4:]]
            feature_carriers = carriers_by_feature.get(feature_id, {})
            for sample_id, column_index in zip(shared_samples, shared_columns, strict=True):
                value = values[column_index]
                if value is None:
                    missing_z_observations += 1
                    continue

                total_observations += 1
                outlier_flags = [
                    classify_outlier(value, threshold, tail) for threshold in thresholds
                ]
                for threshold_index, is_outlier in enumerate(outlier_flags):
                    if is_outlier:
                        outlier_counts[threshold_index] += 1

                for class_index, minimum_distance in feature_carriers.get(sample_id, {}).items():
                    for distance_index, distance_threshold in enumerate(distances):
                        if minimum_distance > distance_threshold:
                            continue
                        carrier_counts[class_index][distance_index] += 1
                        for threshold_index, is_outlier in enumerate(outlier_flags):
                            if is_outlier:
                                carrier_outlier_counts[threshold_index][class_index][distance_index] += 1

    missing_carrier_features = sorted(set(carriers_by_feature) - seen_features)
    if missing_carrier_features:
        raise ValueError(
            "Carrier features are absent from phenotype BED: " + ", ".join(missing_carrier_features)
        )

    rows = _build_rows(
        thresholds,
        distances,
        ac_classes,
        tail,
        total_observations,
        outlier_counts,
        carrier_counts,
        carrier_outlier_counts,
    )
    _write_rows(output_tsv, rows)
    feature_count = len(seen_features)
    features_with_carriers = len(set(carriers_by_feature))
    write_json(
        output_json,
        {
            "ac_classes": [ac_class.label for ac_class in ac_classes],
            "carrier_key_count": len(carrier_minima),
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
        },
    )


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
    tail: str,
    total_observations: int,
    outlier_counts: Sequence[int],
    carrier_counts: Sequence[Sequence[int]],
    carrier_outlier_counts: Sequence[Sequence[Sequence[int]]],
) -> list[list[str]]:
    rows: list[list[str]] = []
    p_values: list[float] = []
    for threshold_index, threshold in enumerate(thresholds):
        outlier_observations = outlier_counts[threshold_index]
        nonoutlier_observations = total_observations - outlier_observations
        for class_index, ac_class in enumerate(ac_classes):
            for distance_index, distance in enumerate(distances):
                n11 = carrier_outlier_counts[threshold_index][class_index][distance_index]
                n10 = outlier_observations - n11
                n01 = carrier_counts[class_index][distance_index] - n11
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
