import csv
from dataclasses import dataclass
import gzip
import logging
import math
import platform
from pathlib import Path
import re
from typing import Mapping, Sequence

import numpy as np

from rare_variant_enrichment import __version__
from rare_variant_enrichment.io import open_text, write_json
from rare_variant_enrichment.phenotypes import (
    _parse_phenotype_value,
    _parse_tss_interval,
    _read_header,
    _reject_duplicates,
)
from rare_variant_enrichment.statistics import (
    benjamini_hochberg,
    fisher_exact_two_sided,
)


CARRIER_DEFINITIONS = ("any_lof", "HC", "HC_or_LC")
EXCLUSION_REASONS = (
    "insufficient_dof",
    "rank_deficiency",
    "invalid_or_zero_residual_sd",
    "other",
)
LOF_REQUIRED_COLUMNS = (
    "sample_id",
    "gene_id",
    "gene_symbol",
    "has_lof_variant",
    "n_lof_variants",
    "variant_ids",
    "lof_classes",
)
RESULT_HEADER = (
    "pc_count",
    "z_threshold",
    "carrier_definition",
    "eligible_gene_count",
    "total_observations",
    "outlier_observations",
    "carrier_observations",
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
GENE_PC_QC_HEADER = (
    "gene_id",
    "pc_count",
    "usable_sample_count",
    "rank",
    "residual_mean",
    "residual_sd",
    "status",
    "exclusion_reason",
)
PROGRESS_INTERVAL_GENES = 500
LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PrincipalComponents:
    sample_ids: tuple[str, ...]
    values: np.ndarray

    @property
    def available_pc_count(self) -> int:
        return int(self.values.shape[1])


@dataclass(frozen=True)
class LofCarriers:
    pairs_by_definition: dict[str, set[tuple[str, str]]]
    qc: dict[str, int]


@dataclass(frozen=True)
class ResidualFit:
    z_scores: np.ndarray
    usable_sample_count: int
    rank: int | None
    residual_mean: float | None
    residual_sd: float | None
    exclusion_reason: str | None


@dataclass(frozen=True)
class CompleteDataProjection:
    """Cached projection terms for complete-data expression residualization."""

    centered_expression: np.ndarray
    orthonormal_pcs: np.ndarray
    component_scores: np.ndarray

    def advance_prediction(
        self, previous_pc_count: int, pc_count: int, prediction: np.ndarray
    ) -> np.ndarray:
        available_pc_count = self.orthonormal_pcs.shape[1]
        if (
            isinstance(previous_pc_count, bool)
            or not isinstance(previous_pc_count, int)
            or isinstance(pc_count, bool)
            or not isinstance(pc_count, int)
            or previous_pc_count < 0
            or pc_count < previous_pc_count
            or pc_count > available_pc_count
        ):
            raise ValueError("PC counts are outside the prepared principal components")
        values = np.asarray(prediction, dtype=float)
        if values.shape != self.centered_expression.shape:
            raise ValueError("Prediction has incompatible expression shape")
        return values + (
            self.orthonormal_pcs[:, previous_pc_count:pc_count]
            @ self.component_scores[previous_pc_count:pc_count]
        )

    def z_scores(self, prediction: np.ndarray) -> np.ndarray:
        values = np.asarray(prediction, dtype=float)
        if values.shape != self.centered_expression.shape:
            raise ValueError("Prediction has incompatible expression shape")
        residuals = self.centered_expression - values
        residual_sd = np.std(residuals, axis=0, ddof=0)
        z_scores = np.full(residuals.shape, np.nan, dtype=float)
        valid_columns = np.isfinite(residual_sd) & (residual_sd > 0)
        z_scores[:, valid_columns] = (
            residuals[:, valid_columns] / residual_sd[valid_columns]
        )
        return z_scores


def prepare_complete_data_projection(
    expression: np.ndarray,
    principal_components: np.ndarray,
    requested_pc_counts: Sequence[int],
) -> CompleteDataProjection:
    """Prepare complete-data residualization through centered PC projections."""

    expression_values = np.asarray(expression, dtype=float)
    pc_values = np.asarray(principal_components, dtype=float)
    if expression_values.ndim != 2 or pc_values.ndim != 2:
        raise ValueError("Expression and principal components must be two-dimensional")
    if expression_values.shape[0] != pc_values.shape[0]:
        raise ValueError("Expression and principal components have incompatible shapes")
    if not np.all(np.isfinite(expression_values)):
        raise ValueError("Complete-data projection requires complete finite expression")
    if not np.all(np.isfinite(pc_values)):
        raise ValueError("Complete-data projection requires finite principal components")

    available_pc_count = pc_values.shape[1]
    validated_counts: list[int] = []
    for pc_count in requested_pc_counts:
        if (
            isinstance(pc_count, bool)
            or not isinstance(pc_count, int)
            or pc_count < 0
            or pc_count > available_pc_count
        ):
            raise ValueError("PC count is outside the available principal components")
        validated_counts.append(pc_count)
    if not validated_counts:
        raise ValueError("At least one PC count is required")

    maximum_pc_count = max(validated_counts)
    centered_expression = expression_values - np.mean(expression_values, axis=0)
    centered_pcs = pc_values[:, :maximum_pc_count] - np.mean(
        pc_values[:, :maximum_pc_count], axis=0
    )
    pc_norms = np.linalg.norm(centered_pcs, axis=0)
    if np.any(~np.isfinite(pc_norms)) or np.any(pc_norms == 0):
        raise ValueError("Principal components must have nonzero finite variance")
    orthonormal_pcs, _ = np.linalg.qr(centered_pcs, mode="reduced")
    component_scores = orthonormal_pcs.T @ centered_expression
    return CompleteDataProjection(
        centered_expression=centered_expression,
        orthonormal_pcs=orthonormal_pcs,
        component_scores=component_scores,
    )


def normalize_ensembl_id(value: str) -> str:
    return re.sub(r"\.[0-9]+$", "", value, count=1)


def prepare_protein_coding_genes(
    gtf_path: Path, genes_output: Path, qc_output: Path
) -> None:
    genes: set[str] = set()
    total_lines = 0
    gene_records = 0
    coding_records = 0
    duplicate_records = 0

    with open_text(gtf_path) as handle:
        for raw_line in handle:
            total_lines += 1
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            gene_records += 1
            attributes = _parse_gtf_attributes(fields[8])
            if attributes.get("gene_type") != "protein_coding":
                continue
            gene_id = normalize_ensembl_id(attributes.get("gene_id", "").strip())
            if not gene_id:
                continue
            coding_records += 1
            if gene_id in genes:
                duplicate_records += 1
            genes.add(gene_id)

    if not genes:
        raise ValueError("No protein-coding genes found in GTF")

    with genes_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id"])
        writer.writerows((gene_id,) for gene_id in sorted(genes))
    write_json(
        qc_output,
        {
            "duplicate_normalized_gene_records": duplicate_records,
            "nine_field_gene_records": gene_records,
            "protein_coding_gene_records": coding_records,
            "total_lines": total_lines,
            "unique_protein_coding_gene_count": len(genes),
        },
    )


def _parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for field in text.split(";"):
        field = field.strip()
        if not field:
            continue
        pieces = field.split(None, 1)
        if len(pieces) != 2:
            continue
        key, value = pieces
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] == '"':
            value = value[1:-1]
        attributes[key] = value
    return attributes


def read_principal_components(path: Path) -> PrincipalComponents:
    sample_ids: list[str] = []
    rows: list[list[float]] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("Principal-components TSV is empty") from error
        if len(header) < 2 or header[0] != "ID":
            raise ValueError("Principal-components TSV header must begin with ID, PC1")
        expected_pcs = [f"PC{index}" for index in range(1, len(header))]
        if header[1:] != expected_pcs:
            raise ValueError(
                "Principal-components TSV header PCs must be consecutive from PC1"
            )
        for line_number, fields in enumerate(reader, start=2):
            if not fields or all(not field.strip() for field in fields):
                continue
            if len(fields) != len(header):
                raise ValueError(
                    f"Principal-components TSV line {line_number} has {len(fields)} "
                    f"columns; expected {len(header)}"
                )
            sample_id = fields[0].strip()
            if not sample_id:
                raise ValueError(
                    f"Principal-components TSV line {line_number} has an empty sample ID"
                )
            if sample_id in seen:
                raise ValueError(f"Duplicate PC sample ID: {sample_id}")
            seen.add(sample_id)
            values: list[float] = []
            for value_text in fields[1:]:
                try:
                    value = float(value_text)
                except ValueError as error:
                    raise ValueError(
                        f"Principal-components TSV line {line_number} has a non-numeric PC value"
                    ) from error
                if not math.isfinite(value):
                    raise ValueError("All principal-component values must be finite")
                values.append(value)
            sample_ids.append(sample_id)
            rows.append(values)
    if not sample_ids:
        raise ValueError("Principal-components TSV must contain at least one sample")
    return PrincipalComponents(tuple(sample_ids), np.asarray(rows, dtype=float))


def build_pc_grid(requested: Sequence[int], available: int) -> list[int]:
    if isinstance(available, bool) or not isinstance(available, int) or available < 0:
        raise ValueError("Available PC count must be a non-negative integer")
    if requested:
        values = list(requested)
        if any(
            isinstance(value, bool)
            or not isinstance(value, int)
            or value < 0
            or value > available
            for value in values
        ):
            raise ValueError(
                "PC counts must be non-negative integers no greater than available PCs"
            )
        if any(left >= right for left, right in zip(values, values[1:])):
            raise ValueError("PC counts must be strictly increasing and unique")
        return values

    values = list(range(0, min(available, 10) + 1))
    values.extend(range(20, min(available, 100) + 1, 10))
    values.extend(range(150, min(available, 500) + 1, 50))
    values.extend(range(600, available + 1, 100))
    if not values or values[-1] != available:
        values.append(available)
    return values


def read_lof_carriers(path: Path) -> LofCarriers:
    pairs = {definition: set() for definition in CARRIER_DEFINITIONS}
    input_rows = 0
    truthy_rows = 0
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("LoF carrier table is empty") from error
        if len(header) != len(set(header)):
            raise ValueError("LoF carrier table header contains duplicate columns")
        missing = [column for column in LOF_REQUIRED_COLUMNS if column not in header]
        if missing:
            raise ValueError("Missing required LoF carrier columns: " + ", ".join(missing))
        indexes = {column: header.index(column) for column in LOF_REQUIRED_COLUMNS}
        for line_number, fields in enumerate(reader, start=2):
            if not fields or all(not field.strip() for field in fields):
                continue
            input_rows += 1
            if len(fields) != len(header):
                raise ValueError(
                    f"LoF carrier table line {line_number} has {len(fields)} columns; "
                    f"expected {len(header)}"
                )
            sample_id = fields[indexes["sample_id"]].strip()
            gene_id = normalize_ensembl_id(fields[indexes["gene_id"]].strip())
            if not sample_id:
                raise ValueError(f"LoF carrier table line {line_number} has an empty sample ID")
            if not gene_id:
                raise ValueError(f"LoF carrier table line {line_number} has an empty gene ID")
            if not _parse_truth_value(
                fields[indexes["has_lof_variant"]], line_number
            ):
                continue
            truthy_rows += 1
            pair = (sample_id, gene_id)
            pairs["any_lof"].add(pair)
            classes = {
                token.strip().upper()
                for token in fields[indexes["lof_classes"]].split(",")
                if token.strip()
            }
            if "HC" in classes:
                pairs["HC"].add(pair)
                pairs["HC_or_LC"].add(pair)
            elif "LC" in classes:
                pairs["HC_or_LC"].add(pair)

    return LofCarriers(
        pairs,
        {
            "input_row_count": input_rows,
            "truthy_row_count": truthy_rows,
            "unique_any_lof_pair_count": len(pairs["any_lof"]),
            "unique_hc_or_lc_pair_count": len(pairs["HC_or_LC"]),
            "unique_hc_pair_count": len(pairs["HC"]),
        },
    )


def _parse_truth_value(value: str, line_number: int) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(
        f"LoF carrier table line {line_number} has invalid has_lof_variant; "
        "expected true/false/1/0/yes/no"
    )


def validate_negative_z_thresholds(thresholds: Sequence[float]) -> list[float]:
    values: list[float] = []
    for threshold in thresholds:
        if isinstance(threshold, bool):
            raise ValueError("Expected finite, unique, negative z-score thresholds")
        try:
            value = float(threshold)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "Expected finite, unique, negative z-score thresholds"
            ) from error
        if not math.isfinite(value) or value >= 0:
            raise ValueError("Expected finite, unique, negative z-score thresholds")
        values.append(value)
    if not values or len(values) != len(set(values)):
        raise ValueError("Expected finite, unique, negative z-score thresholds")
    return values


def residualize_expression(
    expression: np.ndarray, principal_components: np.ndarray, pc_count: int
) -> ResidualFit:
    values = np.asarray(expression, dtype=float)
    pcs = np.asarray(principal_components, dtype=float)
    if values.ndim != 1 or pcs.ndim != 2 or pcs.shape[0] != values.shape[0]:
        raise ValueError("Expression and principal components have incompatible shapes")
    if (
        isinstance(pc_count, bool)
        or not isinstance(pc_count, int)
        or pc_count < 0
        or pc_count > pcs.shape[1]
    ):
        raise ValueError("PC count is outside the available principal components")

    z_scores = np.full(values.shape, np.nan, dtype=float)
    usable = np.isfinite(values)
    if pc_count:
        usable &= np.all(np.isfinite(pcs[:, :pc_count]), axis=1)
    usable_count = int(np.count_nonzero(usable))
    selected_y = values[usable]
    selected_pcs = pcs[usable, :pc_count]
    design = np.column_stack((np.ones(usable_count), selected_pcs))

    try:
        rank = int(np.linalg.matrix_rank(design))
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return ResidualFit(z_scores, usable_count, None, None, None, "other")
    if usable_count <= rank + 1:
        return ResidualFit(
            z_scores, usable_count, rank, None, None, "insufficient_dof"
        )
    if rank < pc_count + 1:
        return ResidualFit(
            z_scores, usable_count, rank, None, None, "rank_deficiency"
        )

    try:
        coefficients = np.linalg.lstsq(design, selected_y, rcond=None)[0]
        residuals = selected_y - design @ coefficients
        residuals = residuals - np.mean(residuals)
        residual_mean = float(np.mean(residuals))
        residual_sd = float(np.std(residuals, ddof=0))
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return ResidualFit(z_scores, usable_count, rank, None, None, "other")

    scale = max(1.0, float(np.max(np.abs(selected_y))))
    numerical_zero = np.finfo(float).eps * scale * 16.0
    if not math.isfinite(residual_sd) or residual_sd <= numerical_zero:
        return ResidualFit(
            z_scores,
            usable_count,
            rank,
            residual_mean,
            residual_sd,
            "invalid_or_zero_residual_sd",
        )
    z_scores[usable] = residuals / residual_sd
    if not np.all(np.isfinite(z_scores[usable])):
        return ResidualFit(
            np.full(values.shape, np.nan, dtype=float),
            usable_count,
            rank,
            residual_mean,
            residual_sd,
            "invalid_or_zero_residual_sd",
        )
    return ResidualFit(
        z_scores, usable_count, rank, residual_mean, residual_sd, None
    )


def calculate_lof_pc_enrichment(
    phenotype_bed: Path,
    lof_carriers_path: Path,
    principal_components_path: Path,
    protein_coding_genes_path: Path,
    negative_z_thresholds: Sequence[float],
    requested_pc_counts: Sequence[int],
    results_output: Path,
    summary_output: Path,
    gene_pc_qc_output: Path,
    analysis_qc_output: Path,
) -> None:
    thresholds = validate_negative_z_thresholds(negative_z_thresholds)
    principal_components = read_principal_components(principal_components_path)
    pc_counts = build_pc_grid(
        requested_pc_counts, principal_components.available_pc_count
    )
    coding_genes = _read_protein_coding_genes(protein_coding_genes_path)
    carriers = read_lof_carriers(lof_carriers_path)
    LOGGER.info(
        "Starting LoF/PC enrichment: thresholds=%s pc_counts=%s coding_genes=%d",
        thresholds,
        pc_counts,
        len(coding_genes),
    )
    pc_sample_indexes = {
        sample_id: index
        for index, sample_id in enumerate(principal_components.sample_ids)
    }

    per_pc = {
        str(pc_count): {
            "carrier_observations": {
                definition: 0 for definition in CARRIER_DEFINITIONS
            },
            "eligible_gene_count": 0,
            "exclusion_counts": {reason: 0 for reason in EXCLUSION_REASONS},
            "total_observations": 0,
        }
        for pc_count in pc_counts
    }
    outlier_observations = {
        pc_count: {threshold: 0 for threshold in thresholds}
        for pc_count in pc_counts
    }
    outlier_carriers = {
        pc_count: {
            threshold: {definition: 0 for definition in CARRIER_DEFINITIONS}
            for threshold in thresholds
        }
        for pc_count in pc_counts
    }
    bed_gene_count = 0
    coding_bed_gene_count = 0
    seen_genes: set[str] = set()

    coding_expression: list[tuple[str, np.ndarray]] = []
    with open_text(phenotype_bed) as bed_handle:
        header = _read_header(bed_handle)
        bed_samples = [sample.strip() for sample in header[4:]]
        if any(not sample for sample in bed_samples):
            raise ValueError("Phenotype BED sample IDs must be non-empty strings")
        _reject_duplicates(bed_samples, "sample")
        shared_samples = [sample for sample in bed_samples if sample in pc_sample_indexes]
        if not shared_samples:
            raise ValueError("No shared phenotype BED and PC samples")
        LOGGER.info(
            "Found %d shared BED/PC samples (%d BED, %d PC)",
            len(shared_samples),
            len(bed_samples),
            len(principal_components.sample_ids),
        )
        shared_bed_columns = [bed_samples.index(sample) for sample in shared_samples]
        shared_pc_rows = [pc_sample_indexes[sample] for sample in shared_samples]
        shared_pc_values = principal_components.values[shared_pc_rows, :]

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
            gene_id = normalize_ensembl_id(feature_id.strip())
            if not gene_id:
                raise ValueError(f"Line {line_number} has an empty feature ID")
            if gene_id in seen_genes:
                raise ValueError(f"Duplicate normalized phenotype gene ID: {gene_id}")
            seen_genes.add(gene_id)
            bed_gene_count += 1
            values = [
                _parse_phenotype_value(value, line_number) for value in fields[4:]
            ]
            if gene_id not in coding_genes:
                continue
            coding_bed_gene_count += 1
            if coding_bed_gene_count % PROGRESS_INTERVAL_GENES == 0:
                LOGGER.info(
                    "Processed %d protein-coding BED genes",
                    coding_bed_gene_count,
                )
            coding_expression.append(
                (
                    gene_id,
                    np.asarray(
                        [
                            np.nan if values[column] is None else values[column]
                            for column in shared_bed_columns
                        ],
                        dtype=float,
                    ),
                )
            )

    if coding_bed_gene_count == 0:
        raise ValueError(
            "No protein-coding genes from the supplied list occur in the phenotype BED"
        )

    with gzip.open(
        gene_pc_qc_output, "wt", encoding="utf-8", newline=""
    ) as qc_handle:
        qc_writer = csv.writer(qc_handle, delimiter="\t", lineterminator="\n")
        qc_writer.writerow(GENE_PC_QC_HEADER)

        def record_fit(gene_id: str, pc_count: int, fit: ResidualFit) -> None:
            rank_text = "NA" if fit.rank is None else str(fit.rank)
            mean_text = _format_optional_float(fit.residual_mean)
            sd_text = _format_optional_float(fit.residual_sd)
            if fit.exclusion_reason is not None:
                reason = (
                    fit.exclusion_reason
                    if fit.exclusion_reason in EXCLUSION_REASONS
                    else "other"
                )
                per_pc[str(pc_count)]["exclusion_counts"][reason] += 1
                qc_writer.writerow(
                    [
                        gene_id,
                        pc_count,
                        fit.usable_sample_count,
                        rank_text,
                        mean_text,
                        sd_text,
                        "excluded",
                        reason,
                    ]
                )
                return

            qc_writer.writerow(
                [
                    gene_id,
                    pc_count,
                    fit.usable_sample_count,
                    rank_text,
                    mean_text,
                    sd_text,
                    "included",
                    "",
                ]
            )
            pc_qc = per_pc[str(pc_count)]
            pc_qc["eligible_gene_count"] += 1
            pc_qc["total_observations"] += fit.usable_sample_count
            finite = np.isfinite(fit.z_scores)
            carrier_masks: dict[str, np.ndarray] = {}
            for definition in CARRIER_DEFINITIONS:
                mask = np.asarray(
                    [
                        finite[index]
                        and (sample_id, gene_id)
                        in carriers.pairs_by_definition[definition]
                        for index, sample_id in enumerate(shared_samples)
                    ],
                    dtype=bool,
                )
                carrier_masks[definition] = mask
                pc_qc["carrier_observations"][definition] += int(
                    np.count_nonzero(mask)
                )
            for threshold in thresholds:
                outlier_mask = finite & (fit.z_scores <= threshold)
                outlier_observations[pc_count][threshold] += int(
                    np.count_nonzero(outlier_mask)
                )
                for definition in CARRIER_DEFINITIONS:
                    outlier_carriers[pc_count][threshold][definition] += int(
                        np.count_nonzero(outlier_mask & carrier_masks[definition])
                    )

        complete_expression = all(
            np.all(np.isfinite(expression)) for _, expression in coding_expression
        )
        rank_deficient_requested_prefix = False
        if complete_expression:
            for pc_count in pc_counts:
                try:
                    rank = int(
                        np.linalg.matrix_rank(
                            np.column_stack(
                                (
                                    np.ones(len(shared_samples)),
                                    shared_pc_values[:, :pc_count],
                                )
                            )
                        )
                    )
                except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                    rank_deficient_requested_prefix = True
                    break
                if rank < pc_count + 1:
                    rank_deficient_requested_prefix = True
                    break

        if complete_expression and not rank_deficient_requested_prefix:
            expression_matrix = np.column_stack(
                [expression for _, expression in coding_expression]
            )
            projection = prepare_complete_data_projection(
                expression_matrix, shared_pc_values, pc_counts
            )
            prediction = np.zeros_like(expression_matrix)
            previous_pc_count = 0

            for pc_count in pc_counts:
                prediction = projection.advance_prediction(
                    previous_pc_count, pc_count, prediction
                )
                z_scores = projection.z_scores(prediction)
                residuals = projection.centered_expression - prediction
                usable_sample_count = len(shared_samples)
                try:
                    rank = int(
                        np.linalg.matrix_rank(
                            np.column_stack(
                                (
                                    np.ones(usable_sample_count),
                                    shared_pc_values[:, :pc_count],
                                )
                            )
                        )
                    )
                except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                    rank = None

                for gene_index, (gene_id, expression) in enumerate(coding_expression):
                    if rank is None:
                        fit = ResidualFit(
                            np.full(usable_sample_count, np.nan),
                            usable_sample_count,
                            None,
                            None,
                            None,
                            "other",
                        )
                    elif usable_sample_count <= rank + 1:
                        fit = ResidualFit(
                            np.full(usable_sample_count, np.nan),
                            usable_sample_count,
                            rank,
                            None,
                            None,
                            "insufficient_dof",
                        )
                    elif rank < pc_count + 1:
                        fit = ResidualFit(
                            np.full(usable_sample_count, np.nan),
                            usable_sample_count,
                            rank,
                            None,
                            None,
                            "rank_deficiency",
                        )
                    else:
                        residual_mean = float(np.mean(residuals[:, gene_index]))
                        residual_sd = float(np.std(residuals[:, gene_index], ddof=0))
                        numerical_zero = (
                            np.finfo(float).eps
                            * max(1.0, float(np.max(np.abs(expression))))
                            * 16.0
                        )
                        exclusion_reason = None
                        if (
                            not math.isfinite(residual_sd)
                            or residual_sd <= numerical_zero
                            or not np.all(np.isfinite(z_scores[:, gene_index]))
                        ):
                            exclusion_reason = "invalid_or_zero_residual_sd"
                        fit = ResidualFit(
                            z_scores[:, gene_index],
                            usable_sample_count,
                            rank,
                            residual_mean,
                            residual_sd,
                            exclusion_reason,
                        )
                    record_fit(gene_id, pc_count, fit)

                pc_qc = per_pc[str(pc_count)]
                LOGGER.info(
                    "Completed PC count %d: eligible_genes=%d observations=%d",
                    pc_count,
                    pc_qc["eligible_gene_count"],
                    pc_qc["total_observations"],
                )
                previous_pc_count = pc_count
        else:
            for pc_count in pc_counts:
                for gene_id, expression in coding_expression:
                    record_fit(
                        gene_id,
                        pc_count,
                        residualize_expression(
                            expression, shared_pc_values, pc_count
                        ),
                    )

                pc_qc = per_pc[str(pc_count)]
                LOGGER.info(
                    "Completed PC count %d: eligible_genes=%d observations=%d",
                    pc_count,
                    pc_qc["eligible_gene_count"],
                    pc_qc["total_observations"],
                )

    rows = _build_result_rows(
        pc_counts, thresholds, per_pc, outlier_observations, outlier_carriers
    )
    _write_result_rows(results_output, rows)
    write_json(
        analysis_qc_output,
        {
            "bed_gene_count": bed_gene_count,
            "bed_sample_count": len(bed_samples),
            "pc_sample_count": len(principal_components.sample_ids),
            "protein_coding_bed_gene_count": coding_bed_gene_count,
            "protein_coding_gene_count": len(coding_genes),
            "shared_bed_pc_sample_count": len(shared_samples),
            "pre_join_carrier_pair_counts": {
                definition: len(carriers.pairs_by_definition[definition])
                for definition in CARRIER_DEFINITIONS
            },
            "lof_carrier_table": carriers.qc,
            "per_pc": per_pc,
        },
    )
    write_json(
        summary_output,
        {
            "available_pc_count": principal_components.available_pc_count,
            "carrier_definitions": list(CARRIER_DEFINITIONS),
            "emitted_result_rows": len(rows),
            "fdr_scope": "global_across_all_emitted_rows",
            "negative_z_thresholds": thresholds,
            "observation_unit": "eligible sample-gene residual",
            "pc_grid_mode": "adaptive" if not requested_pc_counts else "explicit",
            "provenance": {
                "input_files": {
                    "lof_carriers": str(lof_carriers_path),
                    "phenotype_bed": str(phenotype_bed),
                    "principal_components": str(principal_components_path),
                    "protein_coding_genes": str(protein_coding_genes_path),
                },
                "software_versions": {
                    "numpy": np.__version__,
                    "python": platform.python_version(),
                    "rare_variant_enrichment": __version__,
                },
            },
            "residualization": {
                "design": "intercept plus first k principal components",
                "outlier_rule": "residual_z <= z_threshold",
                "residual_standard_deviation_ddof": 0,
            },
            "selected_pc_counts": pc_counts,
            "statistical_limitation": (
                "Pooled Fisher tests are screening statistics because samples and genes recur "
                "as observations; use dependence-aware inference for confirmation."
            ),
        },
    )
    LOGGER.info(
        "Wrote LoF/PC enrichment outputs: results=%s summary=%s gene_pc_qc=%s analysis_qc=%s",
        results_output,
        summary_output,
        gene_pc_qc_output,
        analysis_qc_output,
    )


def _read_protein_coding_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError("Protein-coding gene TSV is empty") from error
        if header != ["gene_id"]:
            raise ValueError("Protein-coding gene TSV header must be gene_id")
        for line_number, fields in enumerate(reader, start=2):
            if not fields or all(not field.strip() for field in fields):
                continue
            if len(fields) != 1:
                raise ValueError(
                    f"Protein-coding gene TSV line {line_number} must have one column"
                )
            gene_id = normalize_ensembl_id(fields[0].strip())
            if not gene_id:
                raise ValueError(
                    f"Protein-coding gene TSV line {line_number} has an empty gene ID"
                )
            if gene_id in genes:
                raise ValueError(f"Duplicate normalized protein-coding gene ID: {gene_id}")
            genes.add(gene_id)
    if not genes:
        raise ValueError("Protein-coding gene TSV must contain at least one gene")
    return genes


def _build_result_rows(
    pc_counts: Sequence[int],
    thresholds: Sequence[float],
    per_pc: Mapping[str, Mapping[str, object]],
    outlier_observations: Mapping[int, Mapping[float, int]],
    outlier_carriers: Mapping[int, Mapping[float, Mapping[str, int]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for pc_count in pc_counts:
        pc_qc = per_pc[str(pc_count)]
        total = int(pc_qc["total_observations"])
        carrier_counts = pc_qc["carrier_observations"]
        for threshold in thresholds:
            outliers = outlier_observations[pc_count][threshold]
            for definition in CARRIER_DEFINITIONS:
                carriers = int(carrier_counts[definition])
                n11 = outlier_carriers[pc_count][threshold][definition]
                n10 = carriers - n11
                n01 = outliers - n11
                n00 = total - n11 - n10 - n01
                outlier_rate = _divide_or_none(n11, n11 + n01)
                nonoutlier_rate = _divide_or_none(n10, n10 + n00)
                rate_ratio = (
                    None
                    if outlier_rate is None or nonoutlier_rate in {None, 0.0}
                    else outlier_rate / nonoutlier_rate
                )
                odds_ratio = _divide_or_none(n11 * n00, n10 * n01)
                corrected = (n11 + 0.5) * (n00 + 0.5) / (
                    (n10 + 0.5) * (n01 + 0.5)
                )
                p_value = fisher_exact_two_sided(n11, n10, n01, n00)
                p_values.append(p_value)
                rows.append(
                    {
                        "pc_count": pc_count,
                        "z_threshold": threshold,
                        "carrier_definition": definition,
                        "eligible_gene_count": pc_qc["eligible_gene_count"],
                        "total_observations": total,
                        "outlier_observations": outliers,
                        "carrier_observations": carriers,
                        "n11": n11,
                        "n10": n10,
                        "n01": n01,
                        "n00": n00,
                        "outlier_carrier_rate": _format_optional_float(outlier_rate),
                        "nonoutlier_carrier_rate": _format_optional_float(nonoutlier_rate),
                        "carrier_rate_ratio": _format_optional_float(rate_ratio),
                        "odds_ratio": _format_optional_float(odds_ratio),
                        "odds_ratio_corrected_0_5": _format_optional_float(corrected),
                        "fisher_p_value": _format_optional_float(p_value),
                        "fisher_fdr_bh": "",
                    }
                )
    for row, adjusted in zip(rows, benjamini_hochberg(p_values)):
        row["fisher_fdr_bh"] = _format_optional_float(adjusted)
    return rows


def _write_result_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=RESULT_HEADER, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _divide_or_none(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else numerator / denominator


def _format_optional_float(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return str(value)
