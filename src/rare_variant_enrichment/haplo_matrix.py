"""Gene-by-sample haplo calls from covariate-adjusted log2-CPM."""

import csv
import gzip
import logging
import math
from pathlib import Path
from typing import Sequence

import numpy as np


LOGGER = logging.getLogger(__name__)


def validate_haplo_drop(drop: float) -> None:
    if not math.isfinite(drop) or drop < 0:
        raise ValueError('Haplo logCPM drop must be finite and non-negative')


def _centered_predictor_basis(predictors: np.ndarray) -> tuple[np.ndarray | None, str | None]:
    """An intercept is removed by centering; retain at least one residual df."""
    n, k = predictors.shape
    if n <= k + 1:
        return None, 'insufficient_dof'
    try:
        design = np.column_stack([np.ones(n), predictors])
        if np.linalg.matrix_rank(design) < k + 1:
            return None, 'rank_deficiency'
        centered = predictors - predictors.mean(axis=0)
        basis = np.linalg.qr(centered, mode='reduced')[0] if k else np.empty((n, 0))
        return basis, None
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None, 'numerical_failure'


def export_haplo_matrix(
    expression: Sequence[tuple[str, np.ndarray]],
    pc_values: np.ndarray,
    pc_count: int,
    sample_ids: Sequence[str],
    output: Path,
    drop: float,
    *,
    additional_covariates: np.ndarray | None = None,
) -> dict:
    """Match the haplo rule without dividing residuals by gene-specific SD.

    Adjusted expression minus its gene mean equals the residual from a
    fit with an intercept, all additional covariates, and selected phenotype PCs.
    The caller supplies the same aligned genes and samples as the Z export.
    """
    validate_haplo_drop(drop)
    pcs = pc_values[:, :pc_count]
    covariates = (
        np.empty((len(sample_ids), 0), dtype=float) if additional_covariates is None
        else np.asarray(additional_covariates, dtype=float)
    )
    if covariates.ndim != 2 or covariates.shape[0] != len(sample_ids):
        raise ValueError('Expression and additional covariates have incompatible shapes')
    predictors = np.column_stack([covariates, pcs])
    predictor_usable = np.all(np.isfinite(predictors), axis=1)
    # Reuse the full basis for complete genes; refit incomplete genes on usable rows.
    common_basis, common_reason = _centered_predictor_basis(predictors[predictor_usable])
    counts = {'positive_cell_count': 0, 'negative_cell_count': 0, 'missing_cell_count': 0}
    exclusions = {'insufficient_dof': 0, 'rank_deficiency': 0, 'numerical_failure': 0}
    LOGGER.info(
        'Starting haplo matrix export: selected_pc_count=%d additional_covariate_count=%d logcpm_drop=%g',
        pc_count, covariates.shape[1], drop,
    )
    with gzip.open(output, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_id', *sample_ids])
        for gene_id, values in expression:
            usable = predictor_usable & np.isfinite(values)
            basis, reason = (
                (common_basis, common_reason) if np.array_equal(usable, predictor_usable)
                else _centered_predictor_basis(predictors[usable])
            )
            calls = np.full(len(sample_ids), 'NA', dtype=object)
            if reason is None:
                observations = values[usable]
                # Identical values have exactly zero drop, including at drop=0.
                # Avoid mean-rounding noise for decimal constants such as 0.1.
                centered = (
                    np.zeros_like(observations) if np.all(observations == observations[0])
                    else observations - observations.mean()
                )
                residuals = centered - basis @ (basis.T @ centered)
                residuals -= residuals.mean()
                if np.all(np.isfinite(residuals)):
                    calls[usable] = np.where(residuals < -drop, '1', '0')
                else:
                    reason = 'numerical_failure'
            if reason is not None:
                exclusions[reason] += 1
            counts['positive_cell_count'] += int(np.count_nonzero(calls == '1'))
            counts['negative_cell_count'] += int(np.count_nonzero(calls == '0'))
            counts['missing_cell_count'] += int(np.count_nonzero(calls == 'NA'))
            writer.writerow([gene_id, *calls])
    LOGGER.info('Completed haplo matrix export: %s; positive_cells=%d', output, counts['positive_cell_count'])
    return {
        'input_scale': 'log2-CPM',
        'logcpm_drop': drop,
        'selected_pc_count': pc_count,
        'design': 'intercept plus all additional covariates plus first k phenotype principal components',
        'additional_covariates_used': covariates.shape[1] > 0,
        'additional_covariate_count': covariates.shape[1],
        'rule': 'adjusted_log2_cpm < gene_mean_adjusted_log2_cpm - logcpm_drop',
        'gene_mean_samples': 'finite observations in the aligned export cohort',
        'missing_value': 'NA',
        'constant_gene_rule': 'zero drop; non-outlier when the fit is valid',
        'exclusion_counts': exclusions,
        **counts,
    }
