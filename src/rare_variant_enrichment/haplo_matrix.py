"""Gene-by-sample haplo calls from covariate-adjusted log2-CPM."""

import csv
import gzip
import logging
import math
from pathlib import Path
from typing import Sequence

import numpy as np

from rare_variant_enrichment.lof_pc import EXCLUSION_REASONS, iter_gene_residual_fits


OUTLIER_Z_THRESHOLD = -3.0


LOGGER = logging.getLogger(__name__)


def validate_haplo_drop(drop: float) -> None:
    if not math.isfinite(drop) or drop < 0:
        raise ValueError('Haplo logCPM drop must be finite and non-negative')


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
    """Require Z <= -3 and a log2-CPM residual strictly below -drop.

    Use the same fits, fixed covariates, and exclusion rules as the Z export.
    Compare the centered residual directly to avoid rounding from Z times SD.
    """
    validate_haplo_drop(drop)
    covariates = (
        np.empty((len(sample_ids), 0), dtype=float) if additional_covariates is None
        else np.asarray(additional_covariates, dtype=float)
    )
    if covariates.ndim != 2 or covariates.shape[0] != len(sample_ids):
        raise ValueError('Expression and additional covariates have incompatible shapes')
    counts = {'positive_cell_count': 0, 'negative_cell_count': 0, 'missing_cell_count': 0}
    exclusions = {reason: 0 for reason in EXCLUSION_REASONS}
    LOGGER.info(
        'Starting haplo matrix export: selected_pc_count=%d additional_covariate_count=%d logcpm_drop=%g',
        pc_count, covariates.shape[1], drop,
    )
    with gzip.open(output, 'wt', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['gene_id', *sample_ids])
        fits = iter_gene_residual_fits(expression, pc_values, [pc_count], covariates)
        for gene_id, _, fit in fits:
            calls = np.full(len(sample_ids), 'NA', dtype=object)
            if fit.exclusion_reason is None:
                usable = np.isfinite(fit.z_scores)
                residuals = fit.centered_residuals[usable]
                calls[usable] = np.where(
                    (fit.z_scores[usable] <= OUTLIER_Z_THRESHOLD) & (residuals < -drop),
                    '1', '0',
                )
            else:
                exclusions[fit.exclusion_reason if fit.exclusion_reason in exclusions else 'other'] += 1
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
        'z_threshold': OUTLIER_Z_THRESHOLD,
        'z_comparison': '<=',
        'rule': 'Z <= -3 and adjusted_log2_cpm < gene_mean_adjusted_log2_cpm - logcpm_drop',
        'gene_mean_samples': 'finite observations in the aligned export cohort',
        'missing_value': 'NA',
        'constant_gene_rule': 'NA because the residual Z score is undefined',
        'exclusion_counts': exclusions,
        **counts,
    }
