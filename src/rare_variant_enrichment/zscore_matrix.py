"""Export gene-by-sample residual Z scores at the selected PC count."""

import csv
import gzip
import json
import logging
import os
from pathlib import Path
import re

import numpy as np

from rare_variant_enrichment.io import write_json
from rare_variant_enrichment.lof_pc import (
    EXCLUSION_REASONS,
    GENE_PC_QC_HEADER,
    PROGRESS_INTERVAL_GENES,
    _format_optional_float,
    build_pc_grid,
    iter_gene_residual_fits,
    read_additional_covariates,
    read_aligned_gene_expression,
    read_principal_components,
)


LOGGER = logging.getLogger(__name__)


def export_zscore_matrix(
    phenotype_bed: Path,
    principal_components_path: Path,
    selection_input: Path,
    matrix_output: Path,
    gene_qc_output: Path,
    summary_output: Path,
    *,
    additional_covariates_path: Path | None = None,
) -> None:
    inputs = [phenotype_bed, principal_components_path, selection_input]
    if additional_covariates_path is not None:
        inputs.append(additional_covariates_path)
    for path in inputs:
        # Path collapses the double slash in a URI. Detect either spelling.
        if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:/", str(path)):
            raise ValueError(f"Input localization error: unresolved cloud URI: {path}")
        if not path.is_file() or not os.access(path, os.R_OK):
            raise ValueError(f"Input localization error: file is not readable: {path}")

    with selection_input.open() as handle:
        payload = json.load(handle)
    selection = payload.get("selection") if isinstance(payload, dict) else None
    pc_count = selection.get("selected_pc_count") if isinstance(selection, dict) else None
    if isinstance(pc_count, bool) or not isinstance(pc_count, int) or pc_count < 0:
        raise ValueError("Selection must contain a nonnegative integer selected_pc_count")
    pcs = read_principal_components(principal_components_path)
    build_pc_grid([pc_count], pcs.available_pc_count)
    covariates = (
        read_additional_covariates(additional_covariates_path)
        if additional_covariates_path is not None else None
    )
    aligned = read_aligned_gene_expression(phenotype_bed, pcs, covariates)
    if not aligned.expression:
        raise ValueError("No genes occur in the phenotype BED")

    LOGGER.info(
        "Exporting Z-score matrix: selected_pc_count=%d genes=%d samples=%d",
        pc_count, len(aligned.expression), len(aligned.shared_samples),
    )
    exclusions = {reason: 0 for reason in EXCLUSION_REASONS}
    included = 0
    with (
        gzip.open(matrix_output, "wt", encoding="utf-8", newline="") as matrix_handle,
        gzip.open(gene_qc_output, "wt", encoding="utf-8", newline="") as qc_handle,
    ):
        matrix_writer = csv.writer(matrix_handle, delimiter="\t", lineterminator="\n")
        qc_writer = csv.writer(qc_handle, delimiter="\t", lineterminator="\n")
        matrix_writer.writerow(["gene_id", *aligned.shared_samples])
        qc_writer.writerow(GENE_PC_QC_HEADER)
        fits = iter_gene_residual_fits(
            aligned.expression, aligned.pc_values, [pc_count], aligned.covariate_values
        )
        for index, (gene_id, _, fit) in enumerate(fits, start=1):
            reason = fit.exclusion_reason
            if reason is None:
                included += 1
            else:
                exclusions[reason if reason in exclusions else "other"] += 1
            matrix_writer.writerow([
                gene_id,
                *(
                    format(float(value), ".17g")
                    if reason is None and np.isfinite(value) else "NA"
                    for value in fit.z_scores
                ),
            ])
            qc_writer.writerow([
                gene_id, pc_count, fit.usable_sample_count,
                "NA" if fit.rank is None else fit.rank,
                _format_optional_float(fit.residual_mean),
                _format_optional_float(fit.residual_sd),
                "included" if reason is None else "excluded", reason or "",
            ])
            if index % PROGRESS_INTERVAL_GENES == 0:
                LOGGER.info("Wrote Z scores for %d genes", index)

    write_json(summary_output, {
        "selected_pc_count": pc_count,
        "gene_count": len(aligned.expression),
        "included_gene_count": included,
        "excluded_gene_count": sum(exclusions.values()),
        "exclusion_counts": exclusions,
        "sample_count": len(aligned.shared_samples),
        "bed_sample_count": len(aligned.bed_samples),
        "pc_sample_count": len(pcs.sample_ids),
        "bed_feature_count": aligned.bed_feature_count,
        "duplicate_feature_count": aligned.bed_feature_count - aligned.bed_gene_count,
        "gene_scope": "all genes in phenotype BED",
        "gene_order": "first occurrence in phenotype BED",
        "sample_order": "BED order within BED/PC/additional-covariate intersection",
        "duplicate_gene_rule": "minimum finite value per sample before residualization",
        "additional_covariate_names": [] if covariates is None else list(covariates.names),
        "design": "intercept plus all additional covariates plus first k principal components",
        "residual_standard_deviation_ddof": 0,
        "missing_value": "NA",
    })
    LOGGER.info(
        "Wrote Z-score matrix: %s; included_genes=%d excluded_genes=%d",
        matrix_output, included, sum(exclusions.values()),
    )
