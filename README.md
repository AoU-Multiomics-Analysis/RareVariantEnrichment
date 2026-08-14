# Rare variant enrichment

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml)
[![Python tests](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This workflow screens for enrichment of molecular low-expression outliers among pooled loss-of-function (LoF) carrier observations. It directly residualizes every protein-coding gene's expression against an intercept plus principal components (PCs), then pools eligible sample–gene observations into Fisher exact 2×2 tables. It does not accept variant call files, annotation tables, genomic index files, or chromosome selections.

## Inputs

The public WDL has four required file inputs.

### Molecular phenotype BED

`phenotype_bed` is a wide BED, plain or gzip-compressed, with exactly the first three headers `#chr`, `start`, and `end`, followed by a non-empty feature-ID header and sample IDs:

```text
#chr  start  end  gene_id       SAMPLE_1  SAMPLE_2  ...
chr1  999    1000 ENSG000001.2 -1.8      NA        ...
```

Each interval must be one base wide. Feature IDs are normalized by removing exactly one terminal `.<digits>` suffix, so `ENSG000001.2` becomes `ENSG000001` and `ENSG000001.2.3` becomes `ENSG000001.2`. Empty values, `.`, `NA`, and `NaN` are missing; non-finite numeric values are rejected.

### LoF carrier table

`lof_carrier_table` is a tab-separated, plain or gzip-compressed table with these required columns:

```text
sample_id  gene_id  gene_symbol  has_lof_variant  n_lof_variants  variant_ids  lof_classes
```

`sample_id` and `gene_id` are required for every row. `has_lof_variant` accepts case-insensitive `true`, `false`, `1`, `0`, `yes`, or `no`. Truthy rows collapse at normalized `(sample_id, gene_id)` level into three carrier definitions:

- `any_lof`: every truthy row.
- `HC`: a truthy row whose comma-separated `lof_classes` contains `HC`.
- `HC_or_LC`: a truthy row whose classes contain `HC` or `LC`.

Classes are trimmed and case-normalized. Unknown or missing classes still contribute to `any_lof`, but not the class-specific definitions.

### Principal components

`principal_components_tsv` is a TSV with a strict consecutive PC header:

```text
ID  PC1  PC2  PC3  ...
SAMPLE_1  -0.14  1.02  0.37  ...
```

`ID` values are unique, non-empty strings. All PC values must be finite, including columns beyond the largest selected PC count. Header-only files and nonconsecutive labels such as `PC1,PC3` are rejected. The analysis uses the string-ID intersection of BED and PC samples.

### Gene annotation GTF

`gene_annotation_gtf` is a plain or gzip-compressed GTF. Only nine-field `gene` records with `gene_type "protein_coding"` enter analysis. The first task writes a sorted normalized `gene_id` list and QC JSON; it fails if no coding genes are found.

## PC grid and thresholds

`negative_z_thresholds` defaults to `[-2.0, -3.0, -4.0, -5.0, -6.0]`. Values must be finite, unique, and strictly negative. An outlier is a residual z-score `<=` its threshold.

`pc_counts` defaults to `[]`, which selects an adaptive grid: every count 0–10, every 10 through 100, every 50 through 500, every 100 above 500, and always the final available PC count. Provide an explicit strictly increasing, unique list of non-negative counts no greater than the available PC count to override it.

`pc_counts_per_job` defaults to `10`. It controls how many selected PC-count settings are evaluated by each analysis job; it does not change the number of PCA columns available to a model. The workflow validates the PC header, partitions the selected grid into jobs of at most this size, and merges their intermediate outputs before publishing the final analysis files.

For each PC count and coding gene, the workflow fits finite expression observations to an intercept plus the first `k` PCs. It requires observations greater than fitted rank plus one, rejects rank-deficient designs and zero/non-finite residual SD, centers residuals, and divides by population SD (`ddof=0`). Missing BED values remain missing rather than becoming residuals.

## Run

Build or choose a container image that includes this package and NumPy 2 or later. The WDL defaults are intentionally high for cohort-scale input localization: preparation uses 2 CPUs, 32 GB RAM, and 500 GB disk; analysis, PC-grid chunk preparation, and shard merging use 8 CPUs, 128 GB RAM, and the 1000 GB analysis-disk baseline. Disk is dynamically raised to `ceil(2 × localized input GiB + 20)` when required: from the GTF for preparation, the PC file for chunk preparation, analysis inputs for each shard, and all localized shard outputs for merging. `max_retries` defaults to 1 for every task.

```bash
miniwdl run workflows/rare_variant_enrichment.wdl \
  -i examples/rare_variant_enrichment.inputs.json
```

Update the four file paths in the example. For reproducibility, override `docker_image` with an immutable image digest in production.

The Python CLI exposes the same two public operations:

```bash
rare-variant-enrichment prepare-protein-coding-genes \
  --gtf genes.gtf.gz --genes-output protein_coding_genes.tsv --qc-output protein_coding_genes.qc.json

rare-variant-enrichment lof-pc-enrichment \
  --phenotype-bed phenotypes.bed.gz --lof-carriers lof_carriers.tsv \
  --principal-components pcs.tsv --protein-coding-genes protein_coding_genes.tsv \
  --negative-z-thresholds=-2,-3,-4,-5,-6 --pc-counts '' \
  --results-output lof_pc_enrichment.tsv --summary-output lof_pc_enrichment.summary.json \
  --gene-pc-qc-output lof_pc_enrichment.gene_pc_qc.tsv.gz \
  --analysis-qc-output lof_pc_enrichment.analysis_qc.json
```

Legacy Python CLI commands remain available for compatibility, but are not part of the WDL or public analysis contract.

## Outputs

The workflow emits exactly six files:

- `results_tsv`: one row per PC count × negative threshold × carrier definition, merged across analysis shards.
- `summary_json`: selected grid/settings, global FDR scope, residualization description, provenance, and the screening limitation.
- `gene_pc_qc_tsv_gz`: compressed per-normalized-gene/per-PC QC with usable samples, rank, residual mean/SD, status, and exclusion reason.
- `analysis_qc_json`: BED/PC overlap counts, pre-join carrier-pair counts, LoF input QC, and per-PC eligibility, actual carrier-observation, and structured exclusion counters.
- `protein_coding_genes_tsv`: the prepared sorted coding-gene list.
- `protein_coding_genes_qc_json`: GTF record and normalization QC.

`results_tsv` contains these analysis-ready raw cells and derived statistics:

| Column | Meaning |
|---|---|
| `pc_count`, `z_threshold`, `carrier_definition` | The tested residualization and carrier stratum. |
| `eligible_gene_count`, `total_observations`, `outlier_observations`, `carrier_observations` | Pooled denominator counts after PC-specific gene exclusions. |
| `n11`, `n10`, `n01`, `n00` | Outlier carrier, nonoutlier carrier, outlier noncarrier, and nonoutlier noncarrier cells. |
| `outlier_carrier_rate`, `nonoutlier_carrier_rate`, `carrier_rate_ratio` | Carrier prevalence by outlier state and their ratio. |
| `odds_ratio`, `odds_ratio_corrected_0_5` | Uncorrected and 0.5-cell-corrected odds ratios. |
| `fisher_p_value`, `fisher_fdr_bh` | Two-sided Fisher exact p-value and global Benjamini–Hochberg FDR across every emitted row. |

The `analysis_qc_json` reconciles each PC-specific `carrier_observations` count with every result row: `carrier_observations = n11 + n10`. Pre-join carrier counts are intentionally reported separately because carriers absent from the BED/PC-gene observation set do not enter a Fisher table.

Per-shard analysis files are workflow intermediates, not public outputs. During merge, `fisher_fdr_bh` is recomputed across every final result row, preserving a global FDR scope rather than retaining shard-local adjustments.

## Interpretation

Each test pools repeated samples and repeated genes, so Fisher p-values and globally adjusted FDR values are screening statistics, not confirmatory person-level inference. Use the raw cells, PC-specific QC, and appropriate dependence-aware models, permutations, or gene-level meta-analysis for confirmation. The workflow deliberately does not export full residual matrices.
