# Rare variant enrichment

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml)
[![Python tests](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This workflow screens for enrichment of molecular low-expression outliers among pooled loss-of-function (LoF) carrier observations. It directly residualizes every protein-coding gene's expression against an intercept, optional additional covariates, and principal components (PCs), then pools eligible sample–gene observations into Fisher exact 2×2 tables. It does not accept variant call files, annotation tables, genomic index files, or chromosome selections.

## Standalone variant carrier extraction

`workflows/extract_variant_carriers.wdl` creates gene-matched carrier tables from
a filtered VCF and the MTtoVCF `TranscriptAnnotations` table. It does not run an
enrichment test. Run the active enrichment workflow separately after you select
the carrier class that you need.

The input VCF must already have the required quality and frequency filters. The
extractor does not apply a quality, AC, AF, or MAF filter. It joins VCF alleles
to transcript rows by exact chromosome, position, REF, and ALT. It then groups
the matched rows by the normalized Ensembl gene ID. It does not use a GTF, a
phenotype position, a gene-distance window, or a feature list.

For each allele and gene, the extractor retains:

- the most severe Ensembl consequence;
- all distinct consequence terms;
- LoFTEE `HC` in preference to `LC`;
- the maximum finite transcript-level `revel` score;
- the maximum valid `gvs_max_af` value.

The initial classes are `lof_hc`, `lof_hc_or_lc`, `missense`, `splice_core`,
and `splice_region`. The detailed audit also keeps variants that do not have an
initial class. You can define a missense subset with a REVEL threshold after
extraction. This does not require a new VCF or transcript extraction run.

Required inputs are `filtered_vcf`, its `.tbi` index, and block-gzipped
`transcript_annotations`. The transcript index is optional. The preparation
task validates a supplied index or creates one. The default scatter includes
`chr1` through `chr22` and uses 10 Mb non-overlapping annotation chunks.

```bash
miniwdl run workflows/extract_variant_carriers.wdl \
  -i examples/extract_variant_carriers.inputs.json
```

The workflow returns these principal tables:

| Output | Row definition |
|---|---|
| `variant_carrier_audit_tsv_gz` | One sample, exact allele, and normalized VAT gene. |
| `variant_carriers_tsv_gz` | One sample, gene, and initial variant class. |

The audit includes ALT dosage, VCF AC and AF, the collapsed consequence fields,
LoFTEE, REVEL, `gvs_max_af`, and the assigned initial classes. The carrier table
includes `sample_id`, `gene_id`, `gene_symbol`, `variant_class`, `n_variants`,
and sorted `variant_ids`.

`variant_carriers_qc_json` embeds preparation and chromosome QC. It records
input and index provenance, row totals, duplicate counts, unique sample, gene,
allele, and allele-gene counts, and class counts. QC files contain aggregate
counts only. They do not contain sample-level audit records.

## Inputs

The public WDL has four required file inputs and one optional file input (`additional_covariates_tsv`).

### Molecular phenotype BED

`phenotype_bed` is a wide BED, plain or gzip-compressed, with exactly the first three headers `#chr`, `start`, and `end`, followed by a non-empty feature-ID header and sample IDs:

```text
#chr  start  end  gene_id       SAMPLE_1  SAMPLE_2  ...
chr1  999    1000 ENSG000001.2 -1.8      NA        ...
```

Each interval must be one base wide. Feature IDs are normalized by extracting exactly one Ensembl gene token and removing its numeric version suffix. This supports the molecular-phenotype IDs used by the Susie merged files:

- expression: `ENSG00000000419.14` → `ENSG00000000419`;
- proteomics: `A0JNW5_ENSG00000111647.13` → `ENSG00000111647`;
- splicing: `chr20:50941209:50942031:clu_63027_-:ENSG00000000419.14` → `ENSG00000000419`.

Empty, unsupported, or ambiguous IDs are rejected with the input line number. If multiple phenotype rows normalize to the same gene, the analysis collapses them to one gene-level vector by retaining the minimum finite z-score independently for each sample. This represents the most extreme negative splice outlier for that gene; missing values are ignored, and a value remains missing only when all rows for that gene are missing. The downstream residualization and enrichment tests remain gene-level. Empty values, `.`, `NA`, and `NaN` are missing; non-finite numeric values are rejected.

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

### Additional covariates

`additional_covariates_tsv` is an optional tab-separated matrix for covariates such as genotype PCs. It must contain a sample-ID column named `sample_id` or `ID` (in any column position); every other column is included as a finite numeric covariate. Sample IDs are treated as strings, so identifiers such as `001` are preserved. The analysis uses the intersection of BED, molecular-PC, and covariate samples, preserving the BED sample order. The sparse LoF carrier table is not included in this intersection because an absent carrier row represents a noncarrier.

### Gene annotation GTF

`gene_annotation_gtf` is a plain or gzip-compressed GTF. Only nine-field `gene` records with `gene_type "protein_coding"` enter analysis. The first task writes a sorted normalized `gene_id` list and QC JSON; it fails if no coding genes are found.

## PC grid and thresholds

`negative_z_thresholds` defaults to `[-2.0, -3.0, -4.0, -5.0, -6.0]`. Values must be finite, unique, and strictly negative. An outlier is a residual z-score `<=` its threshold.

`pc_counts` defaults to `[]`, which selects an adaptive grid: every count 0–10, every 10 through 100, every 50 through 500, and every 100 above 500. Above 500, the adaptive grid rounds the permitted maximum down to the nearest 100 and stops 100 PCs earlier, with a lower limit of 500. For example, a permitted maximum of 5219 stops at 5100. If the permitted maximum is 500 or less, the adaptive grid includes it. If the PC matrix contains `n` PCs, an explicit strictly increasing, unique list can include non-negative counts no greater than `n-1`.

`pc_counts_per_job` defaults to `10`. It controls how many selected PC-count settings are evaluated by each analysis job; it does not change the number of PCA columns available to a model. The workflow validates the PC header, partitions the selected grid into jobs of at most this size, and merges their intermediate outputs before publishing the final analysis files.

For each PC count and coding gene, the workflow fits finite expression observations to an intercept, all supplied additional covariates, and the first `k` molecular PCs. It requires observations greater than fitted rank plus one, rejects rank-deficient designs and zero/non-finite residual SD, centers residuals, and divides by population SD (`ddof=0`). Missing BED values remain missing rather than becoming residuals. When no additional covariate matrix is supplied, the model reduces to the legacy intercept-plus-molecular-PC model.

## Run

Build or choose a container image that includes this package and NumPy 2 or later. The default image also includes R, `data.table`, `ggplot2`, and `ggrepel` for the PC-sweep QC plot. The WDL defaults are intentionally high for cohort-scale input localization: preparation uses 2 CPUs, 32 GB RAM, and 500 GB disk; analysis shards and shard merging use 8 CPUs, 128 GB RAM, and the 1000 GB analysis-disk baseline; header-only PC-grid chunk preparation uses 1 CPU and 4 GB RAM. Disk is dynamically raised to `ceil(2 × localized input GiB + 20)` when required: from the GTF for preparation, the PC file for chunk preparation, analysis inputs (including optional covariates) for each shard, and all localized shard outputs for merging. `pc_preemptible` defaults to 2 and controls the preemptible retry count for the PC-enrichment scatter jobs; `max_retries` defaults to 1 for every task.

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
  --additional-covariates genetic_pcs.tsv \
  --negative-z-thresholds=-2,-3,-4,-5,-6 --pc-counts '' \
  --results-output lof_pc_enrichment.tsv --summary-output lof_pc_enrichment.summary.json \
  --gene-pc-qc-output lof_pc_enrichment.gene_pc_qc.tsv.gz \
  --analysis-qc-output lof_pc_enrichment.analysis_qc.json
```

Legacy Python CLI commands remain available for compatibility, but are not part of the WDL or public analysis contract.

## Outputs

The workflow emits exactly ten files:

- `results_tsv`: one row per PC count × negative threshold × carrier definition, merged across analysis shards.
- `summary_json`: selected grid/settings, global FDR scope, residualization description, provenance, and the screening limitation.
- `gene_pc_qc_tsv_gz`: compressed per-normalized-gene/per-PC QC with usable samples, rank, residual mean/SD, status, and exclusion reason.
- `analysis_qc_json`: BED/PC/covariate overlap counts, covariate names and input sample count when supplied, pre-join carrier-pair counts, LoF input QC, and per-PC eligibility, actual carrier-observation, and structured exclusion counters.
- `pc_selection_json`: median-logOR plateau summaries, excluded/included z thresholds, and the minimum common PC count selected by the 95% plateau rule.
- `enrichment_plot_svg`: threshold-specific enrichment curves for `HC` and `any_lof`, median logOR curves, and reference lines for the selected PC positions.
- `pc_sweep_qc_summary_tsv`: analysis-ready PC-sweep values containing the median log odds ratio across the selected z thresholds, the maximum-enrichment PC and odds ratio, and each PC's percentage of the maximum.
- `pc_sweep_qc_plot_png`: percentage-of-maximum QC plot with exact odds-ratio annotations at selected PC checkpoints and ggrepel-style labels for the maximum enrichment values.
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

The workflow also emits `pc_selection_json` and `enrichment_plot_svg`. PC selection summarizes the enrichment curves after excluding `z = -2` by default, because that threshold is not intended to represent the true outlier set. For each of `HC` and `any_lof`, it computes the median log odds ratio across `z = -3, -4, -5, -6`, finds the maximum median log odds ratio, and identifies the earliest PC count reaching 95% of that maximum. The reported selected PC count is the larger of those two definition-specific plateau-entry counts, so both carrier definitions meet the criterion while the number of PCs remains minimal. The SVG retains each threshold curve, overlays the median log-odds curve, and marks both definition-specific plateau entries and the common selected PC count.

The `pc_sweep_qc_summary_tsv` and `pc_sweep_qc_plot_png` outputs use the same `z = -3, -4, -5, -6` median-logOR summary. Each carrier definition is normalized to its own maximum median odds ratio; the plot displays percentage of maximum on the y-axis, a 95% plateau reference, and exact median odds-ratio labels at selected PC counts.

## Interpretation

Each test pools repeated samples and repeated genes, so Fisher p-values and globally adjusted FDR values are screening statistics, not confirmatory person-level inference. Use the raw cells, PC-specific QC, and appropriate dependence-aware models, permutations, or gene-level meta-analysis for confirmation. The workflow deliberately does not export full residual matrices.
