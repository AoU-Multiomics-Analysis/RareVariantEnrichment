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

rare-variant-enrichment analyze-lof-pc-enrichment \
  --results-input lof_pc_enrichment.tsv \
  --selection-output lof_pc_selection.json --plot-output lof_pc_enrichment.svg

rare-variant-enrichment export-zscore-matrix \
  --phenotype-bed phenotypes.bed.gz --principal-components pcs.tsv \
  --additional-covariates genetic_pcs.tsv --selection-input lof_pc_selection.json \
  --matrix-output selected_pc_z_scores.tsv.gz \
  --gene-qc-output selected_pc_z_scores.gene_qc.tsv.gz \
  --summary-output selected_pc_z_scores.summary.json
```

Legacy Python CLI commands remain available for compatibility, but are not part of the WDL or public analysis contract.

## Outputs

The workflow emits thirteen files, plus a haplo matrix when `ome_name` is `expression`:

- `results_tsv`: one row per PC count × negative threshold × carrier definition, merged across analysis shards.
- `summary_json`: selected grid/settings, global FDR scope, residualization description, provenance, and the screening limitation.
- `gene_pc_qc_tsv_gz`: compressed per-normalized-gene/per-PC QC with usable samples, rank, residual mean/SD, status, and exclusion reason.
- `analysis_qc_json`: BED/PC/covariate overlap counts, covariate names and input sample count when supplied, pre-join carrier-pair counts, LoF input QC, and per-PC eligibility, actual carrier-observation, and structured exclusion counters.
- `pc_selection_json`: median-logOR plateau summaries, excluded/included z thresholds, and the minimum common PC count selected by the 95% plateau rule.
- `selected_pc_haplo_calls_tsv_gz`: optional file, present only when `ome_name` is exactly `expression`; otherwise null. A compressed gene-by-sample haplo indicator matrix (`1`, `0`, or `NA`) at the selected PC count. Requires BED values on the log2-CPM scale; see the haplo rule below.
- `selected_pc_z_scores_tsv_gz`: compressed gene-by-sample matrix of signed residual Z scores at the selected PC count. Includes all genes in the phenotype BED.
- `selected_pc_z_scores_gene_qc_tsv_gz`: per-gene sample count, model rank, residual mean/SD, status, and exclusion reason for the matrix.
- `selected_pc_z_scores_summary_json`: selected PC count, matrix dimensions, covariate names, exclusion counts, and calculation rules.
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

### Selected-PC Z-score matrix

`ExportSelectedPcZScores` runs after PC selection. It reads `selection.selected_pc_count` from `pc_selection_json` and fits each gene with an intercept, all supplied additional covariates, and the first selected number of PCs. It divides residuals by their population standard deviation (`ddof=0`). The task uses the same calculation and exclusion rules as enrichment. It exports all genes in the BED, including noncoding genes and genes without LoF carriers. The enrichment and PC-selection steps still use protein-coding genes.

The gzip-compressed TSV starts with `gene_id`, followed by sample-ID columns. Genes follow their first occurrence in the BED. Samples follow BED order within the intersection of the BED, PC table, and optional covariate table. Multiple BED rows for one gene are collapsed to the minimum finite value per sample before adjustment. Missing observations have `NA` values. Genes that fail the model checks remain in the matrix with an all-`NA` row; the gene QC file gives the reason. No outlier threshold is applied to matrix values.

The WDL passes each input as `File` or `File?` until command rendering. The export command checks that inputs are readable and reports unresolved cloud paths as localization errors. Local tests cover cloud-to-local path substitution, optional covariates, and shell quoting. These tests do not validate a Terra run. The complete workflow with this export task has not been tested on Terra.

Use an image built from this revision before running the updated WDL. The existing GitHub Actions test job runs the complete fixture workflow, including this export task, in its test image. Omit `--additional-covariates` from the CLI commands when no additional covariates were used.

## Interpretation

Each test pools repeated samples and repeated genes, so Fisher p-values and globally adjusted FDR values are screening statistics, not confirmatory person-level inference. Use the raw cells, PC-specific QC, and appropriate dependence-aware models, permutations, or gene-level meta-analysis for confirmation. The workflow exports one Z-score matrix at the selected PC count.

### Haplo matrix at the selected PC count

Set `RareVariantEnrichment.ome_name` to `expression` to enable haplo export from a log2-CPM BED. This optional input defaults to an empty string, which disables haplo export. Matching is case-sensitive: only the exact label `expression` enables it. In the multi-omics workflow, the manifest `ome_name` supplies this label automatically.

For expression jobs, `selected_pc_haplo_calls_tsv_gz` identifies gene–sample pairs with a large expression drop. It has the same genes, sample IDs, and order as the Z-score matrix. The workflow includes all BED genes and uses the same sample intersection and duplicate-feature collapse as the Z-score export.

The haplo calculation fits an intercept, all supplied additional covariates, and the selected phenotype PCs. It uses the same covariates and aligned samples as the Z-score model. Subtracting a gene's mean from its adjusted log2-CPM gives the residual from this fit. The haplo calculation uses that residual in log2-CPM units; it does not divide by the residual SD.

- `1`: adjusted log2-CPM is strictly less than the gene mean minus `haplo_logcpm_drop`.
- `0`: the fit is valid and the value does not meet that rule, including exact equality.
- `NA`: the observation is missing or the fit is invalid. The fit requires a full-rank design and at least one residual degree of freedom.

`haplo_logcpm_drop` defaults to `1.0` and must be finite and non-negative. At zero selected PCs, the model still adjusts for all supplied additional covariates. With no additional covariates, it compares the original log2-CPM to the gene mean. Constant genes have no drop and receive `0` when the fit is valid, even when their Z scores cannot be calculated. Gene means use the finite observations in the aligned export cohort. For expression jobs, the existing export summary contains a `haplo` section with the threshold, adjustment rule, additional-covariate count, exclusion counts, and cell counts. Genes with insufficient observations or a rank-deficient full design receive `NA`.

This interpretation requires log2-CPM input. Applying a one-unit cutoff to Z scores, rank-normalized values, or splicing ratios does not give the same criterion. The workflow does not convert raw counts to log2-CPM. For the standalone export CLI, add `--haplo-matrix-output haplo.tsv.gz`; use `--haplo-logcpm-drop` to change the default drop. PC selection remains unchanged. Intersections that contain expression also require its haplo call to be `1`.

## Multiple matrices and multi-omics outliers

Use `workflows/multiomics_outliers.wdl` to run `RareVariantEnrichment` separately for each ome. Set `ome_manifest` to the TSV manifest file. See `examples/omics_manifest.tsv` and `examples/multiomics_outliers.inputs.json`. Supply the LoF carrier table for each ome in the manifest and one shared gene annotation. Each ome selects its own PC count and exports all BED genes as a Z-score matrix.

| Manifest column | Contents |
|---|---|
| `ome_name` | Unique ome name, such as `expression`, `protein`, or `splicing`. |
| `phenotype_bed` | Path to the phenotype BED for this ome. |
| `principal_components_tsv` | Path to the PC table for this ome. |
| `lof_carrier_table` | Required path to the LoF carrier table for this ome. |
| `additional_covariates_tsv` | Optional path to all fixed covariates for this ome. Use a blank cell or `.` to omit it. The column can also be omitted. |

Each row must supply `lof_carrier_table`. Repeat the same path when several omes share a table. This column replaces the top-level `MultiOmicsOutliers.lof_carrier_table` input; move that path into each row when updating an existing run. Each `matrix_results` record includes the carrier file used for that ome.

There is no separate ome-covariate field. Each row can use a different additional-covariate file, or several rows can use the same file. All supplied fixed covariates stay in the model while the workflow varies the PC count. Without that file, the model uses an intercept and the selected PCs.

Use `gs://` object paths for Terra. Absolute local paths are also accepted for tests; relative paths are rejected because localization changes the manifest directory. The manifest task checks the names, fields, and thresholds. It returns the paths as metadata. The workflow then declares each path as a `File` before it passes the value to an analysis task. The optional covariate path is declared only when supplied.

The final task creates a gene-by-sample indicator matrix for every dataset combination of size two or greater, at every `intersection_z_thresholds` value. The default thresholds are `-2, -3, -4, -5, -6`. These thresholds are separate from the enrichment and PC-selection thresholds. Each combination uses the same threshold in all participating datasets and tests the lower tail (`Z <= threshold`). If the exact label `expression` participates, its haplo call must also be `1`.

For expression, protein, and splicing, the combinations are expression–protein, expression–splicing, protein–splicing, and expression–protein–splicing. A protein–splicing cell is `1` when both values meet the threshold, regardless of expression. The three-way cell is `1` only when all three Z scores meet the threshold and expression also meets the haplo criterion. Expression–protein and expression–splicing pairs require that same expression haplo call. Pairwise and three-way outputs can therefore overlap.

| Cell | Meaning |
|---|---|
| `1` | All participating Z scores are finite and at or below the threshold, and expression has haplo = `1` when it participates. |
| `0` | All required values are present, but at least one Z score exceeds the threshold or the participating expression haplo call is `0`. |
| `NA` | A participating Z score or required expression haplo call is missing, even if another condition fails. |

Each combination uses its own shared genes and sample IDs. Row and column order follows the first dataset in that combination. Alignment uses IDs, not row or column position. A dataset outside the combination does not affect its sample set or calls. The expression haplo file must have the same gene and sample ID sets as the expression Z-score matrix; order can differ. A missing file or mismatched IDs fails validation instead of silently dropping pairs. For direct CLI use with an `expression` dataset, supply `--expression-haplo-file-list` containing exactly one localized haplo path. Omit it or supply an empty list when there is no expression dataset. Empty overlap produces an empty matrix and zero observation counts in the output index. Missing data are not classified as non-outliers.

The workflow outputs are:

- `matrix_results`: one `OmicsResult` per input dataset, in input order. It contains the dataset name, original input files, selected-PC Z-score matrix, optional expression haplo matrix, PC selection, enrichment results, plots, and QC files.
- `dataset_ids`, `z_score_matrices`, and `haplo_matrices`: corresponding arrays in manifest order for direct downstream use. `haplo_matrices` is an `Array[File?]`: only the `expression` entry has a file, and every other entry is null. If there is no expression row, all entries are null. `haplo_logcpm_drop` applies only to expression. Protein, splicing, and other labels produce no haplo file or haplo summary. The intersection task uses all participating Z scores and requires the haplo criterion only for combinations containing expression.
- `intersection_matrices`: gzip-compressed TSV indicator matrices. The first column is `gene_id`.
- `intersection_manifest_tsv`: maps each output basename to its datasets and threshold; gives gene, sample, outlier, non-outlier, and missing-cell counts. The `expression_haplo_required` column identifies combinations that use the haplo criterion. Use this index to identify files rather than relying on array order.
- `intersection_summary_json`: dataset dimensions, thresholds, calculation rules, and the same per-intersection counts.

For `N` datasets and `T` thresholds, the output count is `T * (2^N - N - 1)`. Three datasets at five thresholds produce 20 matrices. Output size grows quickly as datasets are added. The intersection task stores input gene vectors in a temporary SQLite database and processes one combination at a time. Its disk allocation has a 1,000 GB default floor, with a size-based increase; increase `intersection_disk_gb` for large input sets or many outputs. The other intersection resource inputs are `intersection_cpu` and `intersection_memory_gb`.

Dataset names must be unique, start with a letter, and contain only letters, digits, or underscores, up to 64 characters. At least two datasets are required. The workflow validates names and intersection thresholds before it starts the per-dataset analyses. A failed per-dataset analysis stops intersection generation; it does not silently omit that dataset.

After the manifest is read, all analysis input paths remain typed WDL `File` values until command rendering. The intersection task creates its newline-delimited list of local matrix paths at that point. The list contains no unresolved cloud URIs. The CLI rejects unreadable paths, duplicate IDs, invalid values, and inconsistent row widths. The single-matrix workflow remains available on its own.

The existing GitHub Actions Python test job includes workflow smoke tests for expression plus protein, and protein plus splicing without an expression row. That local fixture test enables miniwdl manifest file access (`MINIWDL__FILE_IO__ALLOW_ANY_INPUT=true`) because the local BED and PC paths originate in a trusted manifest. This test-only setting does not change Terra localization. Local checks cover ID alignment, missing values, threshold equality, pairwise and three-way calls, file localization, shell quoting, and WDL syntax. The complete multi-matrix workflow has not been tested on Terra. Build an updated image before running it. No cloud jobs are submitted by these tests.
