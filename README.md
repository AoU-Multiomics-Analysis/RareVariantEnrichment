# Rare variant enrichment

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml)
[![Python tests](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This repository has separate workflows for variant extraction and carrier
enrichment. Extraction creates a reusable variant audit. Generic enrichment
uses that audit to build named carrier definitions for LoF, missense, splice,
or combined variant classes. It then tests molecular low-expression outliers
with pooled Fisher exact 2×2 tables.

The original LoF-only enrichment workflow remains available for compatibility.
It does not accept variant call files, annotation tables, genomic index files,
or chromosome selections.

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

## Generic carrier enrichment

`workflows/carrier_enrichment.wdl` consumes the extraction audit and extraction
QC as a linked pair. It also consumes a carrier-definition JSON file, a
molecular phenotype BED, a PC matrix, and a protein-coding gene GTF. Additional
covariates are optional.

Run extraction first. Then replace the output paths in the enrichment example:

```bash
miniwdl run workflows/extract_variant_carriers.wdl \
  -i examples/extract_variant_carriers.inputs.json

miniwdl run workflows/carrier_enrichment.wdl \
  -i examples/carrier_enrichment.inputs.json
```

### Materialization

Materialization converts the row-level extraction audit into a sparse table of
sample, gene, and named carrier-definition records. It does not extract or
annotate variants again. The table contains these columns:

```text
sample_id  gene_id  gene_symbol  carrier_definition  n_variants  variant_ids
```

The materialization manifest records the ordered definitions, input hashes,
container image, row counts, and definition counts. The workflow checks that
the extraction QC belongs to the supplied audit before it builds the table.

Use `examples/carrier_definitions.json` as a starting point. Each definition
has a unique `name` and one or more `variant_classes`. The supported base
classes are `lof_hc`, `lof_hc_or_lc`, `missense`, `splice_core`, and
`splice_region`.

- Multiple `variant_classes` use OR logic.
- `minimum_revel` is optional. When it is present, the class rule and the REVEL
  rule must both match.
- A missing REVEL value does not pass a REVEL threshold.
- One audit row can enter more than one definition.
- A definition with no carriers stays in the manifest and downstream results.

The legacy LoF table is not expanded to represent other variant classes. The
generic workflow uses its own materialized carrier table. The LoF table and
`workflows/rare_variant_enrichment.wdl` remain compatibility interfaces.

## Enrichment inputs

The generic WDL has six required file inputs and one optional file input
(`additional_covariates_tsv`). The legacy LoF WDL has four required file inputs
and the same optional covariate input.

### Generic extraction and definition inputs

`variant_carrier_audit_tsv_gz` and `variant_carriers_qc_json` must come from the
same extraction run. `carrier_definitions_json` must use schema version 1 and
must contain an ordered, non-empty `definitions` array. Unknown fields, unknown
classes, duplicate names, invalid REVEL thresholds, and duplicate JSON keys are
rejected.

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

### Legacy LoF carrier table

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

Build or select a container image that contains this package and NumPy 2 or
later. The supplied image uses micromamba and exact package pins from
conda-forge and bioconda. It includes R, tidyverse, `ggplot2`, and `ggrepel` for
the PC-sweep QC plot.

The enrichment defaults support cohort-scale input localization. Preparation
uses 2 CPUs, 32 GB RAM, and 500 GB disk. Analysis and merge tasks use 8 CPUs,
128 GB RAM, and a 1000 GB disk baseline. PC chunk preparation uses 1 CPU and
4 GB RAM. Each task increases its disk request when the localized input size
requires more space. `pc_preemptible` controls the preemptible retry count for
the PC scatter. `max_retries` applies to all tasks.

Run the generic workflow with the audit and QC outputs from extraction:

```bash
miniwdl run workflows/carrier_enrichment.wdl \
  -i examples/carrier_enrichment.inputs.json
```

For reproducibility, replace `docker_image` with an immutable image digest in
production.

The Python CLI exposes the generic materialization and enrichment operations:

```bash
rare-variant-enrichment build-carrier-definitions \
  --audit variant_carrier_audit.tsv.gz \
  --extraction-qc variant_carriers.qc.json \
  --definitions examples/carrier_definitions.json \
  --container-image ghcr.io/example/image@sha256:... \
  --output carrier_definitions.tsv.gz \
  --qc-output carrier_definitions.qc.json

rare-variant-enrichment prepare-protein-coding-genes \
  --gtf genes.gtf.gz --genes-output protein_coding_genes.tsv --qc-output protein_coding_genes.qc.json

rare-variant-enrichment carrier-pc-enrichment \
  --phenotype-bed phenotypes.bed.gz \
  --carrier-table carrier_definitions.tsv.gz \
  --carrier-manifest carrier_definitions.qc.json \
  --principal-components pcs.tsv --protein-coding-genes protein_coding_genes.tsv \
  --additional-covariates genetic_pcs.tsv \
  --negative-z-thresholds=-2,-3,-4,-5,-6 --pc-counts '' \
  --results-output carrier_pc_enrichment.tsv \
  --summary-output carrier_pc_enrichment.summary.json \
  --gene-pc-qc-output carrier_pc_enrichment.gene_pc_qc.tsv.gz \
  --analysis-qc-output carrier_pc_enrichment.analysis_qc.json
```

Run the legacy LoF interface only when an existing process requires its table
or output contract:

```bash
miniwdl run workflows/rare_variant_enrichment.wdl \
  -i examples/rare_variant_enrichment.inputs.json
```

Legacy Python CLI commands remain available for compatibility.

## Generic enrichment outputs

The generic workflow emits 12 files:

- `carrier_definitions_tsv_gz`: materialized sample–gene carrier records for
  each configured definition.
- `carrier_definitions_qc_json`: definition order, counts, hashes, and extraction
  provenance.
- `results_tsv`: one row per PC count × negative threshold × carrier
  definition, merged across PC shards.
- `summary_json`: grid settings, global FDR scope, residualization details,
  provenance, and the screening limitation.
- `gene_pc_qc_tsv_gz`: compressed gene-by-PC QC with usable sample counts,
  model rank, residual summaries, status, and exclusion reason.
- `analysis_qc_json`: sample overlap, carrier counts, and structured PC-specific
  exclusion counters.
- `pc_selection_json`: definition-specific plateau summaries and the minimum
  common PC count that meets the configured plateau fraction.
- `enrichment_plot_svg`: threshold and median log-odds curves for the selected
  definitions.
- `pc_sweep_qc_summary_tsv`: PC-sweep values, maxima, and percentages of the
  maximum median odds ratio.
- `pc_sweep_qc_plot_png`: a clean percentage-of-maximum QC plot.
- `protein_coding_genes_tsv`: the sorted coding-gene list.
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

`analysis_qc_json` reconciles each PC-specific `carrier_observations` count with
every result row: `carrier_observations = n11 + n10`. Pre-join carrier counts
are separate because a carrier that is absent from the phenotype and PC sample
set does not enter a Fisher table.

PC-shard files are workflow intermediates. During merge, the workflow recomputes
the global Benjamini–Hochberg FDR across all final PC, threshold, and definition
rows. It does not retain shard-local adjustments.

By default, PC selection uses all materialized definitions and the thresholds
`z = -3, -4, -5, -6`. It finds the first PC count that reaches 95% of each
definition's maximum median log odds ratio. It then reports the smallest common
PC count that satisfies all estimable definitions. A zero-carrier or otherwise
unestimable definition stays in the output with a structured exclusion reason.

The QC summary and plot use the same definition and threshold selections. Each
estimable definition is normalized to its own maximum median odds ratio.

### Legacy output contract

`workflows/rare_variant_enrichment.wdl` emits the same ten downstream analysis
files without the two materialization files. It retains the `HC`, `HC_or_LC`,
and `any_lof` carrier definitions and the legacy file contract.

## Interpretation

Each test pools repeated samples and repeated genes, so Fisher p-values and globally adjusted FDR values are screening statistics, not confirmatory person-level inference. Use the raw cells, PC-specific QC, and appropriate dependence-aware models, permutations, or gene-level meta-analysis for confirmation. The workflow deliberately does not export full residual matrices.
