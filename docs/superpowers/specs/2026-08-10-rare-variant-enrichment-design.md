# Rare Variant Enrichment Design

**Date:** 2026-08-10  
**Status:** Proposed

## Goal

Build a first-pass WDL workflow for testing whether rare variant carriers are enriched among RNA-seq/proteomics outlier observations across the All of Us multi-omics cohort. The workflow will evaluate multiple absolute z-score thresholds and allele-count classes, while pooling evidence across genes to avoid sparse per-gene tests.

## Scope

This version will provide a screening analysis for rare variants already present in a VCF. It will not yet implement functional annotations, broader allele-frequency annotations, watershed causal modeling, covariate adjustment, ancestry-stratified tests, or a final gene-level causal statistic.

The analysis unit is a gene/feature–sample observation, not a sample with any variant anywhere. A variant contributes to an observation only when the variant is mapped to that observation's feature/gene and the sample carries the alternate allele.

## Statistical definition

For each z-score threshold `t` and allele-count class `c`:

1. Read every non-missing `(sample_id, feature_id, z_score)` observation.
2. Classify it as an outlier when `abs(z_score) >= t` by default. The workflow will also expose a tail setting for positive-only or negative-only outliers.
3. Mark the observation as a carrier when the sample carries at least one variant mapped to the same feature and belonging to class `c`.
4. Construct the pooled 2x2 table:

   | | Carrier | Non-carrier |
   |---|---:|---:|
   | Outlier | `a` | `b` |
   | Non-outlier | `c` | `d` |

5. Report the carrier-rate ratio `(a/(a+b))/(c/(c+d))`, the odds ratio `(a*d)/(b*c)`, a two-sided Fisher exact p-value, and Benjamini–Hochberg FDR across the tested threshold × allele-count combinations.

The output will include all underlying cell counts so downstream users can apply alternative filters or statistics. Zero-cell odds ratios will be reported using a clearly labelled 0.5 continuity correction, while the uncorrected cell counts remain available.

Because samples and genes recur in the pooled table, the Fisher test is a first-pass screening statistic rather than a fully independent-observation inference. The report will state this limitation and retain enough counts for a later mixed-model, permutation, or gene-level meta-analysis implementation.

## Allele-count classes

The workflow will support two independently evaluated class families:

- Exact classes: `AC=1`, `AC=2`, `AC=3`, etc.
- Cumulative classes: `AC<=1`, `AC<=2`, `AC<=3`, etc.

The WDL will accept the requested exact counts and cumulative maxima as arrays, so the caller can choose singleton-only, singleton-plus-doubleton, or broader rare-variant analyses without changing code. Each variant is assigned an allele count from `INFO/AC` when available; otherwise the helper calculates AC from the genotype fields in the VCF.

## Inputs

The initial workflow will use these inputs:

- A bgzipped VCF containing the rare variants and all relevant sample genotypes.
- An optional VCF tabix index. If omitted, a preprocessing task will generate `VCF.tbi` using `tabix -p vcf`; if supplied, it will be copied to the expected local name and validated.
- A long-format z-score TSV with header columns `sample_id`, `feature_id`, and `z_score`. Additional columns may be retained as metadata but will not affect the first-pass statistic.
- A chromosome array whose names exactly match the VCF contig names, such as `1`, `2`, ..., `22`, `X` or `chr1`, `chr2`, ..., `chrX`.
- Either a configurable VCF INFO field containing feature/gene identifiers, or an optional variant-to-feature mapping TSV. The mapping-file path will be preferred when supplied, allowing annotations to be generated independently of the VCF.
- Arrays of z-score thresholds, exact AC values, and cumulative AC maxima.
- An optional outlier-tail setting with `absolute` as the default and `positive`/`negative` as alternatives.

The VCF must be bgzip-compressed and coordinate-sorted for tabix. The preprocessing task will fail early with a clear validation error if these conditions are not met.

## Workflow architecture

The WDL will contain four logical stages:

### 1. PrepareVcfIndex

Normalize the optional index input:

- If an index is supplied, copy it beside the VCF and validate it with `tabix`.
- If no index is supplied, run `tabix -p vcf` and emit the generated index.
- Validate that the VCF is bgzip-readable and obtain the available contig names.

### 2. ExtractAndClassifyChromosome

Scatter over the requested chromosome array. Each task will:

- Use `tabix -h indexed.vcf.gz chromosome` to extract only that chromosome.
- Parse the extracted records and sample genotypes.
- Determine each variant's allele count and requested AC classes.
- Resolve each variant to one or more feature IDs using the mapping TSV or configured INFO field.
- Emit a deduplicated carrier-pair table with `(sample_id, feature_id, ac_class)` plus per-chromosome variant/QC counts.

The task will not calculate enrichment independently because a gene–sample pair can have qualifying variants on more than one chromosome. It will emit keys that can be safely deduplicated after gather.

### 3. GatherCarrierPairs

Combine the scattered carrier-pair tables and deduplicate on `(sample_id, feature_id, ac_class)`. This prevents a gene–sample observation with multiple qualifying variants or variants on multiple chromosomes from being counted multiple times.

### 4. CalculateEnrichment

Join the deduplicated carrier pairs with the z-score observations, evaluate all threshold × class combinations, calculate the 2x2 counts and statistics, and write the final report and QC summary.

## Outputs

The workflow will emit:

- `rare_variant_enrichment.tsv`: one row per z-score threshold × AC class, including threshold, tail, class definition, total observations, outlier/non-outlier counts, carrier/non-carrier cells, carrier-rate ratio, continuity-corrected odds ratio, Fisher p-value, and BH FDR.
- `rare_variant_enrichment.json`: run parameters, input summaries, and high-level QC metrics.
- `chromosome_qc.tsv`: per-chromosome records processed, variants retained, variants lacking feature mappings, and carrier pairs emitted.
- `deduplicated_carrier_pairs.tsv`: the gathered carrier keys, useful for auditing and downstream methods.

The final report will also include counts of z-score observations dropped for missing values, features absent from the variant mapping, samples absent from the VCF, malformed/missing genotypes, and variants skipped because their AC could not be determined.

## Error handling and reproducibility

The helper will fail for malformed required columns, duplicate z-score keys with conflicting values, invalid allele-count parameters, and an unusable VCF/index. It will skip and count non-fatal records such as missing genotypes or variants without a feature mapping.

All thresholds, AC class arrays, tail mode, mapping mode, VCF INFO field, and software/container versions will be recorded in the JSON output. The WDL will expose runtime settings for CPU, memory, disk, and task retry behavior rather than hard-coding cohort-specific values.

## Testing strategy

Tests will use tiny synthetic fixtures and will cover:

- Exact and cumulative AC class assignment.
- Positive, negative, and absolute outlier classification.
- Correct construction of the 2x2 table and Fisher p-value on a hand-checkable example.
- Deduplication of the same gene–sample carrier across multiple variants and chromosomes.
- INFO-field versus mapping-file feature resolution.
- Missing genotypes, unmapped variants, and missing z-scores in QC counts.
- WDL static validation and a local miniature end-to-end execution where the required WDL runtime is available.

## Future extension points

The interfaces intentionally leave room for:

- Broader allele-frequency and functional annotation classes.
- Gene-level burden definitions and feature-specific variant windows.
- Ancestry or technical covariate adjustment.
- Permutation or mixed-model inference that accounts for repeated samples and genes.
- Watershed-based variant prioritization beneath significant outlier signals.
