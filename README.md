# Rare variant enrichment

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This WDL 1.0 workflow asks whether molecular outliers are enriched among carriers of rare variants near a feature's transcription start site (TSS). It is a pooled first-pass screen: the observation unit is one non-missing gene/feature–sample z-score, not one person and not one gene. For each z-score threshold, allele-count (AC) class, and TSS-distance threshold, the workflow forms a pooled 2×2 outlier-by-carrier table.

## Inputs

### Molecular phenotype BED

`phenotype_bed` is a wide prepare_QTL BED, optionally gzip- or bgzip-compressed:

```text
#chr  start  end  gene_id  SAMPLE_1  SAMPLE_2  ...
chr1  999    1000 GENE1    2.7       NA        ...
```

The first four columns are chromosome, start, end, and a positional feature identifier; the fourth column's header may differ. Every feature interval must be exactly one base wide. Under the prepare_QTL convention it represents the zero-based half-open interval `[TSS - 1, TSS)`, so the workflow uses the BED `end` as the one-based TSS compared with one-based VCF `POS`.

Columns five onward are sample IDs and scaled-residual z-scores. Empty values, `.`, `NA`, and `NaN` are missing and are excluded from the pooled observation denominator rather than treated as non-outliers. Sample and feature IDs must be unique.

### Rare-variant VCF

`rare_variant_vcf` must be coordinate-sorted, bgzip-compressed, and include contig declarations, a `GT` FORMAT field, and sample genotypes. Chromosome names in `chromosomes` must exactly match both BED and VCF contigs. Each ALT allele of a multiallelic record is evaluated separately.

`rare_variant_vcf_tbi` is optional. When omitted, `PrepareVcfIndex` generates an adjacent tabix index with `tabix -p vcf`; when supplied, it is localized beside the VCF and validated. The example intentionally omits this field to demonstrate automatic index generation.

`INFO/AC` supplies per-ALT allele counts when present. Otherwise AC is derived from genotypes across every VCF sample, including samples absent from the BED. Only sample IDs shared by the BED and VCF can become carrier observations. Missing genotype alleles do not create carrier calls and are counted in chromosome QC; missing phenotype z-scores remove the corresponding gene–sample observation even if that sample is a carrier.

### Thresholds and runtime

The defaults and recommended first-pass grid are:

- `z_thresholds`: `[2.0, 3.0, 4.0, 5.0]`.
- `exact_allele_counts`: `[1, 2, 3, 4, 5]`, producing `AC=N` classes.
- `cumulative_allele_count_maxima`: `[1, 2, 3, 5, 10]`, producing independent `AC<=N` classes.
- `distance_thresholds_bp`: `[1000, 10000, 100000, 1000000]`, applied symmetrically and inclusively as `abs(VCF_POS - TSS) <= distance`.
- `outlier_tail`: `absolute`; `positive` and `negative` are also supported.

Exact and cumulative AC arrays are independent, so either may be empty, but at least one family must be configured. Threshold values must be valid and unique. The workflow queries each chromosome once over merged windows at the largest distance, then applies all smaller thresholds from exact distances.

The default runtime image is `ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main`. CPU, memory, and disk inputs can be overridden separately for preparation, scatter, and gather/calculation tasks.

## Run

Install [miniwdl](https://miniwdl.readthedocs.io/) and a supported container runtime. From a directory containing `phenotypes.scaled.residualized.bed.gz` and `rare_variants.vcf.gz`, run:

```bash
miniwdl run workflows/rare_variant_enrichment.wdl \
  -i examples/rare_variant_enrichment.inputs.json
```

Adjust the two file paths and chromosome naming in the example JSON for the input cohort. To use an existing index, add `"RareVariantEnrichment.rare_variant_vcf_tbi": "rare_variants.vcf.gz.tbi"`.

## Outputs

The workflow publishes:

- `enrichment_tsv`: one row per z threshold × AC class × distance combination.
- `enrichment_json`: run parameters, observation/carrier counts, missing-z count, coordinate convention, and the statistical-limitation statement.
- `chromosome_qc_tsv`: per-chromosome feature, merged-region, extracted-record, ALT, missing-genotype, variant–feature-pair, tabix-query, and emitted-key counts.
- `carrier_minimum_distances_tsv`: deduplicated `(sample_id, feature_id, ac_class, minimum_distance_bp)` audit records. Multiple qualifying variants reduce to the minimum distance, so one gene–sample observation is counted at most once per AC class.
- `generated_or_validated_vcf_tbi`: the generated or validated tabix index.
- `chromosome_query_regions`: one merged maximum-window BED per chromosome.

The enrichment TSV columns are:

| Column | Definition |
|---|---|
| `z_threshold` | Inclusive z-score cutoff used with `tail`. |
| `tail` | `absolute`, `positive`, or `negative` outlier rule. |
| `ac_class` | Human-readable exact (`AC=N`) or cumulative (`AC<=N`) class. |
| `ac_kind` | `exact` or `cumulative`. |
| `ac_value` | Numeric exact count or cumulative maximum. |
| `distance_bp` | Inclusive symmetric maximum absolute TSS distance in base pairs. |
| `total_observations` | Non-missing gene–sample observations among shared BED/VCF samples. |
| `outlier_observations` | Observations meeting the z-score/tail rule. |
| `nonoutlier_observations` | Non-missing observations not meeting the rule. |
| `n11` | Outlier carrier observations. |
| `n10` | Outlier non-carrier observations. |
| `n01` | Non-outlier carrier observations. |
| `n00` | Non-outlier non-carrier observations. |
| `outlier_carrier_rate` | `n11 / (n11 + n10)`, or `NA` if undefined. |
| `nonoutlier_carrier_rate` | `n01 / (n01 + n00)`, or `NA` if undefined. |
| `carrier_rate_ratio` | Outlier carrier rate divided by non-outlier carrier rate, or `NA` if undefined. |
| `odds_ratio` | Uncorrected `(n11 × n00) / (n10 × n01)`, or `NA` when its denominator is zero. |
| `odds_ratio_corrected_0_5` | Odds ratio after adding 0.5 to every 2×2 cell. |
| `fisher_p_value` | Two-sided Fisher exact p-value from the uncorrected cells. |
| `fisher_fdr_bh` | Benjamini–Hochberg adjustment across every emitted threshold combination. |

## Interpretation and limitations

BED-only and VCF-only samples are excluded from analysis. The enrichment JSON reports the shared-sample count, while the `PreparePhenotypes` task records the excluded sample lists in its `phenotype_qc.json` execution artifact. A carrier must have a qualifying ALT genotype and a non-missing phenotype for the same feature–sample observation. A missing genotype is treated as no observed carrier allele, not as evidence of reference homozygosity, and its frequency should be reviewed in chromosome QC.

Samples recur across genes and genes recur across samples, so pooled gene–sample observations are not statistically independent. Fisher p-values and BH FDR values are screening statistics, not confirmatory inference. Use the emitted cell counts and carrier-distance table for dependence-aware mixed models, permutation tests, or gene-level meta-analysis before drawing causal conclusions.

This first pass intentionally does not implement functional consequence or external frequency annotations, ancestry adjustment, gene-specific tests, or Watershed inference. The per-ALT AC classification and deduplicated carrier-distance interface are extension points for future annotation filters; significant screening results can later feed a separate Watershed-based prioritization stage without changing the current observation definition.
