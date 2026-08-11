# Rare variant enrichment

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/docker-image.yml)
[![Python tests](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/RareVariantEnrichment/actions/workflows/python-tests.yml)
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

`INFO/AC` is interpreted independently for every ALT and is authoritative when it is a non-negative integer. A missing value (`AC=.` or a `.` entry such as `AC=.,2`) falls back only that ALT to the called genotype allele indices across every VCF sample, including samples absent from the BED. Negative AC values are rejected. On fully called records, INFO and genotype AC are compared and disagreements are counted while INFO remains authoritative.

The partial-call policy is explicit: a known ALT in a call such as `1/.` is a carrier and contributes one called allele to genotype-derived AC; the unknown allele contributes nothing. `./.` and `.` are fully missing and never create carriers. If all genotype alleles are missing for an ALT without INFO AC, its AC is unavailable and that ALT is skipped. Chromosome QC distinguishes complete, partial, and fully missing calls and reports INFO, genotype-fallback, unavailable, compared, mismatch, and unchecked per-ALT counts. Missing phenotype z-scores remove the corresponding feature–sample observation even when the sample carries a qualifying ALT.

### Thresholds and runtime

The defaults and recommended first-pass grid are:

- `z_thresholds`: `[2.0, 3.0, 4.0, 5.0]`.
- `exact_allele_counts`: `[1, 2, 3, 4, 5]`, producing `AC=N` classes.
- `cumulative_allele_count_maxima`: `[1, 2, 3, 5, 10]`, producing independent `AC<=N` classes.
- `distance_thresholds_bp`: `[1000, 10000, 100000, 1000000]`, applied symmetrically and inclusively as `abs(VCF_POS - TSS) <= distance`.
- `outlier_tail`: `absolute`; `positive` and `negative` are also supported.

Exact and cumulative AC arrays are independent, so either may be empty, but at least one family must be configured. Threshold values must be valid and unique. The workflow queries each chromosome once over merged windows at the largest distance, then applies all smaller thresholds from exact distances.

The default runtime image is `ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main`. The GHCR package must be public, or the execution backend must have credentials that can pull it. For production, override `docker_image` with an immutable digest such as `ghcr.io/aou-multiomics-analysis/rarevariantenrichment@sha256:...`; tags can move and are not sufficient provenance. CPU, memory, and disk inputs can be overridden separately for preparation, scatter, and gather/calculation tasks. `max_retries` defaults to `1` and is applied as WDL `maxRetries` to every task.

`publish_carrier_audit` defaults to `false`. Set it to `true` only when the sample-level carrier-distance audit is needed and its access controls are appropriate.

## Run

Install [miniwdl](https://miniwdl.readthedocs.io/) and a supported container runtime. From a directory containing `phenotypes.scaled.residualized.bed.gz` and `rare_variants.vcf.gz`, run:

```bash
miniwdl run workflows/rare_variant_enrichment.wdl \
  -i examples/rare_variant_enrichment.inputs.json
```

Adjust the two file paths and chromosome naming in the example JSON for the input cohort. To use an existing index, add `"RareVariantEnrichment.rare_variant_vcf_tbi": "rare_variants.vcf.gz.tbi"`.

The selected chromosome array controls the exact feature set in both preparation and calculation. Rows on other BED chromosomes are validated during preparation but do not enter enrichment denominators.

## Outputs

The workflow publishes:

- `enrichment_tsv`: one row per z threshold × AC class × distance combination.
- `enrichment_json`: analysis inputs; phenotype overlap/QC; per-chromosome and summed variant QC; selected chromosomes; index source; retry count; container image; Python, SQLite, package, and workflow versions; observation/carrier counts; coordinate convention; and the statistical-limitation statement.
- `phenotype_qc_json`: selected/input feature counts, selected chromosomes, non-missing/missing and threshold-specific outlier counts, and BED/VCF overlap counts. It intentionally contains counts rather than excluded sample IDs.
- `chromosome_qc_tsv`: per-chromosome feature, merged-region, extracted-record, total/classified ALT, AC-source/unavailable/mismatch, complete/partial/fully-missing genotype, maximum-distance boundary, variant–feature-pair, tabix-query, and emitted-key counts.
- `carrier_minimum_distances_tsv`: optional deduplicated `(sample_id, feature_id, ac_class, minimum_distance_bp)` audit records, emitted only when `publish_carrier_audit=true`. Multiple qualifying variants reduce to the minimum distance, so one feature–sample observation is counted at most once per AC class. This output contains sample IDs and should be treated as controlled data.
- `generated_or_validated_vcf_tbi`: the generated or validated tabix index.
- `vcf_index_provenance`: machine-readable `generated` or `supplied` index provenance.
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

BED-only and VCF-only samples are excluded from analysis. Published QC reports overlap counts without unnecessarily exposing identifiers. A carrier must have a qualifying known ALT allele and a non-missing phenotype for the same feature–sample observation. Unknown genotype alleles are not evidence of reference homozygosity; review partial, fully missing, unavailable-AC, and mismatch counters before interpreting enrichment.

## Scale and determinism

Feature TSS positions are sorted once per chromosome and queried by binary search. Classification, cross-chromosome gather, and final carrier joins use temporary on-disk SQLite tables with minimum-distance upserts rather than cohort-sized Python dictionaries. The phenotype matrix is streamed one selected feature row at a time, and Fisher's exact test uses a mode-relative constant-memory recurrence. Carrier and result TSVs use explicit stable sort orders.

The CI scale regression creates only 40,000 gather keys and 30,000 calculation keys, then checks bounded Python heap rather than constructing a cohort-sized fixture. The current test ceilings are 4 MiB for gather, 8 MiB for calculation, and 1 MiB for a Fisher table with 250,001 feasible cells. SQLite page cache and temporary files consume bounded native memory/disk separately, so production disk settings still need to cover the actual carrier table.

Samples recur across genes and genes recur across samples, so pooled gene–sample observations are not statistically independent. Fisher p-values and BH FDR values are screening statistics, not confirmatory inference. Use the emitted cell counts and carrier-distance table for dependence-aware mixed models, permutation tests, or gene-level meta-analysis before drawing causal conclusions.

This first pass intentionally does not implement functional consequence or external frequency annotations, ancestry adjustment, gene-specific tests, or Watershed inference. The per-ALT AC classification and deduplicated carrier-distance interface are extension points for future annotation filters; significant screening results can later feed a separate Watershed-based prioritization stage without changing the current observation definition.
