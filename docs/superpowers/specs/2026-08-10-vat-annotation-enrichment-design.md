# VAT Annotation Enrichment Design

**Date:** 2026-08-10
**Status:** Proposed; approved in conversation and awaiting written-spec review

## Goal

Extend the pooled rare-variant outlier enrichment workflow so it can test
gene-matched functional variant classes from the All of Us Variant Annotation
Table (VAT). The extension must preserve the existing z-score, allele-count,
and TSS-distance analyses while keeping memory bounded for cohort-scale VCF and
VAT inputs.

The first annotation release will test exact coding consequence terms and
LoFTEE confidence classes. It will also use `gvs_max_af` as a second rarity
filter so variants that are common in an All of Us ancestry group do not enter
the enrichment analysis even if they are present in the rare-variant VCF.

## Source data and row granularity

The workflow accepts a bgzipped, coordinate-sorted VAT TSV plus an optional
tabix index. The official VAT contains one row per variant-transcript
combination, but the workflow must also accept pre-collapsed VAT derivatives
with one row per variant or variant-gene pair.

Header discovery will recognize both official VAT coordinate names and the
aliases used by the SuSiE annotation code:

| Meaning | Official name | Accepted alias |
|---|---|---|
| Chromosome | `contig` | `chrom` |
| Position | `position` | `pos` |
| Reference allele | `ref_allele` | `ref` |
| Alternate allele | `alt_allele` | `alt` |

`gene_id`, `consequence`, and `gvs_max_af` are required. `LoF` is optional
because it is not guaranteed by the published VAT schema. If `LoF` is absent,
the workflow runs consequence and baseline enrichment but emits no LoFTEE
strata. Missing or ambiguous required aliases fail schema validation.

The allele join key is the exact `(chromosome, position, reference,
alternate)` tuple. A multiallelic VCF record is split logically into one key
per ALT before joining. The workflow does not trim or left-normalize alleles;
VCF and VAT inputs must use the same GRCh38 representation. Join-rate QC is
used to detect incompatible representations.

## Gene matching and transcript collapse

The phenotype BED contains versioned Ensembl gene IDs. Both BED `feature_id`
and VAT `gene_id` values are normalized by removing the terminal Ensembl
version suffix before comparison. Functional annotations are assigned only
when the normalized VAT gene matches the phenotype gene. A coding consequence
for an overlapping neighboring gene must not classify the tested gene.

For each `(variant allele, normalized gene_id)`:

1. Union all consequence terms across every matching VAT row and every element
   of each row's `consequence` array.
2. Select exactly one most-severe term using the documented Ensembl consequence
   ordering, from more severe to less severe.
3. Assign `LoFTEE=HC` if any matching row has `LoF=HC`; otherwise assign
   `LoFTEE=LC` if any matching row has `LoF=LC`; otherwise assign no LoFTEE
   class. Values are matched case-insensitively, and HC takes precedence over
   LC.

The selected consequence and LoFTEE class are separate annotation families.
A variant can therefore contribute to one exact consequence stratum and one
LoFTEE stratum. Within the consequence family, a variant-gene pair has only one
class. Within the LoFTEE family, HC and LC are mutually exclusive.

Most-severe collapse is per variant-gene pair, not per sample-gene observation.
If one sample carries two qualifying variants for a gene with different
selected consequences, that sample-gene observation is a carrier in both
corresponding consequence analyses.

Unknown consequence terms are recorded in QC. A variant with no recognized or
configured gene-matched consequence may still enter the baseline analysis if
it passes the frequency filter.

## Consequence classes

The full versioned Ensembl severity hierarchy determines which term is selected.
A configurable `consequence_classes` WDL input controls which selected terms
receive enrichment rows. The coding-focused defaults are:

- `splice_acceptor_variant`
- `splice_donor_variant`
- `stop_gained`
- `frameshift_variant`
- `stop_lost`
- `start_lost`
- `inframe_insertion`
- `inframe_deletion`
- `missense_variant`
- `protein_altering_variant`
- `splice_region_variant`
- `synonymous_variant`
- `coding_sequence_variant`

The workflow also always evaluates a baseline class named
`all_rare_variants`. Broad classes such as pLoF or damaging missense can be
added later without changing allele extraction or carrier-state semantics.

## Frequency normalization and filtering

`gvs_max_af` is a variant-level, transcript-independent field. Every non-missing
value must parse as a number in `[0, 1]`. The workflow converts it to the minor
allele-frequency spectrum per row:

```text
gvs_max_maf = min(gvs_max_af, 1 - gvs_max_af)
```

The default inclusive filter is `gvs_max_maf <= 0.01`, exposed as
`maximum_gvs_maf`. Missing, non-numeric, out-of-range, and above-threshold
values are excluded and counted separately in QC. The observed raw maximum and
the number of values above 0.5 are reported; the workflow does not require the
observed maximum to equal 0.5 because a rare or region-restricted subset may
have a lower maximum.

Repeated transcript rows for one allele are expected to have the same
`gvs_max_af`. Identical values are deduplicated. If parsed values disagree
beyond floating-point tolerance, the allele is excluded as an inconsistent
frequency annotation and counted in QC.

An allele absent from the VAT or lacking valid `gvs_max_af` is excluded from
all enrichment families, including the baseline, because its rarity cannot be
confirmed. Gene-matched consequence is not required for the baseline once the
variant-level frequency passes.

## Public WDL inputs

Add the following workflow inputs:

```wdl
File variant_annotation_table
File? variant_annotation_table_tbi
Float maximum_gvs_maf = 0.01
Int annotation_chunk_size_bp = 10000000
Array[String] consequence_classes = [ ...coding defaults above... ]
```

Existing phenotype, VCF, chromosome, z-threshold, AC-class, distance, tail,
runtime, retry, and carrier-audit inputs remain in effect. The VAT input is
required for the annotation-enabled workflow revision. The optional index
follows the same supplied-versus-generated pattern as the VCF index.

## Workflow architecture

### 1. PrepareVatIndex

Read the VAT header, resolve accepted aliases, and emit a schema manifest for
downstream tasks. If `variant_annotation_table_tbi` is supplied, localize it
beside the VAT and validate it. Otherwise, create a generic tabix index using
the detected chromosome and position columns while skipping the TSV header.

Index generation requires the bgzipped VAT to be coordinate-sorted. An
unsorted file fails with a clear error rather than triggering an implicit
whole-file sort. The task records whether the index was supplied or generated
and emits the generated or validated index.

### 2. PreparePhenotypes and maximum windows

Retain the existing phenotype preparation and one-base TSS interpretation.
Normalize a second internal copy of each feature ID for VAT gene matching while
preserving the original versioned feature ID in carrier and result interfaces.

Build the union of all TSS windows at the largest requested distance. Smaller
distance thresholds are not queried separately.

### 3. Chunked chromosome classification

Continue scattering one WDL task per selected chromosome. Within each task,
split the union of maximum-distance regions into non-overlapping genomic chunks
no larger than `annotation_chunk_size_bp`. Chunk boundaries are inclusive and
non-overlapping so a genomic position is extracted exactly once per input.

For each chunk:

1. Query the VAT with tabix and stream its rows into a temporary disk-backed
   SQLite annotation store.
2. Maintain variant-level frequency state keyed by allele and gene-level
   collapsed consequence/LoFTEE state keyed by allele plus normalized gene ID.
3. Query the VCF with tabix and stream records without materializing all chunk
   genotypes or variants in Python memory.
4. Split each record by ALT, join its variant-level VAT frequency, apply the
   MAF filter, determine existing AC classes, and find every feature TSS within
   the maximum distance using the chromosome's sorted feature index.
5. For each carrier and nearby feature, emit baseline state and any qualifying
   gene-matched consequence or LoFTEE state into the persistent chromosome
   carrier store.
6. Drop the temporary annotation store before moving to the next chunk.

The persistent chromosome carrier store remains disk-backed and performs a
minimum-distance upsert. Querying one gene at a time is intentionally avoided:
large TSS windows overlap, so gene-at-a-time tabix calls would repeatedly read
and parse the same VAT and VCF records. Non-overlapping chunks preserve bounded
memory while reading each covered genomic position once.

### 4. Annotation-aware carrier reduction

Extend the carrier key to:

```text
(sample_id, feature_id, ac_class, annotation_family, annotation_class)
```

Store only `minimum_distance_bp` for each key. This is a compact representation
of carrier status at every distance threshold: a key is a carrier at threshold
`w` exactly when `minimum_distance_bp <= w`. Z-score thresholds remain phenotype
properties and are evaluated later rather than being expanded into the carrier
table.

The gathered carrier audit schema is:

```text
sample_id
feature_id
ac_class
annotation_family
annotation_class
minimum_distance_bp
```

As before, the sample-level audit is published only with explicit opt-in.

### 5. CalculateEnrichment

Extend the existing pooled calculation across:

```text
z threshold x AC class x distance x annotation family x annotation class
```

For each combination, a gene-sample observation is a carrier if at least one
qualifying carrier key has `minimum_distance_bp` at or below the threshold.
Construct the same pooled 2x2 table and retain the existing observation,
outlier-tail, rate-ratio, odds-ratio, Fisher exact, and screening-statistic
limitations.

Benjamini-Hochberg FDR is calculated globally across every emitted combination,
including baseline, consequence, and LoFTEE rows.

## Outputs

Extend `rare_variant_enrichment.tsv` with:

- `annotation_family`: `baseline`, `consequence`, or `loftee`.
- `annotation_class`: `all_rare_variants`, an exact selected consequence, `HC`,
  or `LC`.

All existing threshold, AC definition, distance, cell-count, carrier-rate,
effect-size, p-value, and FDR columns remain. Stable sorting includes annotation
family and class.

Extend the JSON provenance/QC output with VAT path metadata, resolved schema,
MAF threshold, chunk size, severity-order version, configured consequence
classes, LoFTEE availability, and index provenance. Emit a gathered annotation
QC table and the generated or validated VAT index in addition to the existing
outputs.

## QC and failure behavior

Per-chromosome and gathered annotation QC includes:

- chunks and tabix queries processed;
- VAT rows queried and duplicate rows observed;
- unique VAT alleles and allele-gene pairs;
- VCF ALT alleles queried and joined to VAT;
- gene-matched and unmatched allele-feature pairs;
- observed raw `gvs_max_af` maximum and values converted from above 0.5;
- missing, non-numeric, out-of-range, inconsistent, and above-threshold
  frequency exclusions;
- consequence arrays parsed, recognized terms, unknown terms, and configured
  classes selected;
- HC, LC, missing, and unrecognized LoFTEE values;
- baseline, consequence, and LoFTEE carrier keys emitted; and
- VAT-index source.

Fail early for missing/ambiguous required columns, an invalid MAF threshold or
chunk size, unusable compression/indexing, requested contigs absent from the
VAT index, and malformed required coordinates or alleles. If eligible queried
VCF alleles exist across the run but none match VAT allele keys, fail the
gathered analysis with a representation/build mismatch error. Individual
missing or malformed annotations are excluded and reported rather than leaking
variant-level VAT data into ordinary QC outputs.

## Testing strategy

Unit tests will cover:

- official and SuSiE-style header aliases;
- versioned Ensembl gene normalization and gene-specific matching;
- supported consequence-array encodings, null arrays, and duplicate terms;
- full severity ordering and most-severe collapse across transcript rows;
- HC-over-LC LoFTEE collapse and absent `LoF` behavior;
- AF-to-MAF conversion, inclusive 1% filtering, invalid values, and
  transcript-level frequency inconsistency;
- multiallelic VCF allele-to-VAT joins; and
- annotation-aware minimum-distance upserts.

Integration and WDL tests will cover:

- one-row-per-transcript and pre-collapsed VAT fixtures;
- generated and supplied generic VAT tabix indexes;
- non-overlapping chunk boundaries and exact maximum-distance boundaries;
- each covered genomic position being queried once rather than once per gene or
  threshold;
- a hand-checkable baseline, consequence, and LoFTEE enrichment result;
- a sample carrying multiple consequence classes for one gene;
- bounded Python heap with disk-backed annotation and carrier stores;
- WDL input/output contracts, command-argument safety, and static validation;
  and
- Docker-backed miniature end-to-end runs for supplied and generated indexes.

## Non-goals

This revision does not implement annotation-weighted Watershed inference,
ancestry-stratified enrichment, covariate adjustment, regulatory annotation
classes, pathogenicity-score thresholds, gene-specific association tests, or
implicit VAT sorting/normalization. The pooled Fisher tests remain exploratory
screens with repeated genes and participants, as documented by the base
workflow.

## References

- [All of Us Variant Annotation Table schema](https://support.researchallofus.org/hc/en-us/articles/4615256690836-Variant-Annotation-Table)
- [Ensembl calculated consequence severity order](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html)
- [Existing SuSiE VAT loading workflow](https://github.com/AoU-Multiomics-Analysis/susieR/blob/main/workflows/AnnotateSusie.wdl)
