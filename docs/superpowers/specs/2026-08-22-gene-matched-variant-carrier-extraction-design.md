# Gene-Matched Variant Carrier Extraction Design

## Purpose

Create a standalone WDL workflow that extracts gene-matched variant carriers from an already filtered VCF and its transcript-level All of Us Variant Annotation Table (VAT) annotations.

The workflow must produce two reusable data products:

1. A per-sample, per-gene, per-variant carrier audit.
2. An aggregated sample-gene carrier table stratified by variant class.

The workflow must not call the enrichment workflow. Enrichment remains a separate analysis step.

This work also updates MT-to-VCF so that its `TranscriptAnnotations` output retains REVEL scores.

## Scope

### In scope

- Add `revel` to the MT-to-VCF transcript-annotation export.
- Add an independent carrier-extraction WDL to RareVariantEnrichment.
- Join filtered VCF alleles to transcript annotations by exact allele identity.
- Match variants to normalized VAT gene IDs.
- Collapse transcript annotations per variant-gene pair.
- Assign the approved initial variant classes.
- Emit a detailed carrier audit, an aggregated carrier table, and QC outputs.
- Add unit, integration, WDL-contract, and GitHub Actions smoke tests.
- Add explicit logging to all new or changed WDL commands.

### Out of scope

- Calling or changing the enrichment workflow.
- Applying new genotype-quality, variant-quality, allele-frequency, or MAF filters.
- Selecting variants by distance from a transcription start site or phenotype interval.
- Adding a REVEL-based damaging-missense class in this change.
- Adding a SpliceAI-based class in this change.

The input VCF is the authoritative filtered variant set. Frequency values are retained for audit and provenance only.

## Repositories

This design changes two repositories.

### MTtoVCF

The transcript export must retain `revel` for every available transcript annotation. The change must update the exported schema, schema validation, tests, and documentation. It must not change MT filtering or the variant-level VCF and annotation outputs.

### RareVariantEnrichment

The repository receives a standalone extraction workflow and the Python commands that it invokes. Existing parsing and VAT consequence utilities should be reused where their semantics match this design. The active LoF PC-enrichment workflow remains unchanged.

## Workflow interface

The standalone workflow accepts:

- A filtered, bgzipped VCF.
- A tabix VCF index.
- The bgzipped MT-to-VCF `TranscriptAnnotations` table.
- An optional transcript-annotation tabix index.
- A chromosome list, with autosomes as the default.
- Container and task runtime inputs.

The workflow emits:

- A gzip-compressed per-variant carrier audit.
- A gzip-compressed aggregated carrier table.
- Per-chromosome QC JSON files.
- One gathered QC and provenance JSON file.
- The generated or validated transcript-annotation index.

If the transcript-annotation index is not supplied, the workflow creates and validates one. The transcript annotations must be coordinate-sorted and use columns 1 and 2 for chromosome and position.

The workflow does not accept a phenotype BED, GTF, feature table, gene list, distance threshold, quality threshold, or frequency threshold.

## Architecture and data flow

The workflow has four stages.

### 1. Prepare and validate inputs

- Read the VCF header and sample IDs.
- Validate the requested contigs against the VCF index.
- Read the transcript-annotation header.
- Validate all required fields.
- Validate a supplied transcript index or create one when it is absent.
- Record index provenance.

The required transcript fields are:

```text
chrom
pos
ref
alt
gene_id
gene_symbol
transcript
consequence
LoF
gvs_max_af
revel
```

Other existing transcript fields remain available but are not required by the carrier classification logic.

### 2. Extract chromosome carriers

Scatter one task per selected chromosome. Each task reads the VCF and transcript annotations in bounded coordinate chunks.

For each VCF ALT allele:

1. Construct the exact allele key from chromosome, position, REF, and ALT.
2. Find transcript rows with the same exact allele key.
3. Group matching transcript rows by normalized gene ID.
4. Collapse each variant-gene group according to the annotation rules below.
5. Parse each VCF genotype and calculate the sample ALT-allele count.
6. Emit one audit record for each sample with an ALT-allele count greater than zero and each matched gene.
7. Emit zero or more class-carrier records for that audit record.

The task must not use distance or phenotype position when it assigns a variant to a gene.

### 3. Gather the per-variant audit

Gather and deduplicate chromosome audit records. Sort the output deterministically by sample ID, gene ID, chromosome, position, REF, and ALT.

### 4. Aggregate carriers

Reduce audit records by normalized sample ID, normalized gene ID, and variant class. Count distinct variants and retain a deterministic list of distinct variant IDs.

## Annotation collapse rules

### Gene normalization

Remove only a terminal numeric version suffix from an Ensembl gene ID. Use the normalized gene ID as the primary gene key. Preserve the gene symbol as metadata.

During gather, build one gene-symbol mapping for each normalized gene ID. If records for one normalized gene ID contain different non-empty gene symbols, fail with an explicit error. Missing symbols are allowed. Use the one available non-empty symbol when other records for the gene have no symbol.

### Consequence selection

For each exact variant-gene pair:

- Parse all consequence terms from all matching transcript rows.
- Deduplicate the complete recognized term set.
- Select one most-severe consequence using the documented Ensembl consequence severity order already used by RareVariantEnrichment.
- Preserve the full deduplicated term set in the audit.
- Count unknown terms in QC.

Only the selected most-severe consequence controls consequence-class assignment.

### LoFTEE selection

Collapse LoFTEE independently from the consequence selection:

- Select `HC` if any matching transcript row has `LoF=HC`.
- Otherwise select `LC` if any matching transcript row has `LoF=LC`.
- Otherwise leave LoFTEE empty.

A variant can enter one consequence class and one or more LoF-derived classes.

### REVEL selection

Parse all non-missing REVEL values for the variant-gene pair. Store the maximum finite value in the audit.

Missing REVEL values are allowed. A missing result remains empty. Nonnumeric, non-finite, or out-of-range values are invalid. Valid REVEL values must be from 0 through 1, inclusive.

### Frequency selection

Use the ALT-specific `AC` and `AF` values from the input VCF. Preserve missing values when those INFO fields are absent.

For `gvs_max_af`, parse all non-missing values for the variant-gene pair and store the maximum finite value. Missing values are allowed. Nonnumeric, non-finite, or out-of-range values are invalid. Valid frequency values must be from 0 through 1, inclusive. Frequency selection does not filter a variant.

## Initial variant classes

The extraction workflow assigns these initial classes:

| Variant class | Rule |
| --- | --- |
| `lof_hc` | Collapsed LoFTEE is `HC`. |
| `lof_hc_or_lc` | Collapsed LoFTEE is `HC` or `LC`. |
| `missense` | Most-severe consequence is `missense_variant`. |
| `splice_core` | Most-severe consequence is `splice_acceptor_variant` or `splice_donor_variant`. |
| `splice_region` | Most-severe consequence is `splice_donor_5th_base_variant`, `splice_region_variant`, `splice_donor_region_variant`, or `splice_polypyrimidine_tract_variant`. |

An HC variant enters both `lof_hc` and `lof_hc_or_lc`. A variant can also enter its selected consequence class. Each sample-gene-variant combination counts only once within a class.

The workflow does not assign a REVEL-based class in this version. The audit retains REVEL so a later change can define a configurable rule such as `missense AND REVEL >= threshold` without repeating the VCF-VAT join.

## Per-variant carrier audit

The required audit is a gzip-compressed TSV with one row per sample-gene-variant carrier:

```text
sample_id
gene_id
gene_symbol
variant_id
chrom
pos
ref
alt
sample_alt_allele_count
cohort_ac
cohort_af
gvs_max_af
most_severe_consequence
all_gene_consequences
loftee
revel
variant_classes
```

Requirements:

- `variant_id` uses a deterministic `chrom:pos:ref:alt` representation.
- `sample_alt_allele_count` records the count of the matched ALT allele in the sample genotype.
- `cohort_ac` and `cohort_af` come from the input VCF when available.
- `gvs_max_af` comes from the gene-matched VAT rows.
- `all_gene_consequences` is a deterministic comma-separated list.
- `variant_classes` is a deterministic comma-separated list and can be empty.
- The audit retains all exact VCF alleles with at least one gene-matched VAT annotation and at least one carrier, even when no initial variant class applies.

## Aggregated carrier table

The aggregated output is a gzip-compressed TSV:

```text
sample_id
gene_id
gene_symbol
variant_class
n_variants
variant_ids
```

Requirements:

- Emit one row per sample-gene-class carrier.
- Count distinct exact alleles.
- Store a deterministic comma-separated list of distinct variant IDs.
- An absent sample-gene-class row means noncarrier.
- Produce a valid header-only file when no class carriers exist.

This table is designed for a later generalized enrichment interface. This change does not connect it to the current LoF enrichment command.

## Error handling

The workflow must fail for:

- Missing or duplicate required columns.
- An invalid VCF or transcript-annotation index.
- Requested contigs that are absent from either indexed input.
- Malformed VCF records or genotype fields.
- Invalid gene IDs.
- Conflicting non-empty gene symbols for one normalized gene ID in the gathered audit.
- Nonnumeric, non-finite, or out-of-range REVEL values.
- Nonnumeric, non-finite, or out-of-range frequency values when present.
- Allele-key inconsistencies within a queried coordinate region.

The workflow must not fail only because:

- A VCF allele has no VAT match.
- A VAT annotation has no VCF allele.
- A valid VAT row has missing REVEL, LoFTEE, gene symbol, or frequency values.
- A chromosome has no carriers or no class-assigned carriers.

These cases must be counted in QC.

## QC and provenance

Emit per-chromosome QC and one gathered JSON summary. Include at least:

- Input VCF records and ALT alleles examined.
- Called and missing genotypes.
- Carrier audit rows emitted.
- Unique samples, genes, alleles, and variant-gene pairs.
- VCF alleles joined and not joined to VAT.
- VAT rows read and duplicate rows removed.
- Recognized and unknown consequence terms.
- Selected consequence counts.
- LoFTEE HC, LC, and missing counts.
- REVEL present and missing counts, plus the observed minimum and maximum.
- Carrier counts for each initial variant class.
- Header schemas, selected chromosomes, input paths, and index provenance.
- A statement that the input VCF was treated as prefiltered and no quality or frequency filters were applied.

## Logging

Every new or changed WDL command must log clear start, progress, and completion messages. Chromosome tasks must log the chromosome, query chunks, processed allele counts, joined VAT counts, audit row counts, and class-carrier counts. Logs must not print sample-level carrier records.

## Verification

Development must use tests before implementation changes.

### MTtoVCF tests

- Verify that `revel` is part of the transcript schema.
- Verify that transcript export selects and writes REVEL.
- Verify missing REVEL handling.
- Update WDL interface and documentation tests as required.

### RareVariantEnrichment tests

- Unit tests for gene normalization, consequence collapse, LoFTEE collapse, maximum REVEL selection, and genotype ALT dosage.
- Unit tests for each initial class and allowed overlap.
- Tests that repeated transcript and genotype records do not double-count variants.
- Tests that carrier extraction does not use TSS distance or phenotype positions.
- Tests for missing VAT matches and header-only outputs.
- An end-to-end fixture with a bgzipped VCF, tabix index, transcript-annotation BGZ, and annotation index.
- WDL contract and syntax validation tests.
- QC reconciliation tests between the audit and aggregated carrier table.

### Pipeline smoke test

Add a GitHub Actions smoke test that runs the extraction workflow on small fixtures. Do not build a Docker image locally only for the smoke test. Existing published or CI-built images may be used according to the repository CI pattern.

## Compatibility

- Do not change the active RareVariantEnrichment LoF PC workflow or its inputs and outputs.
- Do not change MT-to-VCF filtering behavior.
- Preserve all existing MT-to-VCF transcript columns and add `revel` as a new column.
- Keep deterministic output ordering so reruns with the same inputs produce identical tables.

## Success criteria

The change is complete when:

1. MT-to-VCF emits REVEL in `TranscriptAnnotations` and its tests pass.
2. The standalone WDL accepts a filtered indexed VCF and transcript annotations.
3. Alleles are assigned to genes only through exact allele and normalized VAT gene matching.
4. Transcript rows collapse to the approved consequence, LoFTEE, and REVEL values.
5. Both required tables and QC outputs are correct and deterministic.
6. No enrichment step runs as part of extraction.
7. Unit, integration, WDL-validation, and GitHub Actions smoke tests pass.
