# Generic Variant Extraction and Enrichment Design

## Purpose

Provide separate extraction and enrichment workflows that support configurable
variant carrier definitions. The extraction workflow creates one canonical
variant-carrier audit. The enrichment workflow converts that audit into named
carrier definitions and tests all selected definitions in one PC sweep.

The change must support LoF, missense, REVEL-thresholded missense, splice-core,
splice-region, and union definitions without repeating the VCF-to-transcript
annotation join.

## Scope

### In scope

- Keep `workflows/extract_variant_carriers.wdl` as an independent workflow.
- Use `variant_carrier_audit.tsv.gz` as the biological source table.
- Add configurable carrier-definition construction to the new enrichment
  workflow.
- Add `workflows/carrier_enrichment.wdl` as the generic enrichment entry point.
- Test any ordered set of configured definitions in one enrichment run.
- Keep the current LoF WDL, CLI commands, inputs, outputs, and result values as
  compatibility interfaces.
- Publish carrier-definition, enrichment, plot, and QC outputs.
- Add unit, integration, WDL-contract, and GitHub Actions runtime tests.
- Update the runtime image to a version-pinned micromamba environment because
  the generic plot path changes its R dependencies.

### Out of scope

- A combined extraction-and-enrichment WDL submission.
- New VCF quality, genotype-quality, AC, AF, MAF, or annotation-frequency
  filters.
- New transcript annotations or a new VCF-to-VAT join.
- A general rule-expression language.
- SpliceAI thresholds, regulatory annotations, ancestry-stratified tests, or
  dependence-aware confirmatory statistics.
- Removal of the existing LoF workflow or legacy distance-based CLI commands.

The input VCF remains a prefiltered source. Extraction records frequency values
for audit purposes but does not apply new frequency filters.

## User-facing workflows

### ExtractVariantCarriers

`workflows/extract_variant_carriers.wdl` remains the only workflow that reads
the filtered VCF and transcript annotations. It does not accept phenotype or PC
inputs, and it does not run enrichment.

Its primary outputs remain:

```text
variant_carrier_audit_tsv_gz
variant_carriers_tsv_gz
variant_carriers_qc_json
```

The audit contains exact allele, sample dosage, normalized gene, consequence,
LoFTEE, REVEL, frequency, and initial-class facts. The six-column
`variant_carriers_tsv_gz` remains a base-class compatibility output. It is not
the canonical input for generic enrichment.

Add an `audit_artifact` object to gathered extraction QC. It must contain the
audit header, data-row count, byte size, and SHA-256 digest. Calculate size and
digest after the deterministic gzip stream closes. The generic enrichment
workflow uses this object to bind the localized audit to its extraction QC.

Use this exact object shape:

```json
{
  "logical_name": "variant_carrier_audit.tsv.gz",
  "header": [
    "sample_id", "gene_id", "gene_symbol", "chrom", "pos", "ref", "alt",
    "variant_id", "variant_ac", "variant_af", "sample_alt_allele_count",
    "most_severe_consequence", "all_consequences", "unknown_consequences",
    "loftee", "revel", "gvs_max_af", "variant_classes"
  ],
  "row_count": 0,
  "size_bytes": 0,
  "sha256": "64 lowercase hexadecimal characters"
}
```

Do not change either extraction table header in this compatibility release.

### CarrierEnrichment

Add `workflows/carrier_enrichment.wdl`. It accepts the extraction audit,
extraction QC, carrier-definition configuration, phenotype and covariate data,
and the current PC-sweep settings.

Required file inputs:

```text
File variant_carrier_audit_tsv_gz
File variant_carriers_qc_json
File carrier_definitions_json
File phenotype_bed
File principal_components_tsv
File gene_annotation_gtf
```

Optional file input:

```text
File? additional_covariates_tsv
```

Analysis and runtime inputs retain the current defaults where applicable:

```text
Array[Float] negative_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
Array[Float] selection_z_thresholds = [-3.0, -4.0, -5.0, -6.0]
Float plateau_fraction = 0.95
Array[Int] pc_counts = []
Int pc_counts_per_job = 10
Array[String] pc_selection_carrier_definitions = []
Int pc_preemptible = 2
String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
Int prepare_cpu = 2
Int prepare_memory_gb = 32
Int prepare_disk_gb = 500
Int analysis_cpu = 8
Int analysis_memory_gb = 128
Int analysis_disk_gb = 1000
Int max_retries = 1
```

An empty `pc_selection_carrier_definitions` array selects all configured
definitions in configuration order. A non-empty array must contain unique
configured names. This input controls PC-selection summaries and plots. It does
not remove enrichment rows for other configured definitions.

The workflow publishes:

```text
carrier_definitions_tsv_gz
carrier_definitions_qc_json
results_tsv
summary_json
gene_pc_qc_tsv_gz
analysis_qc_json
pc_selection_json
enrichment_plot_svg
pc_sweep_qc_summary_tsv
pc_sweep_qc_plot_png
protein_coding_genes_tsv
protein_coding_genes_qc_json
```

The workflow contains these stages:

1. `BuildCarrierDefinitions` validates provenance and constructs carrier
   definitions.
2. `PrepareProteinCodingGenes` and `PreparePcChunks` prepare the current gene
   set and PC grid.
3. `CalculateCarrierPcEnrichment` scatters only over PC-count chunks.
4. `MergeCarrierPcEnrichment` merges shards and recalculates global FDR.
5. `AnalyzeCarrierPcEnrichment` creates dynamic PC-selection and plot outputs.

## Carrier-definition construction

Carrier-definition construction is the first task in
`carrier_enrichment.wdl`. The workflow publishes its table and QC so another
analysis can reuse them. WDL call caching can reuse the task when the audit,
extraction QC, definition configuration, and container are unchanged.

The Python CLI command is:

```text
build-carrier-definitions
```

The command accepts the audit, extraction QC, definition configuration,
container reference, output table, and output QC paths.

### Configuration schema

Use this JSON form:

```json
{
  "schema_version": 1,
  "definitions": [
    {"name": "lof_hc", "variant_classes": ["lof_hc"]},
    {"name": "lof_hc_or_lc", "variant_classes": ["lof_hc_or_lc"]},
    {"name": "missense", "variant_classes": ["missense"]},
    {
      "name": "missense_revel_ge_0_75",
      "variant_classes": ["missense"],
      "minimum_revel": 0.75
    },
    {"name": "splice_core", "variant_classes": ["splice_core"]},
    {"name": "splice_region", "variant_classes": ["splice_region"]},
    {
      "name": "splice_any",
      "variant_classes": ["splice_core", "splice_region"]
    }
  ]
}
```

The top-level object contains exactly `schema_version` and `definitions`.
Schema version 1 supports only these base classes:

```text
lof_hc
lof_hc_or_lc
missense
splice_core
splice_region
```

Within one definition, `variant_classes` uses OR. `minimum_revel`, when
present, uses AND with the class rule. A row with missing REVEL does not meet a
REVEL threshold. The threshold must be finite and from 0 through 1, inclusive.

Definition names must match `[A-Za-z][A-Za-z0-9_.-]*`. Definition and base-class
order is significant. Reject duplicate names, duplicate base classes, empty
definitions, unknown keys, unknown classes, Boolean numeric values, unsupported
schema versions, and duplicate JSON keys.

### Materialized table

The gzip-compressed TSV has this exact header:

```text
sample_id
gene_id
gene_symbol
carrier_definition
n_variants
variant_ids
```

For each deduplicated audit row, evaluate every definition. One audit row can
meet more than one definition. Group matches by normalized sample ID,
normalized gene ID, and definition. Count distinct exact `variant_id` values
and store them in deterministic order.

Use `(sample_id, normalized gene_id, chrom, pos, ref, alt)` as the duplicate
identity. Trim surrounding ASCII whitespace from sample IDs and reject empty
values. Require an Ensembl gene ID of `ENSG` followed by digits and an optional
terminal numeric version. Remove the version. Treat empty or `.` gene symbols
as missing. Fail when duplicate identities disagree in any audit fact or when
one normalized gene has two different non-empty symbols. Backfill the one
available non-empty symbol into every output row for that gene.

Sort output rows by sample ID, gene ID, and definition configuration order. Use
deterministic gzip metadata. Produce a valid header-only file when no audit row
meets any definition.

Use bounded SQLite storage for audit deduplication, definition matches,
aggregation, and sorting. Do not retain the complete cohort audit or all
sample-gene carrier sets in memory during construction.

### Definition QC

Write schema `aou.carrier-definitions-manifest.v1`. Include:

- definition order and normalized definitions;
- input audit header, row count, byte size, and SHA-256 digest;
- the embedded extraction QC digest and relevant extraction provenance;
- definition configuration byte size and SHA-256 digest;
- input, deduplicated, and duplicate audit-row counts;
- present and missing REVEL counts;
- matched audit rows, distinct variants, distinct sample-gene pairs, and output
  rows for every definition, including zero-count definitions;
- output header, row count, byte size, and SHA-256 digest;
- package version and complete container reference.

Fail before construction when the audit does not match the `audit_artifact`
object in extraction QC.

The output-artifact object in the manifest uses the same exact keys as
`audit_artifact`: `logical_name`, `header`, `row_count`, `size_bytes`, and
`sha256`. The logical name is `carrier_definitions.tsv.gz`. Store the complete
validated definition configuration and the ordered definition names in the
manifest. JSON maps that hold per-definition counts must include every
definition in configuration order.

## Generic enrichment engine

### Carrier model

Add one internal model:

```python
@dataclass(frozen=True)
class CarrierSets:
    definitions: tuple[str, ...]
    pairs_by_definition: dict[str, set[tuple[str, str]]]
    qc: dict[str, object]
```

The generic carrier reader validates the materialized table against
`aou.carrier-definitions-manifest.v1`. It preserves declared definitions that
have zero carrier rows.

Refactor the calculation, merge, and selection code to use
`CarrierSets.definitions` instead of the fixed `CARRIER_DEFINITIONS` constant.

Add these generic CLI commands:

```text
carrier-pc-enrichment
merge-carrier-pc-enrichment
analyze-carrier-pc-enrichment
```

The commands accept the materialized table and its manifest. They do not read
the full carrier audit. The legacy command names call the same internal engine
through the legacy LoF adapter.

### Residualization and tests

For each gene and PC count:

1. Calculate the residual z-score vector once.
2. Build one carrier mask for each configured definition.
3. Apply all definition masks to the same outlier mask.
4. Emit one contingency row per definition and negative z-score threshold.

Do not scatter by carrier definition. Keep the current PC-count chunk scatter.
The number of definitions must not increase the number of residualization or
projection calls.

The generic result TSV keeps the current result columns. The
`carrier_definition` field contains the configured name. Calculate
Benjamini-Hochberg FDR across all final PC-count, threshold, and definition
rows after shard merge.

## PC selection and plots

Remove fixed `HC` and `any_lof` assumptions from the shared Python and R plot
logic. Pass the selected ordered definition list explicitly.

Use the selected definitions for median log-odds summaries, plateau selection,
plot facets, colors, labels, and output order. Keep the existing LoF colors and
labels when the legacy wrapper selects `HC` and `any_lof`.

If a selected definition has no complete finite odds-ratio curve across the
selection thresholds, keep its contingency results but exclude it from the
plateau calculation. Record `zero_carriers`, `zero_observations`, or
`incomplete_finite_curve` in selection QC. If no selected definition is
estimable, write valid selection and plot outputs with
`selected_pc_count = null`.

R plots must use tidyverse syntax and a clean, minimal aesthetic. Do not add a
plot title or subtitle.

## Runtime image

Update `envs/Dockerfile` to use a pinned micromamba base image. Install pinned
runtime versions from `conda-forge` and `bioconda`, including Python, NumPy,
htslib, R, tidyverse, ggplot2, and ggrepel. Install the local Python package
after the conda environment is complete. Do not leave an unbounded direct
runtime dependency in the image definition.

Build and exercise this image in GitHub Actions. Do not build it locally only
for a smoke test. Production examples continue to recommend an immutable image
digest.

## QC and provenance

Generic analysis QC records the ordered definition list and these counts for
each definition:

- unique input sample-gene pairs;
- pairs with a sample in the phenotype, PC, and optional-covariate
  intersection;
- pairs with a protein-coding gene;
- pairs in the complete analysis universe;
- eligible carrier observations for every PC count.

Include explicit zero values. All PC shards must report the same definition
order, carrier-table metadata, and materializer-manifest digest. Merge fails
when any of those values differ.

The summary records the definition configuration, PC grid mode, covariate
model, outlier rule, global FDR scope, package version, container reference,
and the screening-statistic limitation.

## Legacy LoF compatibility

Keep the seven-column LoF reader. Convert its rows in memory into
`CarrierSets` with this exact order:

```text
any_lof
HC
HC_or_LC
```

Keep these CLI commands as adapters to the generic engine:

```text
lof-pc-enrichment
merge-lof-pc-enrichment
analyze-lof-pc-enrichment
```

Keep `workflows/rare_variant_enrichment.wdl` and its current required
`lof_carrier_table` input. Preserve its task behavior, ten public output types,
result schema, JSON keys and types, definition order, plots, and values. It
continues to use `HC` and `any_lof` for PC selection.

Do not encode missense or splice carriers as rows in the legacy LoF table.

## Error handling

Fail with a clear message for:

- an invalid audit, extraction QC, definition configuration, or materializer
  manifest;
- mismatched audit size, digest, header, or row count;
- conflicting duplicate audit rows or gene symbols;
- an invalid sample, gene, allele, REVEL value, definition, or variant ID;
- an undeclared carrier definition in the materialized table;
- a materialized table that does not match its manifest;
- inconsistent shard definition order or provenance;
- a requested PC-selection definition that is not configured.

Do not fail only because a valid definition has zero variants, carrier pairs,
or estimable odds ratios. Report the condition in QC.

## Logging

Every new or changed WDL command must log start, input counts, selected
definitions, PC chunk, periodic progress, output counts, and completion.
Extraction keeps its chromosome and chunk progress logs. Logs must not contain
sample IDs or carrier rows.

## Documentation and registration

- Add an example carrier-definition JSON file.
- Add an example `CarrierEnrichment` input JSON file.
- Document the two-workflow handoff in the README.
- State that extraction expects a prefiltered VCF.
- State that definition construction applies only class unions and optional
  REVEL thresholds.
- Register `workflows/carrier_enrichment.wdl` in `.dockstore.yml`.
- Keep both existing workflow registrations.

## Verification

Use test-first development for every behavior change.

### Unit and integration tests

- Configuration validation, including duplicate JSON keys and invalid values.
- Values below, equal to, and above a REVEL threshold, plus missing REVEL.
- LoF, missense, splice-core, splice-region, and union definitions.
- One audit row matching multiple definitions.
- Duplicate audit rows and conflicting duplicate facts.
- Deterministic row order and deterministic gzip bytes.
- Header-only materialization with explicit zero QC.
- Manifest size, digest, header, and row-count binding.
- Four or more definitions in one enrichment invocation.
- Complete-data projections independent of definition count.
- Missing-data residualization independent of definition count.
- Arbitrary ordered definitions through calculation, merge, selection, Python
  plots, and R plots.
- Zero-carrier and unestimable-definition selection behavior.
- Exact legacy LoF Fisher cells, QC, JSON schemas, and public WDL contract.

### WDL validation and runtime tests

- Validate every WDL with `miniwdl check`.
- Test the exact `CarrierEnrichment` inputs, task wiring, scatter dimension,
  runtime settings, logging, and outputs.
- Add one hand-checkable end-to-end fixture:

```text
carrier audit + extraction QC + definition JSON + phenotype BED + PCs + GTF
  -> carrier definitions
  -> generic enrichment
```

The fixture contains HC, LC, missense below and at a REVEL threshold,
splice-core, splice-region, overlapping definitions, and hand-checked Fisher
cells.

Run the Docker-backed smoke test in GitHub Actions with the CI-built image. Do
not build a Docker image locally only for the smoke test.

## Success criteria

The change is complete when:

1. Extraction remains an independent workflow and publishes a bound canonical
   audit.
2. `carrier_enrichment.wdl` constructs configured definitions from that audit
   and tests them in one run.
3. LoF, REVEL-thresholded missense, and splice definitions produce
   hand-checked enrichment cells.
4. The number of definitions does not increase residualization calls.
5. Generic QC and manifests reconcile all carrier and enrichment counts.
6. The existing LoF WDL remains value-compatible.
7. Unit, integration, WDL-contract, WDL-runtime, and GitHub Actions tests pass.
