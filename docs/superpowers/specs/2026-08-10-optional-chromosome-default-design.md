# Optional chromosome default design

## Goal

Allow callers to omit the workflow chromosome input. An omitted input runs all
human autosomes using the repository's `chr` naming convention. A supplied
array continues to restrict the workflow to exactly those chromosomes.

## Public WDL contract

Keep the existing input type and name:

```wdl
Array[String] chromosomes
```

Add this WDL default:

```wdl
Array[String] chromosomes = [
  "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
  "chr19", "chr20", "chr21", "chr22"
]
```

A defaulted WDL input is optional for callers while remaining a concrete
`Array[String]` inside the workflow. This avoids nullable-array resolution
and leaves every existing task and scatter interface unchanged.

The default intentionally excludes `chrX`, `chrY`, mitochondrial contigs,
and alternate contigs. Callers may include any of those explicitly when their
BED and VCF use matching contig names.

## Data flow

The resolved `chromosomes` array continues to feed:

1. VCF index/contig validation.
2. Phenotype feature selection.
3. The one-task-per-chromosome scatter.
4. Final selected-chromosome provenance and denominator scoping.

An explicit input such as `["chr7"]` therefore produces a chr7-only run.
An explicit input such as `["chr1", "chr2"]` produces a two-chromosome run.
Omitting the input produces 22 scatter shards.

## Documentation and examples

- Remove `RareVariantEnrichment.chromosomes` from the primary example input
  JSON so the example demonstrates the autosomal default.
- Document the 22-autosome default and show single- and multi-chromosome
  overrides.
- Retain the requirement that chromosome labels exactly match BED and VCF
  contigs.

## Validation

- Update the parsed WDL contract test to assert the full default array.
- Update the miniwdl input-template assertion so `chromosomes` is absent from
  required inputs.
- Preserve existing Docker-backed explicit-chr1 runtime tests to prove
  overrides still restrict execution.
- Add a parsed workflow assertion that the scatter and downstream calls use the
  resolved `chromosomes` value.
- Run the complete test suite, `miniwdl check`, JSON validation, and
  `git diff --check`.

## Compatibility

Existing input JSON files that supply `chromosomes` remain valid and behave
unchanged. The only interface change is that callers may now omit the field.
