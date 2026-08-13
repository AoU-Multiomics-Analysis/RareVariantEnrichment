# Enrichment Progress Logging

## Goal

Make long LoF/PC enrichment runs observable in workflow logs without changing
analysis results or emitting high-volume per-sample or per-gene output.

## Design

`lof-pc-enrichment` will always emit concise progress messages to standard
error. Standard output remains unused so it cannot interfere with tabular or
JSON outputs.

The command will log:

1. Validated analysis configuration: z-score thresholds, PC grid, number of
   protein-coding genes, and carrier-pair counts.
2. Phenotype BED and PC sample overlap after reading the BED header.
3. Progress after every 500 protein-coding BED genes processed.
4. Completion of each PC count after the BED has been processed, including
   eligible-gene and observation totals.
5. Completion with result-row count and the paths of all four analysis output
   files.

## Constraints

- Logging is always enabled; there is no verbosity flag.
- Messages must not alter output file contents or the statistical calculation.
- Progress reporting must be bounded: no individual sample or gene IDs are
  logged during the loop.

## Verification

Tests will capture standard error from a small real enrichment fixture and
assert that configuration, sample-overlap, periodic progress, PC completion,
and output-completion messages are emitted. Existing result and workflow tests
will continue to verify unchanged analysis behavior.
