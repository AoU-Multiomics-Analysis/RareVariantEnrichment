# Scattered LoF/PC Enrichment Design

## Goal

Parallelize the LoF/PC enrichment sweep across groups of selected PC-count settings while preserving the existing six public workflow outputs and global statistical summaries.

## Requirements

- Add a workflow input named `pc_counts_per_job`, defaulting to `10`.
- Interpret `pc_counts_per_job` as the number of selected PC-count settings assigned to each job, not the number of PCA columns.
- Preserve both explicit `pc_counts` and the existing adaptive PC-count grid when `pc_counts` is empty.
- Reject non-positive `pc_counts_per_job` values before scattering.
- Keep one analysis task per PC-count chunk; do not create one task per PC-count setting.
- Merge all shard outputs into the existing public result, summary, gene-PC QC, and analysis-QC files.
- Recompute Fisher FDR globally across all result rows after merging; shard-local FDR values must not be retained.
- Preserve result row definitions, output headers, QC fields, JSON fields, and existing residualization behavior.

## Workflow architecture

The workflow will add a `PreparePcChunks` task before the scatter. It will inspect only the principal-component header to validate the available PC count, construct the explicit or adaptive grid using the same rules as the Python implementation, and write the grid as JSON chunks of at most `pc_counts_per_job` values.

The workflow will scatter `CalculateLofPcEnrichment` over the chunk array. Each shard receives one `Array[Int]` of selected PC-count settings and runs the existing analysis once for that group. The shard outputs are intermediate files and are not exposed as workflow outputs.

The workflow will then call `MergeLofPcEnrichment` with the arrays of shard result, summary, gene-PC QC, and analysis-QC files. The merge task will validate compatible metadata, concatenate shard rows, aggregate per-PC QC counters, recalculate global Benjamini–Hochberg FDR, and emit files using the current public filenames and schemas.

The coding-gene preparation remains independent of PC chunking. The final workflow outputs remain exactly:

- `results_tsv`
- `summary_json`
- `gene_pc_qc_tsv_gz`
- `analysis_qc_json`
- `protein_coding_genes_tsv`
- `protein_coding_genes_qc_json`

## Python interfaces

Add a header-only PC-count reader and a chunking operation that does not parse the full PC matrix. The chunking operation will share the existing `build_pc_grid` validation and adaptive-grid semantics.

Add a merge operation exposed through the CLI. It will accept arrays of shard files, validate that shard metadata agree on thresholds, carrier definitions, available PC count, and input provenance, then write the four merged analysis outputs. Result rows will be ordered by selected PC count and the existing threshold/carrier ordering. Global FDR will be calculated from all merged Fisher p-values.

## Resource behavior

Each analysis shard will retain the current analysis CPU, memory, retry, and dynamic disk inputs. The chunk-preparation task will use a small CPU/memory floor and a disk value sized from the localized PC file. The merge task will use a disk value sized from the localized shard outputs, with a configurable baseline tied to the existing analysis resource defaults.

## Error handling

- Empty PC grids are rejected.
- Non-positive `pc_counts_per_job` values are rejected.
- Invalid or non-consecutive PC headers are rejected before scatter.
- Missing, inconsistent, or duplicated shard metadata cause merge failure rather than silently producing a partial analysis.
- A shard failure propagates through the scatter and prevents merge.

## Testing

- Unit-test explicit and adaptive chunk construction, including the final short chunk and invalid chunk sizes.
- Unit-test merging with multiple synthetic shards, including global FDR values that differ from shard-local FDR values.
- Unit-test per-PC analysis-QC aggregation and gene-PC gzip concatenation.
- Extend WDL contract tests for the chunk task, scatter, merge task, new input, and unchanged public outputs.
- Extend the fixture runtime test to exercise multiple chunks and verify merged result/QC content.
- Run the complete Python suite and MiniWDL validation.
