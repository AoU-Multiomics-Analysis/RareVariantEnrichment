# Multi-matrix enrichment and outlier intersections

## Intended result

Run the existing enrichment and PC-selection workflow separately for each named molecular-phenotype matrix. Keep each original matrix and all per-matrix outputs, including the selected-PC Z-score matrix. Compare matching gene and sample IDs across datasets and export one indicator matrix per dataset combination and threshold.

The user confirmed that every dataset in an intersection uses the same threshold. Intersections are inclusive by default: a three-way outlier also appears in each applicable pairwise intersection. The user confirmed that each intersection tests only its participating datasets. Use the existing lower-tail rule, Z <= threshold.

## Inputs and rules

- WDL 1.0 for Terra managed Cromwell.
- A TSV manifest with ome_name, phenotype_bed, principal_components_tsv, and optional additional_covariates_tsv. Each ome has its own PC table and at most one file of fixed covariates. At least two datasets, with unique safe names.
- Shared LoF carrier table, gene annotation, and enrichment/PC-selection settings. Each dataset selects its own PC count. The automatic PC grid adapts to each PC table.
- Intersection thresholds default to -2, -3, -4, -5, -6 and are independent of the thresholds used for PC selection.
- Generate every combination of at least two datasets, at each threshold. The number of files is T * (2^N - N - 1).
- Align on shared genes and shared sample IDs separately for each combination. Preserve the first participating dataset's order. Do not use matrix position as identity.
- Cell 1: all participating Z scores are finite and <= threshold. Cell 0: all are finite and at least one exceeds threshold. Cell NA: any required Z score is missing, even if another score exceeds threshold.
- Preserve all genes in the per-dataset Z-score exports, including noncoding genes. Reuse the current normalization and duplicate-feature collapse there. Reject duplicate IDs in exported Z-score inputs to the intersection step.
- Empty gene/sample overlap produces an empty matrix with appropriate headers and zero observation counts, with the overlap recorded in QC.

## Implementation

Import rare_variant_enrichment.wdl and scatter its workflow over validated manifest rows. A preparation task checks names, path fields, and thresholds, and returns a normalized metadata TSV. The workflow reads that TSV and declares the BED, PC table, and optional covariate paths as File values before the imported workflow call. Preserve these File values inside results and file arrays until command rendering. Relative manifest paths are rejected; use GCS URIs or absolute local paths. There is no separate ome-specific fixed-covariate field. The intersection task creates its newline-delimited local file list during command rendering; JSON is used only for structured output summaries.

Load each compressed Z-score matrix once into a temporary SQLite table with one numeric BLOB per gene. This supports mismatched row order without holding all matrices in memory. Join genes by ID, align sample columns by ID, and write all threshold outputs for one dataset combination in one pass. Remove the temporary database after use.

Return typed per-dataset results, the Z-score file array, intersection file array, a TSV index that maps filenames to datasets/thresholds and counts, and a JSON summary. Index paths are output metadata, not task input paths.

## Verification and limits

Test exact threshold equality, mixed signs, missing data, reordered identifiers, pairwise-specific overlap, three-way membership, empty overlap, duplicate identifiers, malformed rows, invalid thresholds, unreadable files, and cloud URI rejection. Test WDL file localization with simulated cloud-to-local substitution, safe shell quoting, and static workflow-scope file-writing checks. Extend the existing GitHub Actions fixture smoke test. Do not build a local Docker image or submit cloud jobs. Do not claim Terra validation without an actual Terra run.
