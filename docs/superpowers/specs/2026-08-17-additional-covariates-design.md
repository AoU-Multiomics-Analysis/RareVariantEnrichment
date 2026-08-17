# Additional Covariates for LoF/PC Enrichment Design

## Goal

Allow the LoF/PC enrichment analysis to adjust for an optional sample-by-covariate table, such as genotype principal components, while preserving the current behavior when no additional covariates are supplied.

## Inputs

The WDL and command-line interface will accept an optional tab-delimited additional-covariate file. The file must contain one sample identifier column named `sample_id` or `ID`, plus one or more numeric covariate columns. The identifier column may occur anywhere in the header. Every non-identifier column is included as a covariate, so a table with columns `GENETICPC1` through `GENETICPC67` and `sample_id` is supported without manually listing columns.

All sample identifiers are treated as stripped strings throughout the analysis. Readers will not infer numeric types or cast identifiers to numbers. Empty identifiers, duplicate identifiers, missing values, non-numeric values, and non-finite covariate values are rejected with actionable errors.

## Sample alignment

The analysis sample universe is the intersection of sample identifiers present in the phenotype BED, molecular-PC matrix, and optional additional-covariate matrix. If the optional matrix is absent, the intersection remains the existing BED/phenotype-PC intersection.

The LoF carrier table is not included in this intersection because it is a sparse carrier representation: the absence of a sample/gene row indicates a non-carrier and must not remove that sample from the analysis.

The analysis QC output will report the sample counts for each input and the final shared sample count. It will also report the number and names of additional covariate columns when supplied.

## Residualization model

For each requested molecular-PC count `k`, each gene’s molecular phenotype is residualized using an intercept, all supplied additional covariates, and the first `k` molecular PCs:

```text
phenotype ~ intercept + additional_covariates + molecular_PC1 + ... + molecular_PCk
```

When `k = 0`, the model still includes the intercept and all additional covariates. When no additional covariate file is supplied, the model reduces to the current intercept-plus-molecular-PC design.

The existing complete-data projection and PC-sweep behavior will be extended to include the additional covariate columns in the fixed part of the design. Rank, degrees-of-freedom, residual-scale, and exclusion checks remain unchanged.

## Outputs and provenance

The analysis QC and summary JSON outputs will record:

- whether an additional-covariate file was supplied;
- the number and names of additional covariate columns;
- the number of samples in the additional-covariate file;
- the final BED/molecular-PC/additional-covariate intersection count; and
- the additional-covariate input path in provenance metadata when supplied.

The summary residualization description will state that the design includes all additional covariates followed by the first `k` molecular PCs. Existing result, gene-PC QC, selection, and plot outputs remain unchanged.

## Error handling

The workflow will fail before enrichment calculations if the covariate header has no recognized sample-ID column, contains duplicate column names, contains no covariate columns, has duplicate sample IDs, or contains malformed/non-finite covariate values. If the three required sample-bearing inputs have no shared samples, the analysis will fail with a clear intersection error.

## Compatibility

The additional-covariate input is optional. Existing WDL input files and command invocations remain valid and produce the existing model when the new input is omitted. The WDL disk-sizing expressions will include the optional covariate file when present.

The workflow will add an integer `pc_preemptible` runtime input, defaulting to `2`, to control the number of preemptible attempts for the scattered PC-fitting/enrichment tasks. The value will be passed to the `preemptible` runtime attribute for each `CalculateLofPcEnrichment` shard. Other tasks retain their existing retry/runtime behavior.

## Testing

Tests will cover:

1. Reading a `GENETICPC...` table with `sample_id` as the final column.
2. String-preserving sample-ID alignment across all input readers.
3. Validation of malformed covariate tables.
4. Residualization with fixed additional covariates at `k = 0` and with molecular PCs included.
5. End-to-end sample intersection and QC/provenance output.
6. CLI dispatch and optional WDL input/template compatibility.
7. WDL propagation of the configurable PC-fitting `preemptible` runtime value.
