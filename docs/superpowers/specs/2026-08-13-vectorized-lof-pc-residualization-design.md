# Vectorized LoF/PC Residualization

## Goal

Replace repeated per-gene least-squares fits with a complete-data NumPy
projection path, so the LoF/PC enrichment sweep is tractable for the 8,891
sample cohort and reports progress as each PC count finishes.

## Scope and data assumption

The supplied RNA/protein phenotype BED has no missing expression values for
the shared samples. The vectorized path therefore operates on complete rows.
The existing per-gene regression implementation remains as a compatibility
fallback when a future input contains missing values; its output semantics are
unchanged.

The existing explicit or adaptive `pc_counts` selection remains unchanged.
No PC-count cap is introduced by this change.

## Projection algorithm

For complete expression matrix `Y` with samples in rows and genes in columns,
and the first `k` principal-component score columns `P_k`, the intercept-plus-
PC residual is calculated after centering:

`R_k = Y_centered - U_k @ (U_k.T @ Y_centered)`

where `U_k` contains normalized PC score columns. This is the matrix form of
the existing ordinary-least-squares residualization for the PCA basis.

The implementation will:

1. Read coding BED rows into sample-by-gene blocks and retain their normalized
   Ensembl IDs and carrier metadata.
2. Center each gene once and compute PC projection coefficients by NumPy
   matrix multiplication.
3. Incrementally add PC contributions between requested PC counts rather than
   refitting each gene with `matrix_rank` and `lstsq`.
4. Standardize each residual vector with the existing population SD (`ddof=0`)
   and preserve the existing zero/non-finite residual exclusion rules.
5. Aggregate the Fisher contingency cells and gene-PC QC rows for each PC
   count without changing result or JSON schemas.

## Processing and logging

Processing becomes PC-major. A PC-count completion log is emitted immediately
after all coding genes for that count have been aggregated, with eligible-gene
and observation totals. Gene-block progress is emitted while that PC count is
running. This replaces the current gene-major arrangement, in which a single
gene must be fitted at every PC count before any PC-completion log can appear.

## Numerical and compatibility safeguards

- The fast path validates finite expression and PC inputs.
- It centers PC columns to preserve the intercept term and normalizes each
  nonzero PC column before projection.
- If no missing expression values are present, fast-path residuals must match
  `residualize_expression` on small fixtures within floating-point tolerance.
- Inputs with any missing expression retain the current per-gene fallback.
- The output TSV, QC gzip, and JSON schemas remain unchanged.

## Verification

Tests will compare vectorized and legacy residuals, Fisher cells, QC values,
and output rows for complete small fixtures across multiple PC counts. A
missing-expression fixture will verify fallback compatibility. A logging test
will verify that PC completion is emitted during PC-major processing rather
than only after the entire sweep.
