# Multi-omics outlier intersections implementation plan

**Goal:** Run enrichment on named matrices and export inclusive outlier intersections at common thresholds.
**Architecture:** Import and scatter the existing WDL; combine its Z-score outputs in a separate task using disk-backed gene alignment.
**Tech stack:** WDL 1.0, Python, NumPy, SQLite, pytest, miniwdl.
**Spec:** ../specs/2026-09-21-multiomics-outlier-intersections-design.md

## Constraints

Terra managed Cromwell; typed File values until command rendering; no workflow-scope write functions; task-local newline-delimited file lists; existing negative-tail semantics; no local Docker builds or cloud submissions.

## Review focus

Mismatched identifiers must align by ID. NA must never silently become 0. Triple outliers must appear in pairwise outputs. Empty overlap must have explicit QC. Local file lists must contain localized paths, including optional dataset covariates.

## Tasks

- [x] Add tests/test_multiomics.py with literal expected matrices for three datasets at -2 and -3; assert 8 output matrices, per-combination order, exact equality at thresholds, NA propagation, empty overlap, and input errors. Run before implementation and confirm the missing command fails.
- [x] Add src/rare_variant_enrichment/multiomics.py. Implement validate_multiomics_inputs(names, thresholds) and build_outlier_intersections(matrix_paths, names, thresholds, output_directory, manifest_output, summary_output). Read and validate matrices, store gene vectors in temporary SQLite tables, join by gene ID, align samples, and write matrices plus counts. Add CLI commands validate-multiomics-inputs and multiomics-intersections with named list-file arguments. Run the unit and CLI tests.
- [x] Add workflows/multiomics_outliers.wdl. Define OmicsDataset and OmicsResult structs, validation and intersection tasks, a scatter of the imported RareVariantEnrichment workflow, and typed outputs. Add an example inputs file and Dockstore registration. Test WDL contracts and cloud-to-local command rendering for both new tasks and struct fields; prohibit workflow-scope write functions.
- [x] Add a Docker-backed fixture test to the existing GitHub Actions test suite. It must run two datasets through the wrapper and check per-dataset PC counts, matrices, and intersection values. Document input/output mapping, lower-tail and missing-value rules, inclusive combinations, output counts, and Terra validation status.
- [x] Run the relevant tests, full suite, miniwdl checks, and git diff --check. Request one independent code review, resolve material findings, and report results.

## Execution record

Implemented in an isolated checkout from origin/main. No cloud jobs were submitted and no local Docker image was built.

- Final full suite: 333 passed, 10 skipped (Docker/htslib prerequisites unavailable).
- WDL validation passed for all three workflows; example inputs passed typed WDL input validation.
- Independent review found a smoke-test threshold-order error and a repeated SQLite sort. Corrected the smoke expectation and verified it with real fixture exports. Made the stored row position the integer primary key; the query-plan regression test failed before the change and passed after it.
- Twenty-six focused tests cover the new calculator, CLI, integration, and WDL command localization.
- Complete Terra execution remains untested.

## Manifest interface update

- [x] Replace the dataset-array input with a TSV manifest: ome_name, phenotype_bed, principal_components_tsv, optional additional_covariates_tsv.
- [x] Validate and normalize metadata inside a task. Declare each referenced path as a WDL File before downstream localization. Keep fixed covariates optional and separate from the PC sweep.
- [x] Test missing fields, invalid paths, optional covariates, safe quoting, literal TSV output, File conversion, and cloud-to-local path substitution. Update the GitHub Actions workflow smoke fixture to use the manifest.
- [x] Run the full suite and WDL checks after the manifest update: 344 passed, 10 skipped. All three WDL files and the example input types passed validation. The three manifest String-to-File warnings are expected at the explicit localization declarations.

## PC-sweep covariate localization

PR preparation found that the existing PC-sweep task converted the optional covariate File to a String before localization. A cloud-path regression failed with an unresolved GCS URI. Move optional argument construction to command rendering, with safe quoting, as in the selected-PC export task. The regression covers both present and absent covariates and executes the rendered command. The regression design remains unchanged.

Final PR checks after the covariate path fix: 346 passed, 10 skipped. WDL validation passed with expected manifest String-to-File warnings. The complete workflow has not run on Terra.

## Haplo matrix extension

The user requested this output in PR #26 and confirmed that the phenotype BED contains log2-CPM. Export one haplo matrix per ome at its selected PC count, with the same aligned genes and samples as its Z-score matrix. Match UnderlierPrevelance: fit phenotype PCs only, compare adjusted values to the gene mean with a strict configurable drop (default 1), and retain missing values as NA. Fixed covariates still apply to Z scores and cohort alignment. Include haplo counts and rules in the export summary. Test threshold equality, constant genes, missing values, selected PCs, fixed-covariate separation, and an independent least-squares reference.

The first GitHub smoke run found miniwdl's local manifest file-access restriction. Enable that access only in the trusted local fixture test; task File localization in the production WDL remains unchanged. Verify the fix in GitHub Actions without a local Docker build or cloud submission.

Haplo verification: 358 passed, 10 skipped locally; WDL and example input checks passed. Independent review found and resolved a constant-decimal rounding edge at drop=0, with a failing-then-passing regression. Review found no remaining material defects.

GitHub ran all workflow tasks successfully, including haplo value checks. The remaining test failure compared output alias paths, which miniwdl materializes in separate directories. The test now compares decompressed matrix contents in manifest order.

## Expression-only haplo export

The user restricted haplo output to jobs labelled expression. Forward each manifest ome_name to the child workflow and emit haplo only for an exact expression match. The standalone workflow defaults ome_name to an empty string and skips haplo unless explicitly labelled expression. Keep File? haplo outputs and Array[File?] wrapper outputs so null entries preserve manifest order. Other jobs continue to export Z scores without a haplo file or summary. Test label matching, optional outputs, a mixed-ome workflow, and a workflow with no expression row.

Expression-only local verification: 366 passed, 11 skipped. WDL and both example input types passed. Independent review found no actionable defects.
