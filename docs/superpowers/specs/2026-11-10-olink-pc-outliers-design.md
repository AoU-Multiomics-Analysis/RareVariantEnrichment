# Olink PC outlier workflow

## Intent

Convert the agreed R analysis into WDL 1.0 for the All of Us Researcher Workbench and Terra-managed Cromwell. Accept several PC counts. For each count, produce the residual Z-score matrix and protein outlier calls for NPX and PCNormalizedNPX. The user selected Z-score outputs only.

## Method and constraints

Preserve signed bicor connectivity, sample Z < -3, assay scaling, exact PCA projection, population residual SD, and strict protein Z < -3 or Z > 3. Keep configurable exclusions with GBP1 and MAP2K1 as defaults. Preserve string sample IDs and reject ambiguous normalized IDs. Keep separate missingness filters for each NPX column. Do not change PC counts automatically.

Use typed WDL File inputs, safe named CLI arguments, command-time localization, and task-local text lists. Do not create files at workflow scope or wrap script arguments in JSON. Log each command stage. Use a pinned micromamba image with conda-forge and bioconda. Bake scripts and packages into the image; batch tasks must not download code or packages. Do not submit cloud jobs without approval.

## Approaches

Running the complete R analysis independently for every PC count repeats input loading, sample QC, and PCA. The implementation will instead run preparation once per NPX type and calculate one exact PCA up to the largest requested count. Each PC task uses the corresponding score/loading prefix. This preserves the method while avoiding repeated PCA. The choice is an implementation decision; the requested scientific outputs are unchanged.

## Components

- `scripts/olink_outliers.R`: validate, prepare, adjust, and collect CLI stages.
- `workflows/olink_pc_outliers.wdl`: validate the grid, scatter over both NPX types, then scatter over PC counts. Collect summary files.
- `envs/Dockerfile.olink` and environment definition: independent analysis image.
- Synthetic numerical, CLI, localization, workflow-scope, and container checks.
- AoU instructions and example inputs with workspace-accessible Cloud Storage paths and a required accessible image reference.

## Output contract

Every result records NPX type and PC count. Matrices have Assay rows and sample-ID columns. Protein call files contain only threshold outliers. Sample QC occurs before PC removal and is shared across counts. Intermediate prepared RDS files are explicit task File dependencies; they are not JSON path wrappers. Final summaries describe files as metadata. WDL outputs carry the actual output locations.

## Delivery and validation

Push an isolated branch to the requested repository and use its existing PR checks. Build and smoke-test images in GitHub Actions, not locally. Image rebuild filters include only bundled production code and its build definitions. Run local R tests, Python tests, miniwdl, and the static workflow-scope regression. Add a full synthetic WDL run in Actions. A passing local or Actions test does not establish Terra or AoU validation; report that distinction.
