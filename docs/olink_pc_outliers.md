# Olink PC outliers in All of Us Workbench

`workflows/olink_pc_outliers.wdl` runs the Olink outlier analysis at several PC counts. It returns a residual Z-score matrix and outlier calls for every requested count, separately for `NPX` and `PCNormalizedNPX`. It does not return unscaled residual matrices.

## Inputs

Two workflow inputs are required:

| Input | Type | Meaning |
| --- | --- | --- |
| `npx_file` | `File` | Plain or gzip-compressed TSV containing `Assay`, `SampleID`, `NPX`, and `PCNormalizedNPX`. |
| `pc_counts` | `Array[Int]` | Nonempty list of distinct, nonnegative PC counts, such as `[0, 10, 20, 50, 3499]`. |

Use a real `gs://` input URI from the Workbench resource details. The mounted `/home/rstudio/workspace/...` filename in an RStudio session is not a cloud URI for a batch VM. The workflow service account needs read access to the input object and write access to the destination bucket.

An example is provided in `examples/olink_pc_outliers.inputs.json`. Replace the input URI before use. The example uses the default Docker image and includes 3,499 PCs to retain the original setting. A count that leaves no residual component stops preparation with an error; it is never reduced automatically.

Optional inputs include:

| Input | Default | Meaning |
| --- | --- | --- |
| `docker_image` | `ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink:main` | Image used by all tasks. Override it with a fixed tag, digest, or workspace-accessible registry address when needed. |
| `excluded_assays` | `["GBP1", "MAP2K1"]` | Assays excluded before missingness checks. An empty array excludes none. |
| `strip_sample_suffix` | `true` | Remove the suffix from the first underscore in each sample ID, as in the original R code. |
| `connectivity_cutoff` | `-3.0` | Remove samples with connectivity Z scores strictly below this value. |
| `protein_cutoff` | `3.0` | Flag protein residual Z scores strictly below `-cutoff` or above `cutoff`. |
| `block_size` | `500` | Number of sample profiles in each bicor block. |
| `prepare_cpu`, `prepare_memory_gb`, `prepare_disk_gb` | `8`, `128`, `500` | Resources for each of the two preparation tasks. |
| `adjust_cpu`, `adjust_memory_gb`, `adjust_disk_gb` | `8`, `64`, `200` | Resources for each PC-count adjustment task. |
| `preemptible`, `max_retries` | `0`, `1` | GCP task settings. Preparation can be long, so preemption is disabled by default. |
| `zones` | Four `us-central1` zones | GCP zones for task placement. |

Preparation disk is raised when `ceil(12 × compressed input GiB + 50)` exceeds the supplied baseline. Adjustment disk is raised when `ceil(6 × prepared-state GiB + 20)` exceeds its baseline. These are storage allowances, not estimates of the required RAM. Dense PCA and the input table remain memory-intensive.

## Method

For each NPX type, preparation:

1. Reads the sample and assay IDs as strings. Leading zeros and existing `X` characters are preserved.
2. Applies the assay exclusions and optional suffix removal. Duplicate assay/sample keys and distinct IDs that would be merged are rejected.
3. Removes samples with any missing or nonfinite value for this NPX type.
4. Removes constant assays and standardizes each assay across retained samples.
5. Calculates signed bicor adjacency, `0.5 + 0.5 × bicor`, and sums connections without the diagonal. Connectivity Z scores use sample SD.
6. Removes samples below the connectivity cutoff, checks constants again, and rescales the retained data.
7. Calculates one exact PCA up to the largest requested PC count.

Each adjustment task subtracts the fitted values from the requested prefix of that PCA. It centers each assay's residuals and divides by population SD, as in the supplied R code. Near-zero residual SD produces `NA` scores. PC=0 performs residual standardization without subtracting a PC component.

The sample QC is shared across PC counts because it occurs before PC removal. The two NPX types can have different retained sample and assay sets. Protein outliers are individual assay/sample measurements, rather than whole samples. The outlier threshold does not remove cells from the complete Z-score matrix.

For a grid of `k` counts, the workflow runs two preparation tasks and `2k` adjustment tasks, plus grid validation and summary collection. Shared preparation avoids repeated PCA; every adjustment still reads a prepared state and writes a complete matrix. Large grids increase storage and total compute use.

## Outputs

`results` is an array of typed records. Each record contains `npx_type`, `pcs_removed`, and explicit File outputs for its matrix, protein calls, sample calls, undefined assays, and summary. These records link each matrix to the matching calls without relying on array positions.

Examples for 3,499 PCs:

| File | Contents |
| --- | --- |
| `protein_z_scores.NPX.PCs3499.tsv.gz` | Complete NPX residual Z-score matrix. |
| `protein_z_scores.PCNormalizedNPX.PCs3499.tsv.gz` | Complete PCNormalizedNPX residual Z-score matrix. |
| `protein_outliers.NPX.PCs3499.tsv.gz` | Protein outlier calls with NPX type, PC count, assay, sample ID, Z score, and direction. |
| `sample_outliers.NPX.PCs3499.tsv` | Sample connectivity outliers. The PC suffix identifies the run setting. |
| `undefined_residual_assays.NPX.PCs3499.tsv` | Assays with undefined residual Z scores. |
| `summary.NPX.PCs3499.tsv` | Cutoffs, counts, and output filenames for this result. |

Matrix rows are assays. The first column is `Assay`; remaining headers are sample IDs. Missing residual scores are `NA`.

The workflow also returns per-type sample QC, assay exclusions, preparation summaries, R session information, the PC plan, and a combined `summary.tsv`. Prepared RDS files are intermediate File dependencies. They are not final matrix outputs. All participant-level outputs remain in the workspace destination bucket.

## All of Us Workbench setup

All of Us uses Google Batch for batch execution. Its documentation describes restricted internet access from batch VMs and registry-specific image access. This image contains the R script and packages, so tasks do not install software or fetch code from GitHub. See [batch processing](https://support.researchallofus.org/hc/en-us/articles/4692418691732-Overview-of-Batch-Processing) and [Docker images in Workbench](https://support.researchallofus.org/hc/en-us/articles/21179878475028-Using-Docker-Images-on-the-Workbench).

The dedicated GitHub Actions job publishes the tested Linux AMD64 image as:

```text
ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink:sha-<full-commit-sha>
```

The workflow defaults to `ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink:main`. The `main` tag follows the latest image published by the image workflow. Set the optional `docker_image` input to a fixed tag or digest to retain a specific image version. If the workspace cannot pull from GHCR, set it to a workspace-accessible Google Artifact Registry copy. The documented All of Us Docker Hub proxy is specific to Docker Hub; do not prepend that proxy to a GHCR address.

An administrator with access to the source image and the destination registry can mirror the image without building it again:

```bash
# Run only in an environment authorized to access both registries.
source_image='ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink:sha-REPLACE_WITH_COMMIT'
target_image='us-central1-docker.pkg.dev/YOUR_PROJECT/YOUR_REPOSITORY/olink-outliers:sha-REPLACE_WITH_COMMIT'

gcloud auth configure-docker us-central1-docker.pkg.dev
docker pull "$source_image"
docker tag "$source_image" "$target_image"
docker push "$target_image"
```

The destination repository must exist and the workspace runtime identity must be able to read it. Record the resulting digest and set `docker_image` to that address. These commands transfer the software image only. They do not require access to NPX files or results. No registry transfer or All of Us job is submitted by the repository tests.

### Current Workbench interface

In the Workflows interface, add `workflows/olink_pc_outliers.wdl` from GitHub or upload this single WDL to a workspace bucket. Create a job, supply the two required inputs, and select a workspace bucket and folder for outputs. Use the actual storage URI for `npx_file`. The image defaults to `:main`; no image entry is required. If you enter an override in the form, use plain text without surrounding quotes. JSON input files require JSON string quotes. Review the job settings before starting it. See [Workbench Cromwell instructions](https://support.workbench.verily.com/docs/guides/workflows/cromwell/).

### Earlier Cromwell application interface

For workspaces with a Cromwell cloud application, use the official setup snippet and Cromshell instructions. After preparing the input JSON, the documented submission command is:

```bash
cromshell-alpha submit workflows/olink_pc_outliers.wdl examples/olink_pc_outliers.inputs.json
```

Submission starts billable cloud work. The command is provided for the user to run after configuration. See [All of Us Cromwell instructions](https://support.researchallofus.org/hc/en-us/articles/14428263737620-How-to-use-Cromwell-in-the-All-of-Us-Researcher-Workbench).

## Validation

The R tests compare PCA-prefix residuals with an independent QR regression. They cover both metrics, PC=0, sample QC, independent missingness, constants, ID preservation, rank limits, CLI failures, and output correspondence.

The WDL tests begin with cloud File values, replace them with localized paths, and execute rendered commands. They include spaces, quotes, and shell characters. Production CLI checks reject unresolved cloud URIs before opening files. An AST check rejects workflow-scope `write_lines`, `write_tsv`, `write_map`, and `write_json`, including nested expressions and imported workflows.

GitHub Actions builds the micromamba image, runs the R tests without network access, and runs the complete WDL on synthetic data for both NPX types at PCs 0 and 2. It publishes that tested image on `main`. Image rebuilds are triggered by the bundled R script or its build definitions; documentation, WDL, and test edits do not rebuild the image automatically.

Local checks:

```bash
Rscript tests/olink/smoke_test.R
python -m pytest -q
miniwdl check workflows/olink_pc_outliers.wdl
python tools/check_workflow_file_writes.py workflows/olink_pc_outliers.wdl
```

No local Docker build is required. Neither a local test nor a GitHub Actions synthetic run establishes All of Us or Terra validation. The complete workflow has not been run on All of Us Workbench or Terra, and it has not been run on the controlled NPX dataset.
