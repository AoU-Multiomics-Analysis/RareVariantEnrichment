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
| `docker_image` | `docker.io/epadhi/rarevariantenrichment-olink:root-no-entrypoint` | Image used by all tasks. This is the image that passed the user's minimal Workbench workflow. Override the address to select another image version or registry. |
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

The dedicated GitHub Actions job publishes the tested Linux AMD64 image to both registries:

```text
ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink:sha-<full-commit-sha>
docker.io/epadhi/rarevariantenrichment-olink:sha-<full-commit-sha>
```

Both registries also receive the `main` tag after the tests pass. See the Docker Hub setup below before the first publication.

The image runs as root and has no runtime entrypoint. `Rscript` is on `PATH`, so tasks can call R directly without an environment activation script. Micromamba installs the pinned R environment during the image build.

The WDL defaults to `docker.io/epadhi/rarevariantenrichment-olink:root-no-entrypoint`. The user reported that a minimal Workbench workflow passed with this image. Its digest is `sha256:bd3a69b1aad32966959fb4eb68f883afe48fd1a6ba64a29530a710ef659e97ef`. All four tasks use this workflow input. Set the optional `docker_image` input to select another image version or registry. The All of Us Docker guide documents Docker Hub and Google Artifact Registry image access. It does not require all WDL images to be copied into a workspace GAR repository.

A Workbench Network policy restricts Batch VM internet access. Image access also depends on the registry route used by the workflow service. A generic task failure does not prove that a registry pull failed. Interactive apps and GitHub Actions have different network access, so a successful image pull there does not prove that a Batch VM can pull it. See [Workbench network policy and container access](https://support.workbench.verily.com/docs/guides/governance/network_policy/).

### Docker Hub publication

1. Create a public Docker Hub repository named `epadhi/rarevariantenrichment-olink`.
2. Create a Docker personal access token with **Read & Write** permission. See [Docker access tokens](https://docs.docker.com/security/access-tokens/).
3. In this GitHub repository, open **Settings → Secrets and variables → Actions**. Add two repository secrets: `DOCKERHUB_USERNAME` with value `epadhi`, and `DOCKERHUB_TOKEN` with the token value. See [GitHub repository secrets](https://docs.github.com/en/actions/how-tos/write-workflows/choose-what-workflows-do/use-secrets#creating-secrets-for-a-repository). Do not commit the token to a file or paste it into an issue.
4. After these workflow changes reach `main`, open **Actions → Olink Docker Hub Publish → Run workflow**. Select `main`. Enter the full commit SHA of the Olink image that passed Image CI. Use a run that built the current Dockerfile. There is no default SHA, so this step cannot silently select the previous startup configuration.
5. Check that the run succeeds. It copies the selected image with `main` and `sha-<source_sha>` tags and checks that the published image ID matches the source. It does not build an image or start an All of Us job.

Future Olink image builds publish the same tested image to GHCR and Docker Hub. Changes to documentation, WDL files, tests, or publishing configuration do not trigger an image rebuild. Use the manual copy workflow to retry Docker Hub publication if GHCR publication succeeded but Docker Hub publication failed.

### Optional copy into Google Artifact Registry

Run the following from an interactive environment that can reach both registries and has a working Docker daemon. The destination repository must already exist. First set `OLINK_TESTED_DIGEST` to the `sha256:...` digest from the successful Image CI run that built the current Dockerfile. The earlier image digest uses the previous startup configuration. No image is built by this copy:

```bash
(
  set -euo pipefail
  : "${OLINK_TESTED_DIGEST:?Set OLINK_TESTED_DIGEST to the tested image digest from Image CI}"
  olink_source="ghcr.io/aou-multiomics-analysis/rarevariantenrichment-olink@${OLINK_TESTED_DIGEST}"
  olink_target='us-central1-docker.pkg.dev/wb-perky-tomato-1843/olink-images/rarevariantenrichment-olink:main'

  docker info >/dev/null
  gcloud artifacts repositories describe olink-images \
    --project=wb-perky-tomato-1843 --location=us-central1
  gcloud auth configure-docker us-central1-docker.pkg.dev --quiet

  docker pull --platform linux/amd64 "$olink_source"
  docker tag "$olink_source" "$olink_target"
  docker push "$olink_target"

  gcloud artifacts docker images describe "$olink_target" \
    --project=wb-perky-tomato-1843 \
    --format='value(image_summary.digest)'
)
```

The account that copies the image needs write access to the repository. The workflow runtime identity needs read access. A permission error must be resolved before retrying the copy. These commands transfer the software image only; they do not access NPX files or results. No registry transfer or All of Us job is submitted by the repository tests. A future GHCR image update does not update this GAR copy automatically; repeat the copy with the new tested source digest.

### Current Workbench interface

In the Workflows interface, add `workflows/olink_pc_outliers.wdl` from GitHub or upload this single WDL to a workspace bucket. Create a job, supply the two required inputs, and select a workspace bucket and folder for outputs. Use the actual storage URI for `npx_file`. The image defaults to the Docker Hub `:root-no-entrypoint` tag above. Clear any saved `docker_image` override to use this default, or set it to the same working tag. If you enter an override in the form, use plain text without surrounding quotes. JSON input files require JSON string quotes. Review the job settings before starting it. See [Workbench Cromwell instructions](https://support.workbench.verily.com/docs/guides/workflows/cromwell/).

### Earlier Cromwell application interface

For workspaces with a Cromwell cloud application, use the official setup snippet and Cromshell instructions. After preparing the input JSON, the documented submission command is:

```bash
cromshell-alpha submit workflows/olink_pc_outliers.wdl examples/olink_pc_outliers.inputs.json
```

Submission starts billable cloud work. The command is provided for the user to run after configuration. See [All of Us Cromwell instructions](https://support.researchallofus.org/hc/en-us/articles/14428263737620-How-to-use-Cromwell-in-the-All-of-Us-Researcher-Workbench).

## Validation

The R tests compare PCA-prefix residuals with an independent QR regression. They cover both metrics, PC=0, sample QC, independent missingness, constants, ID preservation, rank limits, CLI failures, and output correspondence.

The WDL tests begin with cloud File values, replace them with localized paths, and execute rendered commands. They include spaces, quotes, and shell characters. Production CLI checks reject unresolved cloud URIs before opening files. An AST check rejects workflow-scope `write_lines`, `write_tsv`, `write_map`, and `write_json`, including nested expressions and imported workflows.

GitHub Actions builds the micromamba image, checks that it runs as root with no entrypoint, and runs validation without an activation script. It then runs the R tests without network access and the complete WDL on synthetic data for both NPX types at PCs 0 and 2. It publishes that tested image on `main`. Image rebuilds are triggered by the bundled R script or its build definitions; documentation, WDL, and test edits do not rebuild the image automatically.

The separate `Olink Published Image Startup` workflow checks an existing image without rebuilding it. It runs `validate` with a `/bin/bash` entrypoint, the image's default user, a writable task directory, one CPU, and a hard 2 GB memory limit without extra swap. It checks both validation output files. This test does not reproduce the Workbench network policy or the complete Batch VM.

Local checks:

```bash
Rscript tests/olink/smoke_test.R
python -m pytest -q
miniwdl check workflows/olink_pc_outliers.wdl
python tools/check_workflow_file_writes.py workflows/olink_pc_outliers.wdl
```

No local Docker build is required. Neither a local test nor a GitHub Actions synthetic run establishes All of Us or Terra validation. No successful complete workflow run on All of Us Workbench or Terra has been recorded, including on the controlled NPX dataset.
