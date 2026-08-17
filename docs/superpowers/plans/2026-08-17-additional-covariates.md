# Additional Covariates for LoF/PC Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add optional genotype-PC covariate adjustment and configurable preemptible attempts to the LoF/PC enrichment workflow while preserving existing behavior when no covariate file is provided.

**Architecture:** Add a strict reader for sample-by-covariate TSVs, normalize all LoF/PC sample IDs as strings, and align BED, molecular-PC, and optional covariate matrices before residualization. Extend both the legacy and vectorized residualization paths so every model uses an intercept, all additional covariates, and the requested molecular-PC prefix. Wire the optional file and a `pc_preemptible` runtime setting through the WDL scatter and record the resolved design and sample alignment in QC/provenance outputs.

**Tech Stack:** Python 3, NumPy, pytest, WDL 1.0, MiniWDL, Cromwell Google runtime attributes.

## Global Constraints

- The additional covariate file is optional; omitting it must retain the current intercept-plus-molecular-PC behavior.
- Every sample identifier is loaded, stripped, and compared as a string; no numeric inference or numeric casting is allowed.
- The LoF carrier table remains a sparse lookup and is excluded from the sample intersection.
- Every non-ID column in the additional covariate file is included as a finite numeric covariate.
- The default `pc_preemptible` value is `2`, and it applies to `CalculateLofPcEnrichment` scatter shards only.
- Use the existing NumPy implementation and do not add pandas, Hail, or another matrix dependency.
- Preserve the existing result TSV and gene-PC QC schemas; add metadata only to JSON QC/provenance outputs.

---

### Task 1: Add failing tests for covariate-table parsing and string sample IDs

**Files:**
- Modify: `tests/test_lof_pc.py`
- Modify: `src/rare_variant_enrichment/lof_pc.py` only after the failing tests are observed

**Interfaces:**
- Add `read_additional_covariates(path: Path) -> AdditionalCovariates`.
- Add `normalize_sample_id(value: object, context: str) -> str` and use it in the LoF/PC readers.
- `AdditionalCovariates` exposes `sample_ids: tuple[str, ...]`, `values: np.ndarray`, `names: tuple[str, ...]`, `sample_count: int`, and `covariate_count: int`.

- [ ] **Step 1: Write the failing reader and identifier tests**

Add tests with these exact cases:

```python
def test_read_additional_covariates_accepts_sample_id_as_final_column(tmp_path):
    path = tmp_path / "genetic-pcs.tsv"
    path.write_text(
        "GENETICPC1\tGENETICPC2\tsample_id\n"
        "0.1\t-0.2\t001\n"
        "0.3\t0.4\t1000291\n"
    )

    matrix = lof_pc_module().read_additional_covariates(path)

    assert matrix.sample_ids == ("001", "1000291")
    assert matrix.names == ("GENETICPC1", "GENETICPC2")
    np.testing.assert_allclose(matrix.values, [[0.1, -0.2], [0.3, 0.4]])


def test_read_additional_covariates_accepts_id_column_in_any_position(tmp_path):
    path = tmp_path / "covariates.tsv"
    path.write_text("ID\tGENETICPC1\n001\t0.1\n")

    matrix = lof_pc_module().read_additional_covariates(path)

    assert matrix.sample_ids == ("001",)
    assert matrix.names == ("GENETICPC1",)


@pytest.mark.parametrize(
    "text",
    [
        "GENETICPC1\tsample_id\n0.1\tS1\n0.2\tS1\n",
        "GENETICPC1\tother\n0.1\tS1\n",
        "sample_id\tGENETICPC1\nS1\tNA\n",
        "sample_id\tGENETICPC1\nS1\tinf\n",
        "sample_id\tGENETICPC1\n\t0.1\n",
    ],
)
def test_read_additional_covariates_rejects_invalid_schema_or_values(tmp_path, text):
    path = tmp_path / "invalid.tsv"
    path.write_text(text)

    with pytest.raises(ValueError):
        lof_pc_module().read_additional_covariates(path)


def test_lof_pc_sample_id_readers_preserve_numeric_looking_ids_as_strings(tmp_path):
    pcs = tmp_path / "pcs.tsv"
    pcs.write_text("ID\tPC1\n001\t0.1\n1000291\t0.2\n")
    assert lof_pc_module().read_principal_components(pcs).sample_ids == (
        "001",
        "1000291",
    )
```

The invalid-schema parametrization must also include duplicate header names and a table with no non-ID columns. Add a carrier-table assertion that a numeric-looking `sample_id` remains the exact string used in its carrier pair.

- [ ] **Step 2: Run the focused tests and verify the expected RED failure**

Run:

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py -k 'additional_covariate or numeric_looking'
```

Expected: collection succeeds and the new tests fail because the reader, dataclass, and shared string-normalization helper do not yet exist.

- [ ] **Step 3: Write the minimal reader and shared string normalization**

In `src/rare_variant_enrichment/lof_pc.py`:

1. Add an `AdditionalCovariates` frozen dataclass containing ordered string sample IDs, a two-dimensional finite NumPy matrix, and ordered covariate names.
2. Add `normalize_sample_id(value, context)` that converts the input to `str`, strips surrounding whitespace, and rejects an empty result. Do not parse integers or normalize values such as `001` to `1`.
3. Add `read_additional_covariates` that checks a non-empty unique header, finds exactly one `sample_id` or `ID` column, requires at least one remaining column, rejects blank rows and wrong field counts, rejects duplicate sample IDs, parses every covariate with `float`, and rejects non-finite values.
4. Apply `normalize_sample_id` to BED sample IDs used by `calculate_lof_pc_enrichment`, molecular-PC sample IDs, and LoF carrier sample IDs. Keep feature and gene normalization separate.

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run the same focused command from Step 2. Expected: all new reader and string-ID tests pass. Then run:

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py -k 'read_principal_components or read_lof_carriers or additional_covariate'
```

Expected: the existing reader tests and the new covariate tests pass without changing the legacy reader contracts.

- [ ] **Step 5: Commit the reader boundary**

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: read optional covariate matrices"
```

### Task 2: Extend residualization with fixed additional covariates

**Files:**
- Modify: `tests/test_lof_pc.py`
- Modify: `src/rare_variant_enrichment/lof_pc.py`

**Interfaces:**
- Extend `residualize_expression(expression, principal_components, pc_count, additional_covariates=None)`.
- Extend `prepare_complete_data_projection(expression, principal_components, requested_pc_counts, additional_covariates=None)`.
- Keep `CompleteDataProjection.advance_prediction(previous_pc_count, pc_count, prediction)` and `z_scores(prediction)` usable by existing callers, and add `initial_prediction() -> np.ndarray` that returns the fixed additional-covariate fitted component before the molecular-PC increments.

- [ ] **Step 1: Write the failing model-equivalence tests**

Add a helper in `tests/test_lof_pc.py` that computes a reference fit with `np.linalg.lstsq` on `np.column_stack((np.ones(n), additional_covariates, principal_components[:, :pc_count]))`, then centers residuals and scales by population SD. Add tests that:

1. Compare `residualize_expression` against this reference for `pc_count=0` and `pc_count=1`.
2. Verify `pc_count=0` still removes the additional covariate effect.
3. Compare the vectorized projection output against `residualize_expression` for counts `[0, 1, 2]` with non-orthogonal molecular PCs and two additional covariates.
4. Verify a collinear additional covariate causes the existing rank-deficiency/fallback behavior rather than silently producing a different model.

- [ ] **Step 2: Run the model tests and verify RED**

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py -k 'residualize or complete_data_projection'
```

Expected: the new additional-covariate equivalence tests fail because the functions do not accept or use the new matrix.

- [ ] **Step 3: Implement the legacy residualization extension**

Update `residualize_expression` to validate the optional matrix shape against the expression vector, require finite covariates for usable rows, and build the design in the exact order `intercept`, all additional covariates, then the selected molecular-PC prefix. Keep rank, insufficient-DOF, residual-SD, and missing-expression handling unchanged.

- [ ] **Step 4: Implement the vectorized fixed-plus-prefix projection**

Update `prepare_complete_data_projection` to accept the optional matrix and construct a centered, orthogonal design with a fixed covariate block followed by a molecular-PC prefix block. Store the fixed covariate prediction separately from the incremental molecular-PC prediction. At each requested count, the prediction must equal the least-squares fitted values from the reference design, including the fixed block when `pc_count=0`. Detect invalid or rank-deficient fixed/PC blocks and let the existing legacy fallback execute.

- [ ] **Step 5: Run the model tests and verify GREEN**

Run the focused model command from Step 2, then the existing projection regression tests:

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py -k 'complete_data_projection or residualize_expression or vectorized_multi_pc'
```

Expected: both new covariate tests and legacy no-covariate tests pass.

- [ ] **Step 6: Commit the model extension**

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: adjust LoF PC models for additional covariates"
```

### Task 3: Integrate sample intersection, QC, and provenance into enrichment

**Files:**
- Modify: `tests/test_lof_pc.py`
- Modify: `src/rare_variant_enrichment/lof_pc.py`

**Interfaces:**
- Extend `calculate_lof_pc_enrichment(..., *, pc_grid_mode=None, additional_covariates_path=None)`.
- The final analysis sample order remains the BED order restricted to the molecular-PC and optional covariate sample sets.

- [ ] **Step 1: Write the failing end-to-end intersection and metadata tests**

Extend `_write_analysis_fixture` with an additional-covariate table whose sample order differs from the BED and PC files and whose sample set omits one PC/BED sample. Call `_run_analysis` with `additional_covariates_path` and assert:

```python
analysis_qc["additional_covariates_supplied"] is True
analysis_qc["additional_covariate_count"] == 2
analysis_qc["additional_covariate_names"] == ["GENETICPC1", "GENETICPC2"]
analysis_qc["additional_covariate_sample_count"] == 5
analysis_qc["shared_bed_pc_sample_count"] == 5
analysis_qc["shared_bed_pc_covariate_sample_count"] == 5
```

Assert that the summary residualization design describes additional covariates and that `summary["provenance"]["input_files"]["additional_covariates"]` records the supplied path. Add a no-covariate assertion that the legacy QC count and design remain unchanged. Add a no-shared-sample test that expects `ValueError` containing `shared`.

- [ ] **Step 2: Run the focused integration tests and verify RED**

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py -k 'enrichment or intersection or provenance'
```

Expected: the new tests fail because the calculation function has no covariate input or metadata wiring.

- [ ] **Step 3: Implement optional matrix loading and sample alignment**

In `calculate_lof_pc_enrichment`, read the optional matrix once. Intersect its sample IDs with the BED/molecular-PC IDs while preserving BED order, reorder all three matrices to that order, and fail if the intersection is empty. When absent, retain the existing two-input intersection and set covariate metadata to an explicit false/zero/empty representation.

- [ ] **Step 4: Pass aligned covariates through both analysis paths**

Pass the aligned covariate matrix to `prepare_complete_data_projection` on the complete-data path and to `residualize_expression` on the fallback path. Include the fixed covariate block in rank calculations and in the rank-deficiency threshold (`number of covariates + pc_count + 1` columns including the intercept).

- [ ] **Step 5: Emit and validate QC/provenance metadata**

Add top-level analysis QC keys for covariate presence, count, names, input sample count, and final shared sample count while retaining the existing `shared_bed_pc_sample_count`. Add the optional path under summary provenance input files when supplied. Set the residualization design string to the covariate-aware form when supplied and retain the legacy string otherwise. Update merge metadata validation/tests so shards with different covariate metadata are rejected and identical covariate metadata merge successfully.

- [ ] **Step 6: Run integration and merge tests and verify GREEN**

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_lof_pc.py
```

Expected: all LoF/PC tests, including legacy and additional-covariate paths, pass.

- [ ] **Step 7: Commit the enrichment integration**

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: integrate covariate sample alignment"
```

### Task 4: Expose the optional covariate file through the CLI and WDL

**Files:**
- Modify: `src/rare_variant_enrichment/cli.py`
- Modify: `workflows/rare_variant_enrichment.wdl`
- Modify: `tests/test_cli.py`
- Modify: `tests/test_wdl_contract.py`
- Modify: `tests/test_wdl_runtime.py`
- Modify: `examples/rare_variant_enrichment.inputs.json`

**Interfaces:**
- CLI option: `--additional-covariates PATH`, optional for `lof-pc-enrichment`.
- WDL workflow input: `File? additional_covariates_tsv` with no required default.
- WDL workflow input: `Int pc_preemptible = 2`.
- `CalculateLofPcEnrichment` task input: `File? additional_covariates_tsv` and `Int preemptible`.

- [ ] **Step 1: Write failing CLI and WDL contract tests**

Update the CLI dispatch test fake to accept `**keyword_arguments`, pass `--additional-covariates genetic-pcs.tsv`, and assert that the calculation call receives `additional_covariates_path=Path("genetic-pcs.tsv")` while the no-option test remains valid. Extend the parsed WDL contract to assert:

```python
contract["inputs"]["additional_covariates_tsv"] == {"type": "File?", "default": None}
contract["inputs"]["pc_preemptible"] == {"type": "Int", "default": 2}
contract["task_inputs"]["CalculateLofPcEnrichment"]["preemptible"] == "Int"
```

Assert that the scatter passes `additional_covariates_tsv` and `preemptible: pc_preemptible`, that the task runtime includes `preemptible`, and that the calculated analysis disk expression includes the optional covariate size through a zero-when-undefined declaration. Update the command-rendering fixture with a safe optional covariate file path and assert it is quoted when supplied.

- [ ] **Step 2: Run the focused interface tests and verify RED**

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_cli.py tests/test_wdl_contract.py
```

Expected: the new assertions fail because the CLI and WDL do not expose the optional file or preemptible runtime value.

- [ ] **Step 3: Implement the CLI option**

Add `--additional-covariates` as an optional `Path` argument and pass it as the keyword-only `additional_covariates_path` argument to `calculate_lof_pc_enrichment`. Preserve the existing positional argument order so current callers remain valid.

- [ ] **Step 4: Implement WDL optional-file wiring and dynamic disk sizing**

Add `File? additional_covariates_tsv` to the workflow and `CalculateLofPcEnrichment` task. Render an empty command fragment when undefined and a shell-quoted `--additional-covariates` argument when defined. Add `Float additional_covariates_size_gb = if defined(additional_covariates_tsv) then size(additional_covariates_tsv, "GiB") else 0.0`, include it in the analysis disk calculation, and pass the optional file to every scatter shard.

- [ ] **Step 5: Implement configurable preemptible attempts**

Add `Int preemptible` to `CalculateLofPcEnrichment`, set `preemptible: preemptible` in its runtime block, add workflow input `Int pc_preemptible = 2`, and pass `pc_preemptible` to each scatter shard. Do not change `maxRetries` or add preemptible settings to preparation, merge, or plotting tasks.

- [ ] **Step 6: Update example inputs and runtime fixture**

Add `RareVariantEnrichment.pc_preemptible: 2` to `examples/rare_variant_enrichment.inputs.json`. Add a small optional genotype-PC fixture to `tests/fixtures` or generate it in `tests/test_wdl_runtime.py`, pass it to the runtime test, and assert that the output QC reports the reduced shared sample count and covariate count. Keep the existing no-covariate runtime fixture to prove backward compatibility.

- [ ] **Step 7: Run interface tests and verify GREEN**

Run:

```bash
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q tests/test_cli.py tests/test_wdl_contract.py
PYTHONPATH="/private/tmp/rve-wdl-py311" /opt/homebrew/bin/miniwdl check workflows/rare_variant_enrichment.wdl
PYTHONPATH="/private/tmp/rve-wdl-py311" /opt/homebrew/bin/miniwdl input_template workflows/rare_variant_enrichment.wdl
```

Expected: CLI and contract tests pass, WDL syntax checks pass, and the input template shows the optional covariate file as optional while retaining the required four existing files.

- [ ] **Step 8: Commit the public-interface changes**

```bash
git add src/rare_variant_enrichment/cli.py workflows/rare_variant_enrichment.wdl tests/test_cli.py tests/test_wdl_contract.py tests/test_wdl_runtime.py tests/fixtures examples/rare_variant_enrichment.inputs.json
git commit -m "feat: expose genotype covariates and PC preemption"
```

### Task 5: Document the analysis model and perform complete verification

**Files:**
- Modify: `README.md`
- [ ] **Step 1: Update README usage and statistical model text**

Document the optional `additional_covariates_tsv` input, accepted `sample_id`/`ID` column placement, string-preserving IDs, automatic inclusion of all numeric columns, the BED/PC/covariate intersection, and the fact that LoF carrier rows do not define the sample universe. Update the residualization equation to include all additional covariates before the molecular-PC prefix. Document `pc_preemptible = 2` and that it applies only to scattered PC-fitting tasks.

- [ ] **Step 2: Run formatting and complete tests**

Run:

```bash
git diff --check
PYTHONPATH="/private/tmp/rve-pytest-deps:/private/tmp/rve-numpy-deps:/private/tmp/rve-wdl-py311:src" \
  /opt/homebrew/bin/python3 -m pytest -q
PYTHONPATH="/private/tmp/rve-wdl-py311" /opt/homebrew/bin/miniwdl check workflows/rare_variant_enrichment.wdl
```

Expected: the full suite passes with the established baseline of `232 passed, 4 skipped` plus the new tests, WDL validation passes, and `git diff --check` is clean. If Docker prerequisites are available, also run the covariate-enabled runtime test; otherwise report it as skipped rather than claiming runtime validation.

- [ ] **Step 3: Review generated metadata and diff**

Inspect one covariate-enabled analysis QC JSON, summary JSON, results TSV, and gene-PC QC file. Confirm the result and gene-PC schemas are unchanged, the shared sample counts are correct, the residualization design names the covariates, and all changed files are limited to this feature.

- [ ] **Step 4: Commit documentation and final verification**

```bash
git add README.md
git commit -m "docs: describe genotype covariate adjustment"
git status --short --branch
git log --oneline -6
```
