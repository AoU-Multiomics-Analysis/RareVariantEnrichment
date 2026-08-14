# Vectorized LoF/PC Residualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Use complete-data NumPy projections instead of repeated least-squares fits for the LoF/PC enrichment sweep, while preserving output semantics and emitting live PC-level progress.

**Architecture:** Add a complete-data projection helper that centers and normalizes the PCA score matrix, then residualizes a sample-by-gene expression block using matrix multiplication. `calculate_lof_pc_enrichment` selects this helper when every coding BED expression value is finite, aggregates one requested PC count at a time, and retains the current per-gene route for missing-expression inputs.

**Tech Stack:** Python 3, NumPy, pytest, WDL.

## Global Constraints

- Existing explicit and adaptive `pc_counts` selection remains unchanged; no PC-count cap.
- Fast path is used only when all shared-sample coding expression values are finite.
- Missing-expression inputs retain the current per-gene regression behavior.
- Fast-path residuals use an intercept, centered and normalized PC columns, population SD (`ddof=0`), and existing exclusion semantics.
- TSV, gzip QC, and JSON output schemas remain unchanged.
- Process and log requested PC counts one at a time; never log individual sample or gene IDs in the inner loop.

---

### Task 1: Implement complete-data projection primitives

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:341-410`
- Modify: `tests/test_lof_pc.py:204-246`

**Interfaces:**
- Consumes: an expression matrix `(samples, genes)`, PC matrix `(samples, available_pcs)`, and requested PC counts.
- Produces: `CompleteDataProjection` with `advance_prediction(previous_pc_count: int, pc_count: int, prediction: np.ndarray) -> np.ndarray` and `z_scores(prediction: np.ndarray) -> np.ndarray`.

- [ ] **Step 1: Write failing equivalence tests**

```python
def test_complete_data_projection_matches_legacy_residuals():
    state = prepare_complete_data_projection(expression, pcs, [0, 1])
    prediction = np.zeros_like(expression)
    previous_pc_count = 0
    for pc_count in [0, 1]:
        prediction = state.advance_prediction(previous_pc_count, pc_count, prediction)
        actual = state.z_scores(prediction)
        expected = np.column_stack([
            residualize_expression(expression[:, column], pcs, pc_count).z_scores
            for column in range(expression.shape[1])
        ])
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-10)
        previous_pc_count = pc_count
```

Also test zero-variance expression exclusion and rejection of any missing
expression value by the fast-path helper.

- [ ] **Step 2: Verify RED**

Run `/private/tmp/rve-rebuild-venv/bin/pytest tests/test_lof_pc.py -k complete_data_projection -v`.

Expected: FAIL because the projection API does not exist.

- [ ] **Step 3: Implement the helper**

```python
@dataclass(frozen=True)
class CompleteDataProjection:
    centered_expression: np.ndarray
    normalized_pcs: np.ndarray
    component_scores: np.ndarray

    def advance_prediction(self, previous_pc_count, pc_count, prediction):
        return prediction + self.normalized_pcs[:, previous_pc_count:pc_count] @ self.component_scores[previous_pc_count:pc_count]

    def z_scores(self, prediction: np.ndarray) -> np.ndarray:
        residuals = self.centered_expression - prediction
        return residuals / np.std(residuals, axis=0, ddof=0)
```

Center expression and PC columns, normalize nonzero PC columns, validate PC
counts, and return NaN columns for zero/non-finite residual SD so the caller
uses existing exclusion categories.

- [ ] **Step 4: Verify GREEN and commit**

Run `/private/tmp/rve-rebuild-venv/bin/pytest tests/test_lof_pc.py -k complete_data_projection -v` and then:

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: add vectorized complete-data PC projection"
```

### Task 2: Integrate PC-major aggregation and logging

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:413-684`
- Modify: `tests/test_lof_pc.py:280-470`

**Interfaces:**
- Consumes: `CompleteDataProjection` from Task 1 plus existing carrier, gene, threshold, and output inputs.
- Produces: unchanged enrichment TSV/QC/JSON files and PC-completion logs emitted during the sweep.

- [ ] **Step 1: Write failing integration tests**

```python
def test_complete_data_analysis_logs_pc_counts_in_processing_order(tmp_path, caplog):
    with caplog.at_level(logging.INFO, logger="rare_variant_enrichment.lof_pc"):
        outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0, 1])
    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert {row["pc_count"] for row in rows} == {"0", "1"}
    completion_indexes = {
        pc_count: next(index for index, message in enumerate(caplog.messages)
                       if f"Completed PC count {pc_count}" in message)
        for pc_count in (0, 1)
    }
    assert completion_indexes[0] < completion_indexes[1]
```

Add a missing-expression fixture asserting the legacy route's existing
results and QC exclusions are retained.

- [ ] **Step 2: Verify RED**

Run `/private/tmp/rve-rebuild-venv/bin/pytest tests/test_lof_pc.py -k 'pc_major or missing_expression_fallback' -v`.

Expected: FAIL because the current loop is gene-major and does not choose the
projection path.

- [ ] **Step 3: Integrate PC-major processing**

Read all coding BED rows for the complete-data path, construct the projection
state, and initialize one zero-valued prediction matrix. For each requested
PC count, add only the PC interval since the preceding requested count with
`advance_prediction`, standardize its residual matrix with `z_scores`, then
aggregate carrier masks, outlier cells, and gene-PC QC rows with existing
headers and definitions. Emit `Completed PC count` immediately after each
count. Retain the current line-by-line implementation as the path selected if
any coding expression entry is missing.

- [ ] **Step 4: Verify GREEN**

Run `/private/tmp/rve-rebuild-venv/bin/pytest tests/test_lof_pc.py -k 'pc_major or missing_expression_fallback or enrichment_emits_hand_checked' -v`.

Expected: PASS with unchanged contingency cells and ordered PC logs.

- [ ] **Step 5: Run full validation and commit**

Run `git diff --check && /private/tmp/rve-rebuild-venv/bin/pytest -q && /private/tmp/rve-rebuild-venv/bin/miniwdl check workflows/rare_variant_enrichment.wdl`.

Expected: no whitespace errors, all Python tests passing, and WDL validation succeeds. Then run:

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: vectorize LoF PC enrichment sweep"
```
