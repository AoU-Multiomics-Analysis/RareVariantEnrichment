# Scattered LoF/PC Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Scatter the LoF/PC enrichment sweep over groups of selected PC-count settings and merge the shard outputs without changing the six public workflow outputs or global FDR semantics.

**Architecture:** Add a header-only PC-grid chunking CLI operation and a `PreparePcChunks` WDL task. Scatter the existing analysis task over `Array[Int]` PC-count chunks, then call a merge task that validates shard metadata, combines QC, concatenates rows, and recomputes global Benjamini–Hochberg FDR.

**Tech Stack:** Python 3, NumPy, pytest, WDL 1.0, MiniWDL.

## Global Constraints

- `pc_counts_per_job` is the number of selected PC-count settings per shard, not the number of PCA columns.
- `pc_counts_per_job` defaults to `10` and must be a positive integer.
- Both explicit `pc_counts` and the existing adaptive grid for `pc_counts=[]` remain supported.
- PC chunk preparation reads only the PC header and must not parse the full PC matrix.
- Each analysis shard receives at least one PC-count setting and runs one complete analysis call.
- Merge recomputes Fisher FDR globally across all merged result rows.
- Results TSV, summary JSON, gene-PC QC gzip, analysis-QC JSON, and the two coding-gene outputs retain their current public schemas and names.
- Existing missing-expression, rank-deficiency, numerical-exclusion, and residualization behavior remains unchanged.
- The pre-existing untracked `src/rare_variant_enrichment.egg-info/` directory is not staged.

---

### Task 1: Add header-only PC grid chunking

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py` near `read_principal_components` and `build_pc_grid`
- Modify: `src/rare_variant_enrichment/cli.py` parser and dispatch
- Modify: `tests/test_lof_pc.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Add `read_principal_component_header(path: Path) -> int`, which validates the first row is `ID`, `PC1`, …, consecutive and returns the available PC count without reading data rows.
- Add `build_pc_chunks(requested_pc_counts: Sequence[int], available_pc_count: int, pc_counts_per_job: int) -> list[list[int]]`, which calls the existing `build_pc_grid`, validates a positive chunk size, and partitions the selected grid in order.
- Add CLI command `pc-chunks` with `--principal-components`, `--pc-counts`, `--pc-counts-per-job`, and `--output`, writing a JSON array of integer arrays.

- [ ] **Step 1: Write failing tests for explicit and adaptive chunks**

```python
def test_build_pc_chunks_partitions_explicit_grid_with_short_final_chunk():
    assert build_pc_chunks([0, 1, 10, 20, 30], 30, 2) == [[0, 1], [10, 20], [30]]


def test_build_pc_chunks_uses_adaptive_grid_when_requested_grid_is_empty():
    assert build_pc_chunks([], 25, 3) == [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 20], [25]]


def test_build_pc_chunks_rejects_nonpositive_chunk_size():
    with pytest.raises(ValueError, match="positive"):
        build_pc_chunks([0], 1, 0)
```

Also add a header-only fixture test that places invalid numeric data after a valid PC header and verifies `read_principal_component_header` still returns the header count. Add a CLI dispatch test asserting the JSON output path and parsed PC-count/chunk-size arguments.

- [ ] **Step 2: Run tests to verify the new APIs fail**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k 'pc_chunks or principal_component_header' tests/test_cli.py -k pc_chunks
```

Expected: FAIL because the header reader, chunk builder, and CLI command do not yet exist.

- [ ] **Step 3: Implement header validation, chunking, and CLI dispatch**

Use the existing `build_pc_grid` as the sole source of adaptive-grid behavior:

```python
def build_pc_chunks(
    requested_pc_counts: Sequence[int],
    available_pc_count: int,
    pc_counts_per_job: int,
) -> list[list[int]]:
    if (
        isinstance(pc_counts_per_job, bool)
        or not isinstance(pc_counts_per_job, int)
        or pc_counts_per_job <= 0
    ):
        raise ValueError("pc_counts_per_job must be a positive integer")
    selected = build_pc_grid(requested_pc_counts, available_pc_count)
    return [
        selected[start : start + pc_counts_per_job]
        for start in range(0, len(selected), pc_counts_per_job)
    ]
```

The CLI must parse `--pc-counts` with the existing CSV integer parser, allow an empty value for adaptive mode, call the header-only reader, and write the chunk list with `write_json`.

- [ ] **Step 4: Run focused tests and commit**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k 'pc_chunks or principal_component_header' tests/test_cli.py -k pc_chunks
git diff --check
git add src/rare_variant_enrichment/lof_pc.py src/rare_variant_enrichment/cli.py tests/test_lof_pc.py tests/test_cli.py
git commit -m "feat: add PC-count chunking"
```

### Task 2: Add merged shard-output operation

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py` after result-writing helpers
- Modify: `src/rare_variant_enrichment/cli.py`
- Modify: `tests/test_lof_pc.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Add `merge_lof_pc_enrichment(results_inputs: Sequence[Path], summary_inputs: Sequence[Path], gene_pc_qc_inputs: Sequence[Path], analysis_qc_inputs: Sequence[Path], results_output: Path, summary_output: Path, gene_pc_qc_output: Path, analysis_qc_output: Path) -> None`.
- Add CLI command `merge-lof-pc-enrichment` taking four newline-delimited input-list files (`--results-input-list`, `--summary-input-list`, `--gene-pc-qc-input-list`, `--analysis-qc-input-list`) and four output paths matching the existing analysis outputs.

- [ ] **Step 1: Write failing merge tests**

Create two synthetic shard outputs with disjoint PC counts and known result rows. Assert that the merged TSV has one header, deterministic PC/threshold/carrier ordering, and FDR values recomputed over all rows rather than copied from shard-local values. Assert that gene-PC QC has one gzip header and all shard rows, and that analysis-QC `per_pc` counters are present for every PC while shared top-level metadata are preserved.

```python
def test_merge_lof_pc_enrichment_recomputes_global_fdr_and_combines_qc(tmp_path):
    shard_one = _write_merge_shard(tmp_path, "one", pc_counts=[0], p_values=[0.01, 0.20])
    shard_two = _write_merge_shard(tmp_path, "two", pc_counts=[1], p_values=[0.02, 0.50])
    outputs = _merge_shards(tmp_path, [shard_one, shard_two])
    rows = list(csv.DictReader(outputs["results"].open(), delimiter="\t"))
    assert [row["pc_count"] for row in rows] == ["0", "0", "1", "1"]
    assert [row["fisher_fdr_bh"] for row in rows] == [
        "0.04", "0.26666666666666666", "0.04", "0.5"
    ]
    with gzip.open(outputs["gene_qc"], "rt", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 2
    assert set(json.loads(outputs["analysis_qc"].read_text())["per_pc"]) == {"0", "1"}
```

The `_write_merge_shard` test helper must write a valid result header, one
row per supplied p-value, one gene-PC QC row per PC, and a summary/analysis-QC
JSON matching the existing schemas. `_merge_shards` must invoke the new merge
function with the four corresponding shard path lists and return the four
output paths.

- [ ] **Step 2: Run the merge test to verify it fails**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k merge_lof_pc tests/test_cli.py -k merge_lof_pc
```

Expected: FAIL because the merge operation and CLI command do not yet exist.

- [ ] **Step 3: Implement validated merge and global FDR**

Read all four input lists, require equal nonzero shard counts, and reject duplicate PC counts. Validate each summary’s thresholds, carrier definitions, available PC count, provenance, and selected PC-count set against the first shard. Parse result rows using `RESULT_HEADER`, concatenate and sort by selected PC order followed by each shard’s existing threshold/carrier order, collect all `fisher_p_value` values, and overwrite `fisher_fdr_bh` using `benjamini_hochberg` across the complete row set.

Concatenate gene-PC QC files after verifying identical headers. For analysis-QC JSON, retain the first shard’s shared top-level metadata and sum each per-PC `eligible_gene_count`, `total_observations`, `carrier_observations`, and `exclusion_counts`. Set the merged summary’s `selected_pc_counts` and `emitted_result_rows` to the union and final row count. Write the four outputs with the existing schemas.

- [ ] **Step 4: Run focused merge tests and commit**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k merge_lof_pc tests/test_cli.py -k merge_lof_pc
git diff --check
git add src/rare_variant_enrichment/lof_pc.py src/rare_variant_enrichment/cli.py tests/test_lof_pc.py tests/test_cli.py
git commit -m "feat: merge scattered LoF PC outputs"
```

### Task 3: Wire WDL scatter, merge, resources, and documentation

**Files:**
- Modify: `workflows/rare_variant_enrichment.wdl`
- Modify: `examples/rare_variant_enrichment.inputs.json`
- Modify: `README.md`
- Modify: `tests/test_wdl_contract.py`
- Modify: `tests/test_wdl_runtime.py`

**Interfaces:**
- Add workflow input `Int pc_counts_per_job = 10`.
- Add task `PreparePcChunks` that consumes the PC file, `pc_counts`, `pc_counts_per_job`, and container/resource inputs, and outputs `File pc_chunks_json` containing `Array[Array[Int]]` JSON.
- Scatter `CalculateLofPcEnrichment` over the parsed chunk array, passing each `Array[Int]` chunk as `pc_counts`.
- Add task `MergeLofPcEnrichment` that consumes arrays of the four shard outputs and emits the four merged analysis outputs.
- Keep `RareVariantEnrichment`’s six public outputs unchanged.

- [ ] **Step 1: Add failing WDL contract tests**

Update contract expectations to require `pc_counts_per_job`, `PreparePcChunks`, `CalculateLofPcEnrichment`, and `MergeLofPcEnrichment`; require the scatter and merge wiring; require `pc_counts_per_job = 10`; and assert that the final output names/types remain unchanged. Add a runtime input with `pc_counts=[0,1]` and `pc_counts_per_job=1` so the fixture executes two analysis shards and one merge.

- [ ] **Step 2: Run WDL tests to verify the new contract fails**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_wdl_contract.py tests/test_wdl_runtime.py -k 'pc_chunk or scatter or merge'
```

Expected: FAIL because the current WDL has no chunk task, scatter, merge task, or chunk-size input.

- [ ] **Step 3: Implement the WDL tasks and workflow wiring**

`PreparePcChunks` will call `rare-variant-enrichment pc-chunks`, materializing an explicit empty PC-count list when adaptive mode is requested. The workflow will parse the JSON with `read_json`, scatter over the resulting `Array[Array[Int]]`, and pass the scatter arrays to `MergeLofPcEnrichment`.

The merge command will receive newline-delimited shard file lists generated with `write_lines`. Quote every localized filename in shell commands. Keep current analysis resources for shards; size chunk-preparation disk from the PC file with a `ceil(2 × GiB + 20)` floor and size merge disk from the sum of localized shard outputs with the same dynamic floor, while retaining the configured analysis disk as the baseline. Propagate `max_retries` to every new task.

- [ ] **Step 4: Update examples and README**

Document `pc_counts_per_job`, its meaning as settings per job, the default of 10, and the fact that final FDR is recomputed after merge. Add the input to the example JSON and describe that intermediate scatter files are not public outputs.

- [ ] **Step 5: Run runtime and full validation**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_wdl_contract.py tests/test_wdl_runtime.py
/private/tmp/rve-rebuild-venv/bin/pytest -q
/private/tmp/rve-rebuild-venv/bin/miniwdl check workflows/rare_variant_enrichment.wdl
git diff --check
```

Expected: all tests pass, the two-shard fixture produces the same known Fisher cells and global FDR scope, and MiniWDL validates the scatter/merge workflow.

- [ ] **Step 6: Commit the WDL integration**

```bash
git add workflows/rare_variant_enrichment.wdl examples/rare_variant_enrichment.inputs.json README.md tests/test_wdl_contract.py tests/test_wdl_runtime.py
git commit -m "feat: scatter and merge LoF PC enrichment"
```
