# Enrichment Progress Logging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every LoF/PC enrichment run emit concise, actionable progress messages to standard error.

**Architecture:** Keep logging inside `calculate_lof_pc_enrichment`, the orchestration boundary that knows the validated inputs, sample overlap, processing loop, aggregate PC statistics, and final output paths. Use Python's standard-library `logging` module with a module logger so the CLI can configure a predictable stderr handler without changing result files.

**Tech Stack:** Python 3, standard-library `logging`, NumPy, pytest.

## Global Constraints

- Logging is always enabled; no verbosity flag is introduced.
- Messages go to standard error and must not change TSV, gzip, or JSON outputs.
- Log no individual sample or gene IDs during the processing loop.
- Log periodic BED processing every 500 protein-coding genes.

---

### Task 1: Add always-on LoF/PC enrichment progress logging

**Files:**
- Modify: `src/rare_variant_enrichment/cli.py:1-140`
- Modify: `src/rare_variant_enrichment/lof_pc.py:1-660`
- Modify: `tests/test_lof_pc.py:270-360`

**Interfaces:**
- Consumes: `calculate_lof_pc_enrichment(...) -> None` and its existing file-path arguments.
- Produces: stderr messages for configuration, BED/PC sample overlap, every 500 coding BED genes, PC completion totals, and output completion.

- [ ] **Step 1: Write the failing stderr regression test**

```python
def test_analysis_logs_configuration_progress_and_outputs(tmp_path: Path, caplog):
    inputs = _write_analysis_fixture(tmp_path)

    with caplog.at_level(logging.INFO, logger="rare_variant_enrichment.lof_pc"):
        outputs = _run_analysis(tmp_path, inputs, thresholds=[-0.8], pc_counts=[0])

    messages = "\\n".join(caplog.messages)
    assert "Starting LoF/PC enrichment" in messages
    assert "shared BED/PC samples" in messages
    assert "Completed PC count 0" in messages
    assert "Wrote LoF/PC enrichment outputs" in messages
    assert str(outputs["results"]) in messages
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest tests/test_lof_pc.py::test_analysis_logs_configuration_progress_and_outputs -v`

Expected: FAIL because `calculate_lof_pc_enrichment` emits no progress messages.

- [ ] **Step 3: Add a standard-error logging configuration to the CLI**

```python
import logging

def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = build_parser().parse_args()
```

Use the existing `main()` entrypoint so WDL command invocations always configure a stderr handler. Do not add a logging argument.

- [ ] **Step 4: Add bounded messages at analysis milestones**

```python
LOGGER = logging.getLogger(__name__)
PROGRESS_INTERVAL_GENES = 500

LOGGER.info("Starting LoF/PC enrichment: thresholds=%s pc_counts=%s coding_genes=%d", thresholds, pc_counts, len(coding_genes))
LOGGER.info("Found %d shared BED/PC samples (%d BED, %d PC)", len(shared_samples), len(bed_samples), len(principal_components.sample_ids))
if coding_bed_gene_count % PROGRESS_INTERVAL_GENES == 0:
    LOGGER.info("Processed %d protein-coding BED genes", coding_bed_gene_count)
for pc_count in pc_counts:
    LOGGER.info("Completed PC count %d: eligible_genes=%d observations=%d", pc_count, per_pc[str(pc_count)]["eligible_gene_count"], per_pc[str(pc_count)]["total_observations"])
LOGGER.info("Wrote LoF/PC enrichment outputs: results=%s summary=%s gene_pc_qc=%s analysis_qc=%s", ...)
```

Place the configuration log after validating inputs and reading reference files. Place overlap logging immediately after the no-shared-samples guard. Increment the coding-gene counter before the periodic log, and log PC completion only after the entire BED loop. Emit the final message after writing all outputs.

- [ ] **Step 5: Run the focused regression test**

Run: `pytest tests/test_lof_pc.py::test_analysis_logs_configuration_progress_and_outputs -v`

Expected: PASS and messages include configuration, sample overlap, PC completion, and output paths.

- [ ] **Step 6: Run the complete validation suite**

Run: `git diff --check && pytest -q && miniwdl check workflows/rare_variant_enrichment.wdl`

Expected: zero whitespace errors, all Python tests passing, and WDL validation succeeding.

- [ ] **Step 7: Commit the implementation**

```bash
git add src/rare_variant_enrichment/cli.py src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "feat: log LoF PC enrichment progress"
```
