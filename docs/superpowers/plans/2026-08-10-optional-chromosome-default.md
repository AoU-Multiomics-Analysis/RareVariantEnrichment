# Optional Chromosome Default Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the workflow chromosome array optional for callers by defaulting it to `chr1` through `chr22`, while preserving explicit single- and multi-chromosome overrides.

**Architecture:** Keep `chromosomes` as a concrete `Array[String]` and add a WDL default, so every existing call, scatter, and provenance path continues to consume the same resolved value. Update the public contract test, example JSON, and README without adding nullable-input handling or changing task interfaces.

**Tech Stack:** WDL 1.0, miniwdl, Python 3.12+ with pytest, JSON, Markdown.

## Global Constraints

- The default chromosome array is exactly `chr1` through `chr22` in numeric order.
- The default excludes `chrX`, `chrY`, mitochondrial contigs, and alternate contigs.
- Explicit chromosome arrays retain their existing behavior and must exactly match BED and VCF contig labels.
- The resolved array continues to control VCF validation, phenotype selection, scatter shards, calculation denominators, and provenance.
- Do not change Python package interfaces, statistical behavior, tabix behavior, or task runtime settings.

---

### Task 1: Default the chromosome array and document overrides

**Files:**
- Modify: `tests/test_wdl_contract.py`
- Modify: `workflows/rare_variant_enrichment.wdl`
- Modify: `examples/rare_variant_enrichment.inputs.json`
- Modify: `README.md`

**Interfaces:**
- Consumes: existing workflow input `Array[String] chromosomes`.
- Produces: `Array[String] chromosomes` with the WDL default `["chr1", ..., "chr22"]`; callers may omit it or override it with any non-empty explicit array accepted by existing validation.

- [ ] **Step 1: Write failing WDL contract and example tests**

Update `test_wdl_public_input_and_default_contract` so the miniwdl input
template contains only the two required file inputs:

```python
assert json.loads(template_result.stdout) == {
    "RareVariantEnrichment.phenotype_bed": "File",
    "RareVariantEnrichment.rare_variant_vcf": "File",
}
```

Define the expected default once in the test:

```python
expected_autosomes = [f"chr{chromosome}" for chromosome in range(1, 23)]
```

Change the parsed input assertion for `chromosomes` to:

```python
"chromosomes": {
    "type": "Array[String]",
    "default": expected_autosomes,
},
```

Add a focused example-input assertion:

```python
def test_example_uses_default_autosomes():
    inputs = json.loads(
        (ROOT / "examples" / "rare_variant_enrichment.inputs.json").read_text()
    )
    assert "RareVariantEnrichment.chromosomes" not in inputs
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run:

```bash
.venv/bin/python -m pytest +  tests/test_wdl_contract.py::test_wdl_public_input_and_default_contract +  tests/test_wdl_contract.py::test_example_uses_default_autosomes -v
```

Expected: both assertions fail because `chromosomes` is still required and
the example still supplies it.

- [ ] **Step 3: Add the minimal WDL default**

Replace the required workflow declaration with:

```wdl
Array[String] chromosomes = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
    "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22"
]
```

Do not change any downstream references: every call, scatter, and provenance
input continues to use `chromosomes`.

- [ ] **Step 4: Update the example and README**

Remove the `RareVariantEnrichment.chromosomes` property from
`examples/rare_variant_enrichment.inputs.json`.

In the README input section, state:

```markdown
`chromosomes` defaults to all autosomes, `chr1` through `chr22`. Omit
the input for that default. To restrict a run, supply an explicit array such as
`["chr7"]` or `["chr1", "chr2"]`. Non-autosomal contigs such as `chrX`
must be requested explicitly, and every label must match both BED and VCF
contigs.
```

Update the run section to explain that the example omits `chromosomes` to
demonstrate the autosomal default.

- [ ] **Step 5: Verify GREEN and unchanged explicit override behavior**

Run:

```bash
.venv/bin/python -m pytest tests/test_wdl_contract.py tests/test_wdl_runtime.py -v -rs
```

Expected: all contract tests pass; existing Docker-backed chr1 runtime cases
pass when prerequisites are available and continue to report
`selected_chromosomes == ["chr1"]`.

- [ ] **Step 6: Run complete verification**

Run:

```bash
.venv/bin/python -m pytest -v -rs
miniwdl check workflows/rare_variant_enrichment.wdl
.venv/bin/python -m json.tool examples/rare_variant_enrichment.inputs.json
git diff --check
```

Expected: the full suite has zero failures, WDL validation exits 0, the example
is valid JSON, and the diff has no whitespace errors.

- [ ] **Step 7: Review scope and commit**

Confirm the diff changes only the chromosome default, its public contract,
documentation, and example. Then run:

```bash
git add \
  workflows/rare_variant_enrichment.wdl \
  tests/test_wdl_contract.py \
  examples/rare_variant_enrichment.inputs.json \
  README.md
git commit -m "feat: default workflow to autosomes"
```
