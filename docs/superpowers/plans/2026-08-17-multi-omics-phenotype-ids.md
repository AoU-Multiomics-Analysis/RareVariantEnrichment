# Multi-Omics Phenotype IDs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the LoF/PC enrichment pipeline accept expression, proteomics, and splicing phenotype IDs by normalizing them to unversioned Ensembl gene IDs and collapsing duplicate splice features to the most extreme negative z-score per gene and sample.

**Architecture:** Keep the WDL public interface and downstream gene-level Fisher testing unchanged. Add phenotype-specific ID parsing and a small gene-row aggregation helper in `lof_pc.py`; the existing residualization code will consume one canonical gene vector per protein-coding gene. Extend analysis QC and README documentation to make feature-to-gene collapsing explicit.

**Tech Stack:** Python 3.12+, NumPy, pytest, miniwdl, existing WDL workflow.

## Global Constraints

- Molecular phenotype IDs must map to exactly one Ensembl gene token of the form `ENSG` followed by digits, with an optional numeric version suffix.
- Expression, proteomics, and splicing phenotype rows must use the same sample-column and z-score parsing rules already used by `lof-pc-enrichment`.
- Duplicate phenotype rows for one canonical gene must be collapsed sample-wise with the minimum finite z-score; missing values are ignored.
- The existing carrier matching, residualization, negative-z thresholding, Fisher tests, output filenames, and WDL public inputs remain unchanged.
- The existing `normalize_ensembl_id` behavior for GTF and carrier-table IDs must not change.
- Every implementation change must have a regression test written and observed failing before production code is added.

---

### Task 1: Add failing tests for multi-omics IDs and gene-level collapsing

**Files:**
- Modify: `tests/test_lof_pc.py`

**Interfaces:**
- Test the new function `normalize_molecular_phenotype_id(value: str, line_number: int) -> str`.
- Test the new helper `_collapse_gene_expression_rows(rows: Sequence[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]`.

- [ ] **Step 1: Add canonicalization tests**

Add a parameterized test using these exact Susie-derived cases:

```python
@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("ENSG00000000419.14", "ENSG00000000419"),
        ("A0JNW5_ENSG00000111647.13", "ENSG00000111647"),
        (
            "chr20:50941209:50942031:clu_63027_-:ENSG00000000419.14",
            "ENSG00000000419",
        ),
        (" ENSG1 ", "ENSG1"),
    ],
)
def test_normalize_molecular_phenotype_id(raw: str, expected: str):
    assert lof_pc_module().normalize_molecular_phenotype_id(raw, 7) == expected
```

Add a test asserting that an empty ID, an ID without an Ensembl token, and an ID containing two Ensembl tokens raise `ValueError` with `Line 7` in the message.

- [ ] **Step 2: Add the duplicate-row aggregation test**

Add a test for `_collapse_gene_expression_rows`:

```python
rows = [
    ("ENSG1", np.array([1.0, np.nan, -2.0])),
    ("ENSG1", np.array([0.5, -3.0, np.nan])),
    ("ENSG2", np.array([4.0, 5.0, 6.0])),
]
collapsed = lof_pc_module()._collapse_gene_expression_rows(rows)
assert [gene_id for gene_id, _ in collapsed] == ["ENSG1", "ENSG2"]
np.testing.assert_allclose(
    collapsed[0][1], np.array([0.5, -3.0, -2.0]), equal_nan=True
)
np.testing.assert_allclose(collapsed[1][1], np.array([4.0, 5.0, 6.0]))
```

This verifies that finite values win over missing values and that first-seen gene order is stable.

- [ ] **Step 3: Add an end-to-end Susie-style BED test**

Add a test around the existing `_run_analysis` helper that writes a phenotype BED containing an expression/protein-style ID, two splice-style rows for one gene, and one splice-style row for a second gene:

```text
chr1  0  1  A0JNW5_ENSG1.7                         ...
chr1  0  1  chr1:100:101:clu_1_+:ENSG1.9          ...
chr1  1  2  chr1:200:201:clu_2_-:ENSG2.3          ...
```

Use the existing `ENSG1`/`ENSG2` coding-gene fixture and values with one more-negative duplicate row plus a missing value. Assert that the analysis succeeds and that gene-PC QC contains exactly one row per normalized gene and PC count.

- [ ] **Step 4: Run the new tests before implementation**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k 'molecular_phenotype or collapse or susie'
```

Expected result: FAIL because the new normalizer and collapse helper do not exist, and duplicate normalized IDs are still rejected.

- [ ] **Step 5: Commit the failing tests**

```bash
git add tests/test_lof_pc.py
git commit -m "Test multi-omics phenotype ID handling"
```

### Task 2: Implement canonicalization and duplicate phenotype aggregation

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:197-205, 598-675`
- Test: `tests/test_lof_pc.py`

**Interfaces:**
- Add `normalize_molecular_phenotype_id(value: str, line_number: int) -> str`.
- Add `_collapse_gene_expression_rows(rows: Sequence[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]`.
- Keep `normalize_ensembl_id(value: str) -> str` unchanged.

- [ ] **Step 1: Implement the canonicalizer**

Add a compiled pattern matching one Ensembl token: `(?<![A-Za-z0-9])(ENSG[0-9]+)(?:\.[0-9]+)?`. Strip whitespace, require exactly one match, and return the unversioned captured token. Raise `ValueError` containing `Line {line_number}` and the original ID for empty, unsupported, or ambiguous values.

- [ ] **Step 2: Implement minimum-z row collapsing**

Implement `_collapse_gene_expression_rows` with an insertion-ordered dictionary. Copy the first vector for each gene; combine later vectors with `np.fmin`, retaining `NaN` only when all contributing values are missing. Return `(gene_id, vector)` pairs in first-seen order.

- [ ] **Step 3: Replace duplicate rejection in the BED reader**

In `calculate_lof_pc_enrichment`, replace the current `normalize_ensembl_id(feature_id)` call and duplicate-gene exception with `normalize_molecular_phenotype_id(feature_id, line_number)`. Collect coding rows as `(gene_id, shared-sample NumPy vector)` pairs, then call `_collapse_gene_expression_rows` before checking whether any protein-coding genes remain. Keep existing coordinate, sample, missing-value, and coding-gene validation.

- [ ] **Step 4: Preserve unique-gene downstream behavior**

Use the collapsed list for both complete-data projection and missing-data fallback. Ensure `coding_bed_gene_count` represents unique normalized protein-coding genes after collapsing, so gene-PC QC still contains one row per gene and PC count.

- [ ] **Step 5: Run focused tests and commit**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py -k 'molecular_phenotype or collapse or susie'
```

Expected result: PASS. Then commit:

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py
git commit -m "Support multi-omics phenotype IDs"
```

### Task 3: Add explicit phenotype aggregation QC and documentation

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:882-896`
- Modify: `README.md:18-35, 90-105`
- Test: `tests/test_lof_pc.py`

**Interfaces:**
- Keep `bed_gene_count` and `protein_coding_bed_gene_count` as unique normalized-gene counts.
- Add `bed_feature_count`, `duplicate_feature_count`, `protein_coding_bed_feature_count`, and `protein_coding_duplicate_feature_count` as non-negative integer QC fields.

- [ ] **Step 1: Add QC assertions to the Susie-style test**

For the fixture with three input rows, two unique genes, and one duplicate row, assert:

```python
assert analysis_qc["bed_feature_count"] == 3
assert analysis_qc["bed_gene_count"] == 2
assert analysis_qc["duplicate_feature_count"] == 1
assert analysis_qc["protein_coding_bed_feature_count"] == 3
assert analysis_qc["protein_coding_bed_gene_count"] == 2
assert analysis_qc["protein_coding_duplicate_feature_count"] == 1
```

- [ ] **Step 2: Add the fields to analysis QC**

Track input feature rows and coding feature rows during BED parsing. After collapsing, emit the six counts above. Do not add these fields to `per_pc`; its counters remain PC-specific and merge unchanged.

- [ ] **Step 3: Document supported IDs and collapse semantics**

Update the README molecular phenotype section to document the observed expression, proteomics, and splicing formats, extraction of the Ensembl token, version removal, duplicate splice-row aggregation by sample-wise minimum z-score, and the fact that downstream tests remain gene-level.

- [ ] **Step 4: Run focused tests and commit**

Run:

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc.py
git diff --check
```

Then commit:

```bash
git add src/rare_variant_enrichment/lof_pc.py tests/test_lof_pc.py README.md
git commit -m "Report multi-omics phenotype aggregation QC"
```

### Task 4: Verify the full package and unchanged WDL contract

**Files:**
- Read-only validation: `workflows/rare_variant_enrichment.wdl`
- Read-only validation: `tests/test_wdl_contract.py`

- [ ] **Step 1: Run the full Python test suite**

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q
```

Expected result: all tests pass with only the repository’s existing skips.

- [ ] **Step 2: Validate WDL syntax and the input template**

```bash
/private/tmp/rve-rebuild-venv/bin/miniwdl check workflows/rare_variant_enrichment.wdl
/private/tmp/rve-rebuild-venv/bin/miniwdl input_template workflows/rare_variant_enrichment.wdl
```

Confirm that no new public input was introduced and that the existing ten outputs remain present.

- [ ] **Step 3: Run the existing full fixture CLI path**

```bash
/private/tmp/rve-rebuild-venv/bin/pytest -q tests/test_lof_pc_fixture_end_to_end.py
```

Confirm the original expression fixture still produces its hand-checked Fisher cells and carrier QC counts.

- [ ] **Step 4: Inspect the final diff and working tree**

```bash
git diff --check
git status --short --branch
git log --oneline -4
```

Confirm that only the intended source, tests, README, specification, and plan files changed and that the working tree is clean after commits.
