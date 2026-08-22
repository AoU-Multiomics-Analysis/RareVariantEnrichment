# MT-to-VCF Transcript REVEL Retention Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add gene-matched REVEL scores to the MT-to-VCF `TranscriptAnnotations` output without changing variant filtering.

**Architecture:** Extend the existing transcript schema constants and Hail projection in `filter_and_write_mt.py`. Lock the interface with the lightweight contract tests, verify the actual exported BGZ schema with the pinned Hail integration test, and run both suites in GitHub Actions.

**Tech Stack:** Python 3, Hail, `unittest`, GitHub Actions, WDL 1.0

**Spec:** `../specs/2026-08-22-gene-matched-variant-carrier-extraction-design.md`

## Global Constraints

- Preserve all current `TranscriptAnnotations` columns and add `revel` immediately after `aa_change`.
- Treat `revel` as a numeric Hail `float64` value and export a missing value as an empty TSV field.
- Do not change MT filtering, the VCF schema, or the variant-level annotation schema.
- Do not build a local Docker image for this change.
- Use ASD-STE100-style technical text in documentation and messages.

---

### Task 1: Lock and implement the transcript schema contract

**Files:**
- Modify: `tests/test_filter_and_write_mt_contract.py:14-42`
- Modify: `scripts/filter_and_write_mt.py:83-95`

**Interfaces:**
- Consumes: Existing `TRANSCRIPT_ANNOTATION_FIELDS` and `REQUIRED_TRANSCRIPT_VAT_FIELDS` tuples.
- Produces: A transcript export contract in which `revel` is selected and validated as a required source VAT field.

- [ ] **Step 1: Update the contract test so it requires REVEL**

Change the expected tuple and add a focused source-field assertion:

```python
def test_transcript_output_fields_are_complete(self):
    self.assertEqual(
        MODULE.TRANSCRIPT_ANNOTATION_FIELDS,
        (
            "rsid", "gene_id", "gene_symbol", "transcript",
            "is_canonical_transcript", "consequence", "aa_change", "revel",
            "LoF", "LoF_filter", "LoF_flags", "LoF_info",
            "gvs_max_af", "gvs_max_subpop",
        ),
    )

def test_revel_source_field_is_required(self):
    self.assertIn("revel", MODULE.REQUIRED_TRANSCRIPT_VAT_FIELDS)
```

- [ ] **Step 2: Run the focused contract test and verify failure**

Run:

```bash
python3 -m unittest tests.test_filter_and_write_mt_contract.TranscriptVatContractTests.test_transcript_output_fields_are_complete tests.test_filter_and_write_mt_contract.TranscriptVatContractTests.test_revel_source_field_is_required -v
```

Expected: FAIL because neither tuple contains `revel` in the required position.

- [ ] **Step 3: Add REVEL to both schema tuples**

Update the constants:

```python
TRANSCRIPT_ANNOTATION_FIELDS = (
    "rsid", "gene_id", "gene_symbol", "transcript",
    "is_canonical_transcript", "consequence", "aa_change", "revel",
    "LoF", "LoF_filter", "LoF_flags", "LoF_info",
    "gvs_max_af", "gvs_max_subpop",
)

REQUIRED_TRANSCRIPT_VAT_FIELDS = (
    "vid", "dbsnp_rsid", "gene_id", "gene_symbol", "transcript",
    "is_canonical_transcript", "consequence", "aa_change", "revel",
    "LoF", "LoF_filter", "LoF_flags", "LoF_info",
    "gvs_max_af", "gvs_max_subpop",
)
```

- [ ] **Step 4: Run all lightweight contract tests**

Run:

```bash
python3 -m unittest tests.test_filter_and_write_mt_contract -v
```

Expected: PASS.

- [ ] **Step 5: Commit the schema contract**

```bash
git add scripts/filter_and_write_mt.py tests/test_filter_and_write_mt_contract.py
git commit -m "feat: require REVEL in transcript annotations"
```

---

### Task 2: Export numeric and missing REVEL values through Hail

**Files:**
- Modify: `tests/test_filter_and_write_mt_hail_integration.py:17-207`
- Modify: `scripts/filter_and_write_mt.py:228-244`

**Interfaces:**
- Consumes: The Task 1 `TRANSCRIPT_ANNOTATION_FIELDS` contract and the existing `_cast_to_float` helper inside `_prepare_vat_tables`.
- Produces: `transcript_ht` rows with `revel: tfloat64`, exported in the approved column order.

- [ ] **Step 1: Add REVEL values to the Hail integration fixture**

Insert `"revel"` after `"aa_change"` in `EXPECTED_COLUMNS`. Add `revel="0.82"` to the `ENST0003` missense row and `revel="0.61"` to the `ENST0001` missense row. Leave the `ENST0002` synonymous row empty.

After reading the exported BGZ, add these assertions:

```python
by_transcript = {row["transcript"]: row for row in rows}
self.assertEqual(by_transcript["ENST0003"]["revel"], "0.82")
self.assertEqual(by_transcript["ENST0001"]["revel"], "0.61")
self.assertEqual(by_transcript["ENST0002"]["revel"], "")
```

- [ ] **Step 2: Run the integration test and verify failure in a Hail environment**

Run in the repository's pinned Hail environment:

```bash
python3 -m unittest tests.test_filter_and_write_mt_hail_integration.TranscriptVatHailIntegrationTests.test_composite_key_transcript_export_contract -v
```

Expected: FAIL because `_prepare_vat_tables` does not project `revel` into `transcript_ht`. If Hail is unavailable locally, record the SKIP and rely on the GitHub Actions container job in Task 3. Do not build a local image only for this test.

- [ ] **Step 3: Project REVEL into the transcript Hail table**

Add one field to the existing `transcript_ht` selection expression:

```python
revel=_cast_to_float(vat_source_ht.revel),
```

Place it immediately after `aa_change` so the selected field and exported header order match the contract.

- [ ] **Step 4: Run the lightweight and Hail integration suites**

Run:

```bash
python3 -m unittest tests.test_filter_and_write_mt_contract -v
python3 -m unittest tests.test_filter_and_write_mt_hail_integration -v
```

Expected: Contract tests PASS. Hail integration PASS in the pinned environment or SKIP when Hail is not installed.

- [ ] **Step 5: Commit the Hail projection**

```bash
git add scripts/filter_and_write_mt.py tests/test_filter_and_write_mt_hail_integration.py
git commit -m "feat: export REVEL in transcript annotations"
```

---

### Task 3: Add CI verification and document the schema

**Files:**
- Create: `.github/workflows/python-tests.yml`
- Modify: `README.md:81-124`

**Interfaces:**
- Consumes: The contract and integration suites from Tasks 1 and 2.
- Produces: Pull-request CI that verifies the checked-out script with both standard Python and the published Hail runtime, plus user-facing schema documentation.

- [ ] **Step 1: Add a failing documentation contract**

Add this method to `TranscriptVatContractTests`:

```python
def test_readme_documents_transcript_revel(self):
    readme = (ROOT / "README.md").read_text()
    transcript_section = readme.split("TranscriptAnnotations", 1)[1]
    self.assertIn("`revel`", transcript_section)
    self.assertIn("REVEL", transcript_section)
```

- [ ] **Step 2: Run the documentation contract and verify failure**

Run:

```bash
python3 -m unittest tests.test_filter_and_write_mt_contract.TranscriptVatContractTests.test_readme_documents_transcript_revel -v
```

Expected: FAIL because the transcript-specific section does not list `revel`.

- [ ] **Step 3: Document REVEL in the transcript output**

Update the transcript schema list in `README.md` so it explicitly includes `revel`. State that it is the transcript-level VAT REVEL score, numeric values are exported as decimals, and missing scores are empty fields.

- [ ] **Step 4: Add GitHub Actions tests without a local image build**

Create `.github/workflows/python-tests.yml` with two jobs:

```yaml
name: Python tests

on:
  workflow_dispatch:
  push:
    paths:
      - ".github/workflows/python-tests.yml"
      - "scripts/**"
      - "tests/**"
  pull_request:
    paths:
      - ".github/workflows/python-tests.yml"
      - "scripts/**"
      - "tests/**"

jobs:
  contract:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - run: python -m unittest tests.test_filter_and_write_mt_contract -v

  hail-integration:
    runs-on: ubuntu-latest
    container: ghcr.io/aou-multiomics-analysis/mttovcf:main
    steps:
      - uses: actions/checkout@v4
      - run: python3 -m unittest tests.test_filter_and_write_mt_hail_integration -v
```

This job uses the published Hail image. It does not build an image for the smoke test.

- [ ] **Step 5: Run all locally available tests**

Run:

```bash
python3 -m unittest discover -s tests -v
```

Expected: All available tests PASS. The Hail suite can SKIP only when Hail is absent.

- [ ] **Step 6: Commit CI and documentation**

```bash
git add .github/workflows/python-tests.yml README.md tests/test_filter_and_write_mt_contract.py
git commit -m "ci: verify transcript REVEL export"
```

---

### Task 4: Verify the MT-to-VCF change set

**Files:**
- Verify only; no planned source changes.

**Interfaces:**
- Consumes: All MT-to-VCF tasks above.
- Produces: A verified REVEL-bearing `TranscriptAnnotations` contract for the carrier extractor plan.

- [ ] **Step 1: Run static Python compilation**

```bash
python3 -m compileall -q scripts tests
```

Expected: exit code 0.

- [ ] **Step 2: Run the complete unit-test discovery**

```bash
python3 -m unittest discover -s tests -v
```

Expected: PASS, with only the documented Hail skip allowed outside the pinned environment.

- [ ] **Step 3: Validate WDL files without running Docker**

```bash
miniwdl check main.wdl
miniwdl check workflow/FilterMT.wdl
```

Expected: both checks PASS.

- [ ] **Step 4: Inspect the final diff and repository state**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors and no unintended files.
