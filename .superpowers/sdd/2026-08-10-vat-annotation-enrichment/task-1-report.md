# Task 1 Report: Define VAT Schema and Annotation Semantics

## Implementation

- Added `src/rare_variant_enrichment/annotations.py` using only the Python standard library.
- Added immutable (`frozen=True`) `VatSchema`, `FrequencyValue`, `AnnotationClass`, `VariantKey`, and `GeneAnnotation` value types.
- Implemented exact MTtoVCF required-column resolution, optional `LoF` detection, deterministic missing/duplicate-column validation, schema JSON round-tripping, and duplicate JSON-key rejection.
- Implemented terminal numeric Ensembl gene-version normalization, consequence delimiter parsing/deduplication, the complete specified Ensembl severity order, unknown-term reporting, HC-over-LC LoFTEE collapse, finite `[0, 1]` AF parsing, MAF conversion, and stable annotation-class construction.
- Defined concise canonical error messages and locked them in tests, including deterministic offending columns/classes and specific gene/AF categories.

## Files changed

- `src/rare_variant_enrichment/annotations.py`
- `tests/test_annotations.py`
- `.superpowers/sdd/2026-08-10-vat-annotation-enrichment/task-1-report.md`

## Tests and results

- Focused RED: `.venv/bin/pytest tests/test_annotations.py -q` failed during collection with `ModuleNotFoundError: No module named 'rare_variant_enrichment.annotations'` (exit 2).
- Focused GREEN: `.venv/bin/pytest tests/test_annotations.py -q` passed, `16 passed`.
- Full suite before commit: `.venv/bin/pytest -q` passed, `91 passed, 4 skipped`.
- `git diff --check` passed with no whitespace errors.

## TDD evidence

Tests were added before the production module and observed failing for the expected missing-module reason. The minimal implementation was then added and the focused suite was rerun until all 16 annotation tests passed. The full suite was run once after the focused suite was green.

## Self-review

- Completeness: all interfaces and exact constants in the Task 1 brief are present; schema/value types are immutable; JSON stores the complete header, indices, `lof`, and `lof_enabled`.
- Quality: deterministic ordering is used for required-column errors and duplicate consequence-class errors; allele identity is not transformed; gene normalization is limited to terminal numeric `ENSG` versions; AF non-finite and Boolean inputs have explicit errors.
- Scope: no indexing, storage, workflow, CLI, or non-standard runtime dependency was added.
- Test validity: tests assert public behavior and exact error messages without mocks; focused and full-suite results are fresh.

## Concerns

None identified for Task 1. The canonical error messages are local API decisions authorized by the clarified requirements.

## Commit result

The focused commit was created successfully with message `feat: define VAT annotation semantics`. No commit issue occurred; the report was then included in the amended focused commit.
