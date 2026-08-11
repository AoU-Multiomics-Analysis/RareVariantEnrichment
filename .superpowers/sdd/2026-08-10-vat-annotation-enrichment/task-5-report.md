# Task 5 Report: Chunked VAT-aware chromosome classification

## Status

Implemented and verified.

## Delivered

- Changed `classify_chromosome` to consume indexed VAT data and its `VatSchema`, and to stream VAT and VCF records once per non-overlapping `QueryChunk`.
- Kept transcript annotations disk-backed per chunk with `VatChunkStore`; VCF records, ALT alleles, and genotypes remain record-streamed rather than chromosome-materialized.
- Joined on exact `(chromosome, position, ref, alt)`, required a valid threshold-passing MAF before baseline emission, and matched consequence/LoFTEE annotations to version-stripped feature gene IDs while preserving original feature IDs in output.
- Preserved INFO/AC preference, genotype AC fallback, shared-sample carrier IDs, independent exact/cumulative AC membership, and minimum-distance reduction per annotation family/class.
- Added aggregate chromosome QC for chunk/query counts, VAT ingestion and frequency state, exact allele joins, gene matches, consequence/LoFTEE parsing, and final emitted keys by family. No variant, gene, or sample identifiers are written to ordinary QC.
- Added required classification CLI arguments and transcript-granular indexed VAT fixtures with versioned Ensembl feature IDs.

## TDD and verification

- Red: controlled classifier tests initially failed because the old 10-argument classifier rejected the new 15-argument VAT-aware contract.
- Green: focused Task 5 suite: `18 passed, 2 skipped`.
- Full suite: `120 passed, 5 skipped`.
- Python compilation and `git diff --check` passed.

The two focused skips are the prescribed htslib integration tests in `test_chromosome_classification.py` and `test_end_to_end.py`. The local environment has neither `tabix` nor `bgzip`; controlled pure-Python dual-stream tests ran and passed. The full suite additionally skipped one htslib-dependent VAT test and two WDL runtime tests because the `rare-variant-enrichment:test` Docker image is not built locally.

## Self-review

- Verified every covered coordinate belongs to one chunk and each chunk issues exactly one VAT and one VCF tabix query.
- Verified neighboring-gene frameshift/HC annotations cannot classify the tested gene.
- Verified most-severe transcript collapse, HC-over-LC collapse, AF-to-MAF conversion, common/missing-frequency exclusion, exact-allele mismatch exclusion, multiple classes from different variants, and independent per-class minima.
- Verified empty-feature chromosomes produce header-only outputs and zero aggregate annotation QC without invoking tabix, while strict VCF contig validation remains first.

## Concerns

- Native htslib integration could not run locally because `tabix` and `bgzip` are unavailable; its skip is explicit and controlled coverage exercises the same chunk/join path.
