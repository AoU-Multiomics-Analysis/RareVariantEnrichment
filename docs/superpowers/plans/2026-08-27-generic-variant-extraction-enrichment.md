# Generic Variant Extraction and Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a standalone generic carrier-enrichment WDL that constructs configurable LoF, missense, REVEL, and splice carrier definitions from the existing extraction audit while preserving the legacy LoF workflow.

**Architecture:** Keep `ExtractVariantCarriers` independent and bind its canonical audit to gathered QC with artifact metadata. Add a disk-backed definition builder and refactor the current PC-enrichment engine around ordered `CarrierSets`; the new `CarrierEnrichment` WDL builds definitions first, while the existing LoF WDL calls compatibility adapters to the same engine.

**Tech Stack:** Python 3.12, NumPy 2, SQLite, WDL 1.0, miniwdl, R with tidyverse/ggplot2/ggrepel, deterministic gzip, pytest, GitHub Actions, micromamba.

**Spec:** `docs/superpowers/specs/2026-08-27-generic-variant-extraction-enrichment-design.md`

## Global Constraints

- Keep extraction and enrichment as separate user-facing WDL submissions.
- Do not change either extraction table header or any existing LoF WDL input or output type.
- Do not encode missense or splice carriers in the legacy LoF table.
- Apply no new VCF quality, genotype-quality, AC, AF, MAF, or annotation-frequency filter.
- Definition schema version 1 supports only `lof_hc`, `lof_hc_or_lc`, `missense`, `splice_core`, and `splice_region`, with OR across classes and optional `minimum_revel` AND logic.
- Residualize each gene once per PC count, independent of definition count.
- Recalculate Benjamini-Hochberg FDR across all merged PC, threshold, and definition rows.
- Use bounded SQLite storage for cohort-scale definition construction.
- Use deterministic gzip output and SHA-256 artifact binding.
- Use tidyverse syntax for changed R code. Do not add plot titles or subtitles.
- Every new or changed WDL command must log start, progress or counts, and completion without sample identifiers.
- Use a pinned micromamba base and pinned conda-forge/bioconda runtime packages when changing the Dockerfile.
- Do not build a Docker image locally only for a smoke test. Run the Docker-backed smoke test in GitHub Actions.

---

### Task 1: Bind the Canonical Extraction Audit to Gathered QC

**Files:**
- Create: `src/rare_variant_enrichment/artifacts.py`
- Create: `tests/test_artifacts.py`
- Modify: `src/rare_variant_enrichment/carrier_aggregation.py:32-193`
- Modify: `tests/test_carrier_aggregation.py:41-83`
- Modify: `tests/test_carrier_end_to_end.py:13-77`

**Interfaces:**
- Produces: `file_artifact(path: Path, logical_name: str, header: Sequence[str], row_count: int) -> dict[str, object]`.
- Produces: gathered extraction QC key `audit_artifact` with `logical_name`, `header`, `row_count`, `size_bytes`, and `sha256`.
- Consumes: existing deterministic gathered audit written by `_write_audit`.

- [ ] **Step 1: Write failing artifact tests**

```python
def test_file_artifact_records_exact_bytes_and_schema(tmp_path: Path):
    source = tmp_path / "audit.tsv.gz"
    source.write_bytes(b"abc")
    assert file_artifact(source, "variant_carrier_audit.tsv.gz", ("a", "b"), 3) == {
        "logical_name": "variant_carrier_audit.tsv.gz",
        "header": ["a", "b"],
        "row_count": 3,
        "size_bytes": 3,
        "sha256": hashlib.sha256(b"abc").hexdigest(),
    }

def test_file_artifact_rejects_negative_row_count(tmp_path: Path):
    source = tmp_path / "audit.tsv.gz"
    source.write_bytes(b"")
    with pytest.raises(ValueError, match="row_count"):
        file_artifact(source, "audit", ("a",), -1)
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `python -m pytest tests/test_artifacts.py -v`

Expected: collection fails because `rare_variant_enrichment.artifacts` does not exist.

- [ ] **Step 3: Implement streaming SHA-256 artifact metadata**

```python
def file_artifact(
    path: Path, logical_name: str, header: Sequence[str], row_count: int
) -> dict[str, object]:
    if not logical_name or not header or row_count < 0:
        raise ValueError("Artifact metadata is invalid")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "logical_name": logical_name,
        "header": list(header),
        "row_count": row_count,
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }
```

- [ ] **Step 4: Add the gathered-QC assertion before implementation**

Extend `test_gather_deduplicates_audit_and_aggregates_classes`:

```python
artifact = payload["audit_artifact"]
assert artifact["logical_name"] == "variant_carrier_audit.tsv.gz"
assert artifact["header"] == list(AUDIT_HEADER)
assert artifact["row_count"] == 2
assert artifact["size_bytes"] == audit.stat().st_size
assert artifact["sha256"] == hashlib.sha256(audit.read_bytes()).hexdigest()
```

- [ ] **Step 5: Run the aggregation test and verify RED**

Run: `python -m pytest tests/test_carrier_aggregation.py::test_gather_deduplicates_audit_and_aggregates_classes -v`

Expected: FAIL because `audit_artifact` is absent.

- [ ] **Step 6: Add `audit_artifact` after the audit stream closes**

Call `file_artifact(audit_output, "variant_carrier_audit.tsv.gz", AUDIT_HEADER, audit_count)` when constructing gathered QC. Keep all current QC keys unchanged.

- [ ] **Step 7: Run focused and extraction tests**

Run: `python -m pytest tests/test_artifacts.py tests/test_carrier_aggregation.py tests/test_carrier_end_to_end.py -v`

Expected: PASS, with the indexed end-to-end test allowed to skip outside CI.

- [ ] **Step 8: Commit**

```bash
git add src/rare_variant_enrichment/artifacts.py src/rare_variant_enrichment/carrier_aggregation.py tests/test_artifacts.py tests/test_carrier_aggregation.py tests/test_carrier_end_to_end.py
git commit -m "feat: bind carrier audit to extraction QC"
```

---

### Task 2: Parse and Validate Carrier-Definition Configuration

**Files:**
- Create: `src/rare_variant_enrichment/carrier_definitions.py`
- Create: `tests/test_carrier_definitions.py`

**Interfaces:**
- Produces: `CarrierDefinition(name: str, variant_classes: tuple[str, ...], minimum_revel: float | None)`.
- Produces: `CarrierDefinitionConfig(schema_version: int, definitions: tuple[CarrierDefinition, ...])`.
- Produces: `read_carrier_definition_config(path: Path) -> CarrierDefinitionConfig`.
- Consumes: JSON schema version 1 from the approved specification.

- [ ] **Step 1: Write failing valid-configuration tests**

```python
def test_definition_config_preserves_order_and_threshold(tmp_path: Path):
    config = tmp_path / "definitions.json"
    config.write_text(json.dumps({
        "schema_version": 1,
        "definitions": [
            {"name": "lof_hc", "variant_classes": ["lof_hc"]},
            {
                "name": "missense_revel_ge_0_75",
                "variant_classes": ["missense"],
                "minimum_revel": 0.75,
            },
        ],
    }))
    parsed = read_carrier_definition_config(config)
    assert parsed.names == ("lof_hc", "missense_revel_ge_0_75")
    assert parsed.definitions[1].minimum_revel == 0.75
```

- [ ] **Step 2: Run the valid test and verify RED**

Run: `python -m pytest tests/test_carrier_definitions.py::test_definition_config_preserves_order_and_threshold -v`

Expected: collection fails because the module does not exist.

- [ ] **Step 3: Implement immutable configuration types and strict JSON loading**

```python
SUPPORTED_BASE_CLASSES = frozenset({
    "lof_hc", "lof_hc_or_lc", "missense", "splice_core", "splice_region"
})

@dataclass(frozen=True)
class CarrierDefinition:
    name: str
    variant_classes: tuple[str, ...]
    minimum_revel: float | None = None

    def matches(self, classes: frozenset[str], revel: float | None) -> bool:
        class_match = any(value in classes for value in self.variant_classes)
        threshold_match = self.minimum_revel is None or (
            revel is not None and revel >= self.minimum_revel
        )
        return class_match and threshold_match

@dataclass(frozen=True)
class CarrierDefinitionConfig:
    schema_version: int
    definitions: tuple[CarrierDefinition, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.definitions)
```

Load JSON with an `object_pairs_hook` that rejects duplicate keys. Reject unknown top-level and definition keys, Boolean schema versions and thresholds, invalid names, empty definitions or class lists, duplicates, unsupported classes, and values outside 0 through 1.

- [ ] **Step 4: Add parametrized invalid-configuration tests**

```python
@pytest.mark.parametrize("payload, message", [
    ({"schema_version": True, "definitions": []}, "schema_version"),
    ({"schema_version": 2, "definitions": []}, "Unsupported"),
    ({"schema_version": 1, "definitions": []}, "at least one"),
    ({"schema_version": 1, "definitions": [
        {"name": "1bad", "variant_classes": ["missense"]}
    ]}, "name"),
    ({"schema_version": 1, "definitions": [
        {"name": "bad", "variant_classes": ["unknown"]}
    ]}, "variant class"),
    ({"schema_version": 1, "definitions": [
        {"name": "bad", "variant_classes": ["missense"], "minimum_revel": 1.1}
    ]}, "minimum_revel"),
])
def test_definition_config_rejects_invalid_values(tmp_path: Path, payload, message):
    config = tmp_path / "definitions.json"
    config.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match=message):
        read_carrier_definition_config(config)
```

Add a separate literal duplicate-key fixture because `json.dumps` cannot retain duplicate keys.

- [ ] **Step 5: Run definition tests and verify GREEN**

Run: `python -m pytest tests/test_carrier_definitions.py -v`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/rare_variant_enrichment/carrier_definitions.py tests/test_carrier_definitions.py
git commit -m "feat: validate carrier definition configuration"
```

---

### Task 3: Materialize Generic Carrier Definitions with Bounded Storage

**Files:**
- Modify: `src/rare_variant_enrichment/carrier_definitions.py`
- Modify: `tests/test_carrier_definitions.py`
- Modify: `tests/test_cli.py:1-355`
- Modify: `src/rare_variant_enrichment/cli.py:1-361`

**Interfaces:**
- Consumes: `AUDIT_HEADER`, extraction QC `audit_artifact`, and `CarrierDefinitionConfig` from Tasks 1 and 2.
- Produces: `CARRIER_DEFINITION_HEADER = ("sample_id", "gene_id", "gene_symbol", "carrier_definition", "n_variants", "variant_ids")`.
- Produces: `build_carrier_definitions(audit_path: Path, extraction_qc_path: Path, config_path: Path, output_path: Path, qc_path: Path, *, container_image: str) -> None`.
- Produces: CLI command `build-carrier-definitions`.

- [ ] **Step 1: Write a failing hand-checkable materialization test**

Use audit rows for HC, missense REVEL 0.74, missense REVEL 0.75, splice core, splice region, a duplicate row, and a missing-REVEL row. Assert:

```python
assert [(row["sample_id"], row["gene_id"], row["carrier_definition"], row["n_variants"])
        for row in rows] == [
    ("S1", "ENSG1", "lof_hc", "1"),
    ("S1", "ENSG1", "missense", "2"),
    ("S1", "ENSG1", "missense_revel_ge_0_75", "1"),
    ("S1", "ENSG1", "splice_any", "2"),
]
assert manifest["schema"] == "aou.carrier-definitions-manifest.v1"
assert manifest["definition_order"] == [
    "lof_hc", "missense", "missense_revel_ge_0_75", "splice_any"
]
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python -m pytest tests/test_carrier_definitions.py -k materializes -v`

Expected: FAIL because `build_carrier_definitions` is absent.

- [ ] **Step 3: Implement audit and provenance validation**

Read extraction QC with duplicate-key rejection. Require the exact `audit_artifact` keys and compare logical name, full header, data-row count, byte size, and SHA-256 against the localized audit before opening SQLite.

- [ ] **Step 4: Implement SQLite deduplication and grouping**

Create `audit` and `definition_variant` tables with `WITHOUT ROWID` primary keys. Use this normalized identity:

```sql
PRIMARY KEY (sample_id, gene_id, chrom, pos, ref, alt)
```

For every accepted audit row, evaluate every definition with `definition.matches(classes, revel)`. Insert `(sample_id, gene_id, definition_index, variant_id)` with `INSERT OR IGNORE`. Resolve one non-empty symbol per normalized gene and fail on conflicting non-empty values.

- [ ] **Step 5: Write deterministic gzip output and manifest**

Query by sample ID, gene ID, and definition index. Count distinct variant IDs and join them in lexical order. Calculate the output artifact after close. Include zero counts for every definition and hashes for the audit, extraction QC, configuration, and output.

- [ ] **Step 6: Add failure and determinism tests**

```python
def test_materialization_rejects_audit_digest_mismatch(...):
    audit.write_bytes(audit.read_bytes() + b"x")
    with pytest.raises(ValueError, match="SHA-256"):
        build_carrier_definitions(...)

def test_materialization_is_byte_deterministic(...):
    build_carrier_definitions(..., first, first_qc, container_image="image@sha256:abc")
    build_carrier_definitions(..., second, second_qc, container_image="image@sha256:abc")
    assert first.read_bytes() == second.read_bytes()
```

Also test duplicate conflicts, symbol backfill, zero-match header-only output, missing REVEL, exact threshold inclusion, malformed audit values, and manifest reconciliation.

- [ ] **Step 7: Add CLI dispatch tests before CLI code**

```python
def test_build_carrier_definitions_cli_dispatches(monkeypatch):
    monkeypatch.setattr(cli, "build_carrier_definitions", recorder)
    monkeypatch.setattr(sys, "argv", [
        "rare-variant-enrichment",
        "build-carrier-definitions", "--audit", "audit.tsv.gz",
        "--extraction-qc", "extract.json", "--definitions", "defs.json",
        "--container-image", "image@sha256:abc",
        "--output", "carriers.tsv.gz", "--qc-output", "carriers.json",
    ])
    assert cli.main() == 0
```

- [ ] **Step 8: Run the CLI test and verify RED**

Run: `python -m pytest tests/test_cli.py -k build_carrier_definitions -v`

Expected: FAIL because the command is unknown.

- [ ] **Step 9: Register and dispatch the CLI command**

Add the command to `COMMANDS`, define its exact arguments, and call `build_carrier_definitions` with `container_image` as a keyword-only argument. Log definition and output counts in the implementation without sample IDs.

- [ ] **Step 10: Run materialization and CLI tests**

Run: `python -m pytest tests/test_carrier_definitions.py tests/test_cli.py -v`

Expected: PASS.

- [ ] **Step 11: Commit**

```bash
git add src/rare_variant_enrichment/carrier_definitions.py src/rare_variant_enrichment/cli.py tests/test_carrier_definitions.py tests/test_cli.py
git commit -m "feat: build configurable carrier definitions"
```

---

### Task 4: Refactor PC Enrichment Around Ordered Carrier Sets

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:27-1786`
- Modify: `tests/test_lof_pc.py:313-1580`
- Modify: `src/rare_variant_enrichment/cli.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Consumes: the materialized table and manifest from Task 3.
- Produces: `CarrierSets(definitions: tuple[str, ...], pairs_by_definition: dict[str, set[tuple[str, str]]], qc: dict[str, object])`.
- Produces: `read_generic_carriers(table_path: Path, manifest_path: Path) -> CarrierSets`.
- Produces: `calculate_carrier_pc_enrichment(...) -> None` with the same analysis files as the LoF function plus `carrier_table_path` and `carrier_manifest_path`.
- Preserves: `calculate_lof_pc_enrichment(...) -> None` as a compatibility wrapper.

- [ ] **Step 1: Write failing generic-reader tests**

```python
def test_read_generic_carriers_keeps_declared_zero_definition(tmp_path: Path):
    carriers = read_generic_carriers(table, manifest)
    assert carriers.definitions == ("lof_hc", "missense", "splice_any")
    assert carriers.pairs_by_definition["splice_any"] == set()
    assert carriers.pairs_by_definition["missense"] == {("S1", "ENSG1")}
```

Also assert failure for an undeclared TSV definition, duplicate sample-gene-definition rows, invalid `n_variants`, empty `variant_ids`, output hash mismatch, and manifest count mismatch.

- [ ] **Step 2: Run reader tests and verify RED**

Run: `python -m pytest tests/test_lof_pc.py -k generic_carriers -v`

Expected: FAIL because `read_generic_carriers` is absent.

- [ ] **Step 3: Introduce `CarrierSets` and adapt the LoF reader**

Replace `LofCarriers` with `CarrierSets`. Make `read_lof_carriers` return fixed definitions in exact legacy order and retain its current QC keys and values.

- [ ] **Step 4: Implement manifest-bound generic reading**

Validate the exact six-column table header and artifact metadata. Normalize sample and gene IDs with the existing functions. Initialize every declared definition before reading rows.

- [ ] **Step 5: Extract one shared calculation core**

Use this structure:

```python
def _calculate_pc_enrichment(
    phenotype_bed: Path,
    carriers: CarrierSets,
    carrier_source: Mapping[str, object],
    principal_components_path: Path,
    protein_coding_genes_path: Path,
    negative_z_thresholds: Sequence[float],
    requested_pc_counts: Sequence[int],
    results_output: Path,
    summary_output: Path,
    gene_pc_qc_output: Path,
    analysis_qc_output: Path,
    *,
    pc_grid_mode: str | None,
    additional_covariates_path: Path | None,
    legacy_lof: bool,
) -> None:
    ...
```

Replace every calculation-loop use of `CARRIER_DEFINITIONS` with `carriers.definitions`. Preserve the current complete-data and missing-data residualization paths. Keep legacy JSON keys when `legacy_lof=True`; write generic carrier-table and manifest provenance when false.

- [ ] **Step 6: Prove definition count does not repeat projections**

Add a test that supplies four definitions and monkeypatches `prepare_complete_data_projection` to count calls. Assert the same call count as a one-definition run and assert one result per PC, threshold, and definition.

- [ ] **Step 7: Add hand-checked generic Fisher-cell tests**

Use a small materialized table with LoF, missense, and splice definitions. Assert exact `n11`, `n10`, `n01`, and `n00` values and global FDR across all emitted rows.

- [ ] **Step 8: Add generic CLI tests and commands**

Register `carrier-pc-enrichment` with `--carrier-table` and `--carrier-manifest`, plus the same phenotype, PC, covariate, gene, threshold, grid, and output options as the LoF command. Dispatch to `calculate_carrier_pc_enrichment`.

- [ ] **Step 9: Run generic and legacy calculation tests**

Run: `python -m pytest tests/test_lof_pc.py tests/test_lof_pc_fixture_end_to_end.py tests/test_cli.py -v`

Expected: PASS with unchanged legacy fixture cells.

- [ ] **Step 10: Commit**

```bash
git add src/rare_variant_enrichment/lof_pc.py src/rare_variant_enrichment/cli.py tests/test_lof_pc.py tests/test_lof_pc_fixture_end_to_end.py tests/test_cli.py
git commit -m "feat: run PC enrichment for dynamic carrier definitions"
```

---

### Task 5: Make Merge, PC Selection, SVG, and R Outputs Dynamic

**Files:**
- Modify: `src/rare_variant_enrichment/lof_pc.py:1307-1786`
- Modify: `src/rare_variant_enrichment/pc_selection.py:18-392`
- Modify: `src/rare_variant_enrichment/cli.py`
- Modify: `scripts/pc_sweep_qc.R`
- Modify: `tests/test_lof_pc.py:832-1202`
- Modify: `tests/test_pc_selection.py`
- Modify: `tests/test_cli.py`
- Modify: `tools/lint_r.R`

**Interfaces:**
- Produces: `merge_carrier_pc_enrichment(...) -> None` with dynamic definitions read from shard metadata.
- Produces: `analyze_carrier_pc_enrichment(..., carrier_definitions: Sequence[str], ...) -> None`.
- Preserves: LoF merge and analysis wrappers with fixed legacy behavior.
- Produces: R argument `--carrier-definitions` as an ordered comma-separated list.

- [ ] **Step 1: Write failing dynamic-merge tests**

Create two shards with definitions `("lof_hc", "missense", "splice_any")`. Assert merged ordering, all combinations, global FDR, definition metadata, and rejection when a shard changes definition order or carrier-manifest digest.

- [ ] **Step 2: Run merge tests and verify RED**

Run: `python -m pytest tests/test_lof_pc.py -k 'generic_merge or dynamic_merge' -v`

Expected: FAIL because generic merge is absent or fixed definitions are rejected.

- [ ] **Step 3: Generalize merge validation**

Read the ordered definition list from the first shard. Require a non-empty unique list. Require identical source artifact and manifest metadata for all shards. Use the dynamic list for expected combinations, sorting, merged QC, and carrier counts. Keep the LoF wrapper check for `("any_lof", "HC", "HC_or_LC")`.

- [ ] **Step 4: Write failing selection tests for arbitrary and unestimable definitions**

```python
def test_selection_excludes_zero_carrier_definition_without_dropping_rows(...):
    result = analyze_carrier_pc_enrichment(
        results, selection, plot,
        carrier_definitions=["lof_hc", "missense", "splice_any"],
        selection_z_thresholds=[-3.0, -4.0],
    )
    payload = json.loads(selection.read_text())
    assert payload["selection"]["excluded_definitions"] == {
        "splice_any": "zero_carriers"
    }
```

Add a case where every definition is unestimable and assert `selected_pc_count is None` with valid files.

- [ ] **Step 5: Run selection tests and verify RED**

Run: `python -m pytest tests/test_pc_selection.py -k 'dynamic or unestimable or zero_carrier' -v`

Expected: FAIL because current selection requires fixed complete curves.

- [ ] **Step 6: Generalize Python selection and SVG layout**

Require an explicit non-empty ordered definition list in the generic entry point. Calculate panel count from list length, preserve fixed LoF labels for the legacy wrapper, omit selection lines when the selected PC is null, and emit exclusion reasons.

- [ ] **Step 7: Convert the R script to tidyverse syntax**

Parse `--carrier-definitions`, validate order and uniqueness, use `readr::read_tsv`, `dplyr::filter/group_by/summarise/arrange`, and `ggplot2`. Facet in supplied definition order. Retain the minimal plot theme and blank title/subtitle/caption.

- [ ] **Step 8: Add CLI commands and tests**

Register `merge-carrier-pc-enrichment` and `analyze-carrier-pc-enrichment`. The analysis command requires `--carrier-definitions`. Keep all LoF commands and defaults unchanged.

- [ ] **Step 9: Run merge, selection, CLI, and R lint tests**

Run: `python -m pytest tests/test_lof_pc.py tests/test_pc_selection.py tests/test_cli.py -v`

Run: `Rscript tools/lint_r.R scripts/pc_sweep_qc.R`

Expected: PASS.

- [ ] **Step 10: Commit**

```bash
git add src/rare_variant_enrichment/lof_pc.py src/rare_variant_enrichment/pc_selection.py src/rare_variant_enrichment/cli.py scripts/pc_sweep_qc.R tests/test_lof_pc.py tests/test_pc_selection.py tests/test_cli.py tools/lint_r.R
git commit -m "feat: generalize carrier merge and PC selection"
```

---

### Task 6: Add the Generic CarrierEnrichment WDL

**Files:**
- Create: `workflows/carrier_enrichment.wdl`
- Create: `tests/test_carrier_enrichment_wdl_contract.py`
- Create: `tests/test_carrier_enrichment_wdl_runtime.py`
- Modify: `.dockstore.yml`

**Interfaces:**
- Consumes: extraction audit, extraction QC, definition JSON, phenotype BED, PCs, optional covariates, GTF, and runtime settings from the spec.
- Produces: twelve public files from `CarrierEnrichment`.
- Uses: `build-carrier-definitions`, generic calculate, generic merge, generic analyze, and dynamic R plotting commands.

- [ ] **Step 1: Write a failing exact WDL-contract test**

Load `workflows/carrier_enrichment.wdl` with `WDL.load`. Assert the exact public input names and types from the spec, these tasks, and twelve outputs:

```python
assert sorted(task.name for task in document.tasks) == [
    "AnalyzeCarrierPcEnrichment",
    "BuildCarrierDefinitions",
    "CalculateCarrierPcEnrichment",
    "MergeCarrierPcEnrichment",
    "PreparePcChunks",
    "PrepareProteinCodingGenes",
]
assert "scatter (pc_count_chunk in pc_count_chunks)" in source
assert "scatter (carrier" not in source
```

Assert that commands quote files, materialize arrays through `write_lines`, pass full container provenance, and contain start/count/completion logging.

- [ ] **Step 2: Run the contract test and verify RED**

Run: `python -m pytest tests/test_carrier_enrichment_wdl_contract.py -v`

Expected: FAIL because the WDL does not exist.

- [ ] **Step 3: Implement the WDL tasks and workflow**

Start from the established task patterns in `rare_variant_enrichment.wdl`. Add `BuildCarrierDefinitions` before PC chunk preparation. Calculate materialization disk from audit, extraction QC, and definition JSON sizes. Calculate analysis disk from phenotype, generic carrier table, manifest, PCs, GTF, and optional covariates. Scatter only on PC chunks and merge globally.

- [ ] **Step 4: Add the runtime fixture test**

Build an inputs JSON with fixture audit, extraction QC, definition JSON, phenotype BED, PCs, and GTF. Run miniwdl with the CI test image. Assert all twelve outputs and exact Fisher cells for at least LoF, missense, and splice definitions.

- [ ] **Step 5: Register the workflow**

Append this Dockstore entry without changing the existing two:

```yaml
  - subclass: WDL
    primaryDescriptorPath: /workflows/carrier_enrichment.wdl
```

- [ ] **Step 6: Validate and run WDL contract tests**

Run: `miniwdl check workflows/carrier_enrichment.wdl`

Run: `python -m pytest tests/test_carrier_enrichment_wdl_contract.py -v`

Expected: PASS. The runtime test can skip locally when the CI image is unavailable.

- [ ] **Step 7: Re-run the legacy WDL contract**

Run: `python -m pytest tests/test_wdl_contract.py tests/test_carrier_wdl_contract.py -v`

Expected: PASS with no legacy interface change.

- [ ] **Step 8: Commit**

```bash
git add workflows/carrier_enrichment.wdl tests/test_carrier_enrichment_wdl_contract.py tests/test_carrier_enrichment_wdl_runtime.py .dockstore.yml
git commit -m "feat: add generic carrier enrichment WDL"
```

---

### Task 7: Pin the Micromamba Runtime and Update Examples and Documentation

**Files:**
- Modify: `envs/Dockerfile`
- Create: `examples/carrier_definitions.json`
- Create: `examples/carrier_enrichment.inputs.json`
- Modify: `README.md`
- Modify: `.github/workflows/python-tests.yml`
- Modify: `tests/test_carrier_enrichment_wdl_contract.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Produces: a pinned micromamba runtime with Python, NumPy, htslib, R tidyverse, ggplot2, and ggrepel.
- Produces: runnable example contracts for the two-workflow handoff.
- Preserves: CI as the only Docker-backed smoke-test location.

- [ ] **Step 1: Add failing static runtime and example tests**

Assert that the Dockerfile starts from a pinned micromamba tag, uses both `conda-forge` and `bioconda`, pins every direct runtime package, and contains no `apt-get`. Parse both new JSON examples and assert exact workflow-qualified keys. Assert README contains extraction-to-enrichment commands and the prefiltered-VCF rule.

- [ ] **Step 2: Run static tests and verify RED**

Run: `python -m pytest tests/test_carrier_enrichment_wdl_contract.py tests/test_cli.py -k 'runtime_image or example or readme' -v`

Expected: FAIL because the image and examples are not updated.

- [ ] **Step 3: Replace the runtime image definition**

Use `mambaorg/micromamba:2.8.1` as the base. Configure strict channel priority with `conda-forge` before `bioconda`. Install `python=3.12.14`, `numpy=2.5.2`, `htslib=1.24`, `r-base=4.6.1`, `r-tidyverse=2.0.0`, `r-ggplot2=4.0.3`, and `r-ggrepel=0.9.8`. Install the repository with `/opt/conda/bin/python -m pip install --no-deps .` so the Python package does not start a second dependency solve.

- [ ] **Step 4: Add definition and WDL input examples**

`examples/carrier_definitions.json` uses the approved seven definitions. `examples/carrier_enrichment.inputs.json` supplies audit, extraction QC, definition JSON, phenotype BED, PCs, optional covariates, GTF, PC settings, resources, and a container reference.

- [ ] **Step 5: Update README and CI**

Document two separate submissions and the carrier audit handoff. Define “build carrier definitions” in plain language. Document base classes, REVEL inclusion, overlapping definitions, all twelve outputs, global FDR, and compatibility behavior. Ensure the Python workflow validates every WDL and runs both legacy and generic Docker-backed runtime tests after its one CI image build.

- [ ] **Step 6: Run static, CLI, and documentation tests**

Run: `python -m pytest tests/test_carrier_enrichment_wdl_contract.py tests/test_cli.py -v`

Run: `python -m json.tool examples/carrier_definitions.json`

Run: `python -m json.tool examples/carrier_enrichment.inputs.json`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add envs/Dockerfile examples/carrier_definitions.json examples/carrier_enrichment.inputs.json README.md .github/workflows/python-tests.yml tests/test_carrier_enrichment_wdl_contract.py tests/test_cli.py
git commit -m "build: add generic enrichment runtime and examples"
```

---

### Task 8: Complete End-to-End Compatibility and Verification

**Files:**
- Create: `tests/test_generic_carrier_enrichment_end_to_end.py`
- Modify: `tests/fixtures/transcript_annotations.tsv`
- Modify: `tests/fixtures/rare_variants.vcf`
- Modify: `tests/fixtures/lof_pc_phenotypes.bed` only if the generic hand-checked cells require an additional gene row
- Modify: relevant implementation and tests only when a failing end-to-end assertion exposes a defect

**Interfaces:**
- Consumes: extraction audit contract, definition builder, generic engine, dynamic selection, and the new WDL.
- Produces: one hand-checkable test from carrier audit through all generic enrichment outputs.
- Verifies: unchanged legacy LoF values and schemas.

- [ ] **Step 1: Write the failing complete Python end-to-end test**

Create an audit and extraction QC from fixtures, build definitions, run generic enrichment for PC counts 0 and 1, merge if needed, and analyze selection. Include HC, LC, missense below and equal to 0.75, splice core, splice region, and overlapping `splice_any`. Assert exact materialized rows, exact Fisher cells, result count, global FDR, definition order, manifest digests, QC reconciliation, and plot files.

- [ ] **Step 2: Run the end-to-end test and verify RED or expose gaps**

Run: `python -m pytest tests/test_generic_carrier_enrichment_end_to_end.py -v`

Expected: FAIL on the first missing or inconsistent integrated behavior.

- [ ] **Step 3: Reconcile the integrated contracts**

For each failing assertion, assign it to the owning contract: audit provenance in `artifacts.py` or `carrier_aggregation.py`; definition parsing or materialization in `carrier_definitions.py`; calculation or merge in `lof_pc.py`; selection or SVG in `pc_selection.py`; or command dispatch in `cli.py`. Preserve the failing assertion, change only the owning implementation, and rerun that assertion before the complete end-to-end test. Do not change the approved schemas.

- [ ] **Step 4: Run all focused suites**

Run:

```bash
python -m pytest \
  tests/test_artifacts.py \
  tests/test_carrier_aggregation.py \
  tests/test_carrier_definitions.py \
  tests/test_lof_pc.py \
  tests/test_pc_selection.py \
  tests/test_cli.py \
  tests/test_carrier_wdl_contract.py \
  tests/test_carrier_enrichment_wdl_contract.py \
  tests/test_generic_carrier_enrichment_end_to_end.py -v
```

Expected: PASS.

- [ ] **Step 5: Run repository-wide verification**

Run: `git diff --check`

Run: `python -m pytest -q -rs`

Run: `miniwdl check workflows/extract_variant_carriers.wdl`

Run: `miniwdl check workflows/carrier_enrichment.wdl`

Run: `miniwdl check workflows/rare_variant_enrichment.wdl`

Run: `Rscript tools/lint_r.R scripts/pc_sweep_qc.R`

Expected: all available local checks pass. Docker-backed runtime tests can skip locally and must run in GitHub Actions.

- [ ] **Step 6: Confirm the legacy public contract explicitly**

Run: `python -m pytest tests/test_lof_pc_fixture_end_to_end.py tests/test_wdl_contract.py tests/test_wdl_runtime.py -q -rs`

Expected: legacy Python and contract tests pass; local Docker runtime tests can skip only for the documented missing image or daemon condition.

- [ ] **Step 7: Commit final integration fixes**

```bash
git add tests/test_generic_carrier_enrichment_end_to_end.py tests/fixtures src workflows scripts README.md examples envs .github .dockstore.yml
git commit -m "test: verify generic carrier enrichment end to end"
```

- [ ] **Step 8: Push and inspect GitHub Actions**

Run: `git push -u origin codex/generic-variant-enrichment`

Run: `gh run list --branch codex/generic-variant-enrichment --limit 10`

Wait for the branch workflows. If a job fails, inspect its logs, reproduce the defect with a failing local test when possible, fix it test-first, push the fix, and wait again.

- [ ] **Step 9: Final branch audit**

Run: `git status --short --branch`

Run: `git log --oneline origin/main..HEAD`

Expected: clean worktree, branch tracks `origin/codex/generic-variant-enrichment`, and all intended commits are present.
