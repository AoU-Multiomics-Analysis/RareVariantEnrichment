# Gene-Matched Variant Carrier Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone chromosome-scattered WDL that converts a filtered VCF and transcript VAT annotations into a gene-matched carrier audit, aggregated variant-class table, and reconciled QC.

**Architecture:** Add focused annotation, chunk-storage, extraction, and gather modules beside the existing distance-based code. The new workflow performs an exact allele join, collapses transcript annotations per normalized gene, emits all carrier variants to a detailed audit, and derives the approved initial classes without calling enrichment.

**Tech Stack:** Python 3.12, SQLite, gzip TSV, `tabix`/`bgzip`, WDL 1.0, miniwdl, pytest, GitHub Actions

**Spec:** `../specs/2026-08-22-gene-matched-variant-carrier-extraction-design.md`

## Global Constraints

- Treat the input VCF as prefiltered. Do not apply another quality, AC, AF, or MAF filter.
- Match VCF and transcript rows only by exact chromosome, position, REF, and ALT.
- Match carriers to normalized VAT gene IDs. Do not use TSS distance, phenotype position, a GTF, or a feature list.
- Select one most-severe consequence per variant-gene pair, LoFTEE `HC` over `LC`, maximum finite REVEL, and maximum valid `gvs_max_af`.
- Preserve the active `workflows/rare_variant_enrichment.wdl` interface and behavior.
- Make outputs deterministic and gzip-compressed.
- Add explicit start, progress, count, and completion messages to every new WDL command.
- Do not log sample-level audit records.
- Do not build a local Docker image for smoke testing. Run the Docker-backed smoke test in GitHub Actions.
- Preserve unrelated existing workspace changes and stage explicit paths for every commit.
- Use ASD-STE100-style technical text in documentation, logs, and errors.

---

### Task 1: Add the transcript schema and pure annotation collapse rules

**Files:**
- Create: `src/rare_variant_enrichment/carrier_annotations.py`
- Create: `tests/test_carrier_annotations.py`

**Interfaces:**
- Consumes: `VariantKey`, `ENSEMBL_CONSEQUENCE_ORDER`, `normalize_gene_id`, `parse_consequence_terms`, `most_severe_consequence`, and `collapse_loftee` from `rare_variant_enrichment.annotations`.
- Produces: `TranscriptCarrierSchema`, `TranscriptCarrierRow`, `CollapsedCarrierAnnotation`, `collapse_transcript_rows()`, and `initial_variant_classes()` for storage and extraction tasks.

- [ ] **Step 1: Write failing schema and numeric-value tests**

Create `tests/test_carrier_annotations.py` with a complete header fixture and these first tests:

```python
from rare_variant_enrichment.carrier_annotations import (
    TranscriptCarrierSchema,
    parse_optional_unit_interval,
)

HEADER = (
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "revel", "LoF", "LoF_filter", "LoF_flags", "LoF_info",
    "gvs_max_af", "gvs_max_subpop",
)

def test_transcript_carrier_schema_resolves_required_columns():
    schema = TranscriptCarrierSchema.from_header(HEADER)
    assert schema.gene_id == HEADER.index("gene_id")
    assert schema.gene_symbol == HEADER.index("gene_symbol")
    assert schema.revel == HEADER.index("revel")

def test_optional_unit_interval_accepts_missing_and_bounds():
    assert parse_optional_unit_interval("", "REVEL") is None
    assert parse_optional_unit_interval(".", "REVEL") is None
    assert parse_optional_unit_interval("0", "REVEL") == 0.0
    assert parse_optional_unit_interval("1", "REVEL") == 1.0
```

Add parameterized failures for duplicate required columns, missing `revel`, `nan`, `-0.1`, `1.1`, and nonnumeric values.

- [ ] **Step 2: Run the schema tests and verify import failure**

Run:

```bash
python3 -m pytest -q tests/test_carrier_annotations.py
```

Expected: collection FAILS because `carrier_annotations` does not exist.

- [ ] **Step 3: Implement the schema and primitive parsers**

Create these public types and signatures:

```python
@dataclass(frozen=True)
class TranscriptCarrierSchema:
    header: tuple[str, ...]
    chromosome: int
    position: int
    ref: int
    alt: int
    gene_id: int
    gene_symbol: int
    consequence: int
    lof: int
    gvs_max_af: int
    revel: int

    @classmethod
    def from_header(cls, header: Sequence[str]) -> "TranscriptCarrierSchema":
        required = (
            "chrom", "pos", "ref", "alt", "gene_id", "gene_symbol",
            "consequence", "LoF", "gvs_max_af", "revel",
        )
        missing = [name for name in required if name not in header]
        duplicates = [name for name in required if header.count(name) > 1]
        if missing:
            raise ValueError("Missing required transcript columns: " + ", ".join(missing))
        if duplicates:
            raise ValueError("Duplicate required transcript columns: " + ", ".join(duplicates))
        index = {name: header.index(name) for name in required}
        return cls(
            tuple(header), index["chrom"], index["pos"], index["ref"], index["alt"],
            index["gene_id"], index["gene_symbol"], index["consequence"],
            index["LoF"], index["gvs_max_af"], index["revel"],
        )

    def as_dict(self) -> dict[str, object]:
        return {"header": list(self.header), **{
            name: getattr(self, name) for name in (
                "chromosome", "position", "ref", "alt", "gene_id",
                "gene_symbol", "consequence", "lof", "gvs_max_af", "revel",
            )
        }}

    def write_json(self, path: Path) -> None:
        write_json(path, self.as_dict())

    @classmethod
    def read_json(cls, path: Path) -> "TranscriptCarrierSchema":
        payload = read_json_object_without_duplicate_keys(path)
        resolved = cls.from_header(payload["header"])
        if resolved.as_dict() != payload:
            raise ValueError("Transcript schema indices do not match its header")
        return resolved

@dataclass(frozen=True)
class TranscriptCarrierRow:
    key: VariantKey
    gene_id: str
    gene_symbol: str | None
    consequences: tuple[str, ...]
    loftee: str | None
    revel: float | None
    gvs_max_af: float | None

def parse_optional_unit_interval(value: str, label: str) -> float | None:
    stripped = value.strip()
    if stripped.casefold() in {"", ".", "na", "null"}:
        return None
    try:
        parsed = float(stripped)
    except ValueError as error:
        raise ValueError(f"{label} must be numeric or missing") from error
    if not math.isfinite(parsed) or not 0.0 <= parsed <= 1.0:
        raise ValueError(f"{label} must be from 0 through 1")
    return parsed

def parse_transcript_carrier_row(
    fields: Sequence[str], schema: TranscriptCarrierSchema
) -> TranscriptCarrierRow:
    if len(fields) != len(schema.header):
        raise ValueError("Transcript row column count does not match its header")
    key = VariantKey(
        fields[schema.chromosome], int(fields[schema.position]),
        fields[schema.ref], fields[schema.alt],
    )
    gene_symbol = fields[schema.gene_symbol].strip()
    raw_loftee = fields[schema.lof].strip().upper()
    return TranscriptCarrierRow(
        key=key,
        gene_id=normalize_gene_id(fields[schema.gene_id].strip()),
        gene_symbol=gene_symbol if gene_symbol not in {"", "."} else None,
        consequences=parse_consequence_terms(fields[schema.consequence]),
        loftee=raw_loftee if raw_loftee in {"HC", "LC"} else None,
        revel=parse_optional_unit_interval(fields[schema.revel], "REVEL"),
        gvs_max_af=parse_optional_unit_interval(fields[schema.gvs_max_af], "gvs_max_af"),
    )
```

Use exact required-column names. Reject duplicate required columns. Accept `""`, `"."`, `"NA"`, and `"null"` as missing. Reject `"NaN"`, infinities, and other non-finite numeric values.

- [ ] **Step 4: Run schema tests to verify they pass**

```bash
python3 -m pytest -q tests/test_carrier_annotations.py
```

Expected: PASS for the schema and parser cases added so far.

- [ ] **Step 5: Write failing collapse and class-assignment tests**

Add tests that construct `TranscriptCarrierRow` values for one exact allele and normalized gene. Cover:

```python
def test_collapse_selects_most_severe_hc_and_maximum_scores():
    collapsed = collapse_transcript_rows([
        row("ENSG1.1", "GENE1", "missense_variant", "LC", 0.42, 0.001),
        row("ENSG1.2", "GENE1", "splice_acceptor_variant", "HC", 0.81, 0.004),
        row("ENSG1.3", "GENE1", "synonymous_variant", None, None, None),
    ])
    assert collapsed.gene_id == "ENSG1"
    assert collapsed.most_severe_consequence == "splice_acceptor_variant"
    assert collapsed.all_consequences == (
        "splice_acceptor_variant", "missense_variant", "synonymous_variant"
    )
    assert collapsed.loftee == "HC"
    assert collapsed.revel == 0.81
    assert collapsed.gvs_max_af == 0.004
    assert initial_variant_classes(collapsed) == (
        "lof_hc", "lof_hc_or_lc", "splice_core"
    )
```

Add separate cases for LC-only, missense, each `splice_region` term, no initial class, conflicting symbols, and rows from different alleles or genes.

- [ ] **Step 6: Run the collapse tests and verify failure**

```bash
python3 -m pytest -q tests/test_carrier_annotations.py -k 'collapse or class'
```

Expected: FAIL because the collapse types and functions do not exist.

- [ ] **Step 7: Implement collapse and initial class assignment**

Add these interfaces:

```python
@dataclass(frozen=True)
class CollapsedCarrierAnnotation:
    key: VariantKey
    gene_id: str
    gene_symbol: str | None
    most_severe_consequence: str | None
    all_consequences: tuple[str, ...]
    unknown_consequences: tuple[str, ...]
    loftee: str | None
    revel: float | None
    gvs_max_af: float | None

def collapse_transcript_rows(
    rows: Sequence[TranscriptCarrierRow],
) -> CollapsedCarrierAnnotation:
    if not rows:
        raise ValueError("At least one transcript row is required")
    keys = {row.key for row in rows}
    genes = {row.gene_id for row in rows}
    symbols = {row.gene_symbol for row in rows if row.gene_symbol is not None}
    if len(keys) != 1 or len(genes) != 1:
        raise ValueError("Transcript collapse requires one variant-gene pair")
    if len(symbols) > 1:
        raise ValueError("Conflicting gene symbols for normalized gene ID")
    selected, unknown = most_severe_consequence(
        term for row in rows for term in row.consequences
    )
    recognized = sorted(
        {term for row in rows for term in row.consequences if term not in unknown},
        key=ENSEMBL_CONSEQUENCE_ORDER.index,
    )
    revel_values = [row.revel for row in rows if row.revel is not None]
    frequency_values = [row.gvs_max_af for row in rows if row.gvs_max_af is not None]
    return CollapsedCarrierAnnotation(
        key=next(iter(keys)), gene_id=next(iter(genes)),
        gene_symbol=next(iter(symbols), None),
        most_severe_consequence=selected,
        all_consequences=tuple(recognized),
        unknown_consequences=tuple(sorted(set(unknown))),
        loftee=collapse_loftee(row.loftee or "" for row in rows),
        revel=max(revel_values, default=None),
        gvs_max_af=max(frequency_values, default=None),
    )

def initial_variant_classes(
    annotation: CollapsedCarrierAnnotation,
) -> tuple[str, ...]:
    classes = []
    if annotation.loftee == "HC":
        classes.append("lof_hc")
    if annotation.loftee in {"HC", "LC"}:
        classes.append("lof_hc_or_lc")
    if annotation.most_severe_consequence == "missense_variant":
        classes.append("missense")
    if annotation.most_severe_consequence in SPLICE_CORE_TERMS:
        classes.append("splice_core")
    if annotation.most_severe_consequence in SPLICE_REGION_TERMS:
        classes.append("splice_region")
    return tuple(classes)
```

Sort recognized consequence terms by `ENSEMBL_CONSEQUENCE_ORDER`. Sort unknown terms lexically. Emit classes in this fixed order: `lof_hc`, `lof_hc_or_lc`, `missense`, `splice_core`, `splice_region`.

Define the duplicate-key JSON reader used by `TranscriptCarrierSchema.read_json()` in the same module:

```python
def read_json_object_without_duplicate_keys(path: Path) -> dict[str, object]:
    def reject_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Transcript schema JSON contains duplicate key: {key}")
            result[key] = value
        return result

    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle, object_pairs_hook=reject_duplicates)
    if not isinstance(payload, dict):
        raise ValueError("Transcript schema JSON must contain an object")
    return payload
```

- [ ] **Step 8: Run annotation tests and regression tests**

```bash
python3 -m pytest -q tests/test_carrier_annotations.py tests/test_annotations.py
```

Expected: PASS.

- [ ] **Step 9: Commit the annotation model**

```bash
git add src/rare_variant_enrichment/carrier_annotations.py tests/test_carrier_annotations.py
git commit -m "feat: define carrier annotation collapse rules"
```

---

### Task 2: Extend VCF parsing with carrier ALT dosage and AF

**Files:**
- Modify: `src/rare_variant_enrichment/variants.py:49-165,500-618`
- Modify: `tests/test_variants.py:1-165`
- Modify: `tests/fixtures/rare_variants.vcf:1-8`

**Interfaces:**
- Consumes: Existing `parse_variant_alleles()` behavior used by the distance-based workflow.
- Produces: Backward-compatible `VariantAllele.carriers`, plus `VariantAllele.carrier_alt_counts` and `VariantAllele.af` for the new extractor.

- [ ] **Step 1: Add failing dosage and AF tests**

Add tests with `GT` values `0/1`, `1|1`, `0/2`, partial calls, and ALT-specific `AC`/`AF` arrays:

```python
def test_parse_variant_alleles_retains_dosage_and_alt_specific_af():
    fields = (
        "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=3,1;AF=0.3,0.1\tGT\t"
        "0/1\t1|1\t0/2"
    ).split("\t")
    alleles = parse_variant_alleles(fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"})
    assert alleles[0].carrier_alt_counts == (("S1", 1), ("S2", 2))
    assert alleles[0].carriers == ("S1", "S2")
    assert alleles[0].af == 0.3
    assert alleles[1].carrier_alt_counts == (("S3", 1),)
    assert alleles[1].af == 0.1
```

Add failures for wrong INFO/AF cardinality, negative AF, AF above 1, and nonnumeric AF.

- [ ] **Step 2: Run the focused parser tests and verify failure**

```bash
python3 -m pytest -q tests/test_variants.py -k 'dosage or alt_specific_af'
```

Expected: FAIL because `VariantAllele` has no dosage or AF fields.

- [ ] **Step 3: Extend `VariantAllele` without breaking `.carriers`**

Use this public shape:

```python
@dataclass(frozen=True)
class VariantAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    ac: int
    af: float | None
    carrier_alt_counts: tuple[tuple[str, int], ...]

    @property
    def carriers(self) -> tuple[str, ...]:
        return tuple(sample_id for sample_id, _ in self.carrier_alt_counts)
```

Change `_count_genotypes()` so it returns one ordered `(sample_id, dosage)` list for each ALT. Add `_info_af_values(info: str, alt_count: int) -> list[float | None]`. Keep all current AC fallback and QC behavior.

- [ ] **Step 4: Update the fixture with production-like metadata**

Change the contig header to `##contig=<ID=chr1,length=1000>`, add an INFO/AF declaration, and add ALT-specific `AF` values to each fixture record. Do not change genotypes or AC values.

- [ ] **Step 5: Run VCF and legacy classification tests**

```bash
python3 -m pytest -q tests/test_variants.py tests/test_chromosome_classification.py tests/test_end_to_end.py
```

Expected: PASS. Existing `.carriers` assertions remain unchanged.

- [ ] **Step 6: Commit the reusable VCF carrier representation**

```bash
git add src/rare_variant_enrichment/variants.py tests/test_variants.py tests/fixtures/rare_variants.vcf
git commit -m "feat: retain carrier dosage and allele frequency"
```

---

### Task 3: Add bounded transcript-annotation chunk storage

**Files:**
- Create: `src/rare_variant_enrichment/carrier_annotation_storage.py`
- Create: `tests/test_carrier_annotation_storage.py`

**Interfaces:**
- Consumes: `TranscriptCarrierSchema`, `parse_transcript_carrier_row()`, and `collapse_transcript_rows()` from Task 1.
- Produces: `CarrierAnnotationChunkStore.ingest()`, `.finalize()`, `.has_allele()`, and `.annotations_for()` for one coordinate chunk.

- [ ] **Step 1: Write failing storage lifecycle and collapse tests**

Create tests that ingest duplicate transcript rows, two genes on one allele, multiple transcript consequences for one gene, missing REVEL, and a conflicting symbol. Assert this interface:

```python
with CarrierAnnotationChunkStore(tmp_path, schema) as store:
    store.ingest(fields_one)
    store.ingest(fields_one)
    store.ingest(fields_two)
    qc = store.finalize()
    annotations = store.annotations_for(VariantKey("chr1", 100, "A", "C"))

assert qc["transcript_rows"] == 3
assert qc["duplicate_transcript_rows"] == 1
assert [item.gene_id for item in annotations] == ["ENSG1", "ENSG2"]
```

Also assert that ingest after finalize and query before finalize fail with explicit messages.

- [ ] **Step 2: Run storage tests and verify import failure**

```bash
python3 -m pytest -q tests/test_carrier_annotation_storage.py
```

Expected: collection FAILS because the module does not exist.

- [ ] **Step 3: Implement the SQLite chunk store**

Create `CarrierAnnotationChunkStore` with constructor arguments `(directory: Path, schema: TranscriptCarrierSchema)`. Implement context-manager methods, `ingest(fields: Sequence[str]) -> None`, `finalize() -> dict[str, int | float | str]`, `has_allele(key: VariantKey) -> bool`, and `annotations_for(key: VariantKey) -> tuple[CollapsedCarrierAnnotation, ...]`.

Use this table and index:

```sql
CREATE TABLE transcript_row (
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    gene_symbol TEXT,
    consequences_json TEXT NOT NULL,
    loftee TEXT,
    revel REAL,
    gvs_max_af REAL,
    row_fingerprint TEXT NOT NULL UNIQUE
);
CREATE INDEX transcript_row_allele_gene
ON transcript_row (chromosome, position, ref, alt, gene_id);
```

`ingest()` must parse with `parse_transcript_carrier_row()`, serialize a canonical JSON fingerprint, and use `INSERT OR IGNORE`. `finalize()` commits and freezes ingestion. `annotations_for()` selects rows ordered by gene ID and fingerprint, groups one gene at a time, and calls `collapse_transcript_rows()`. Delete the temporary database on close.

- [ ] **Step 4: Add QC assertions**

Cover unique alleles, unique allele-gene pairs, recognized and unknown consequence terms, LoFTEE counts, REVEL present/missing counts and range, `gvs_max_af` present/missing counts and range, and selected consequence counts.

- [ ] **Step 5: Run storage and annotation tests**

```bash
python3 -m pytest -q tests/test_carrier_annotation_storage.py tests/test_carrier_annotations.py
```

Expected: PASS.

- [ ] **Step 6: Commit bounded annotation storage**

```bash
git add src/rare_variant_enrichment/carrier_annotation_storage.py tests/test_carrier_annotation_storage.py
git commit -m "feat: add bounded carrier annotation storage"
```

---

### Task 4: Extract one chromosome without distance filtering

**Files:**
- Create: `src/rare_variant_enrichment/carrier_extraction.py`
- Create: `tests/test_carrier_extraction.py`
- Modify: `tests/fixtures/transcript_annotations.tsv:1-7`
- Modify: `tests/conftest.py:1-85`

**Interfaces:**
- Consumes: Task 1 annotation types, Task 2 `VariantAllele`, Task 3 chunk store, existing tabix streaming conventions, and `write_json()`.
- Produces: `prepare_carrier_inputs()`, `build_chromosome_chunks()`, and `extract_chromosome_carriers()`.

- [ ] **Step 1: Add REVEL to the transcript fixture**

Insert `revel` after `aa_change`. Use values that test maximum selection: give the two `chr1:100:A:C` `GENE1` transcripts `0.42` and `0.81`; use `0.65` for the `chr1:150:G:T` missense row; leave non-missense values empty.

Update `PreparedFixture` so it creates and returns a `TranscriptCarrierSchema` JSON file in addition to the existing legacy `VatSchema` file. Keep the legacy fixture fields so existing tests do not change behavior.

- [ ] **Step 2: Write failing preparation and chunk tests**

Test preparation with a filtered VCF, transcript annotations, `chromosomes=["chr1"]`, VCF index provenance `"supplied"`, and transcript index provenance `"generated"`. Implement chunk construction with this exact logic:

```python
def build_chromosome_chunks(
    chromosome: str, chromosome_length: int, chunk_size_bp: int
) -> tuple[QueryChunk, ...]:
    if isinstance(chunk_size_bp, bool) or chunk_size_bp < 1:
        raise ValueError("chunk_size_bp must be a positive integer")
    if chromosome_length < 1:
        raise ValueError("VCF contig length must be positive")
    return tuple(
        QueryChunk(chromosome, start, min(chromosome_length, start + chunk_size_bp - 1))
        for start in range(1, chromosome_length + 1, chunk_size_bp)
    )
```

The preparation function signature is `prepare_carrier_inputs(vcf_path: Path, annotation_path: Path, chromosomes: Sequence[str], vcf_index_provenance: str, annotation_index_provenance: str, schema_output: Path, qc_output: Path) -> None`. Assert a 1,000-base contig with a chunk size of 400 gives regions `1-400`, `401-800`, and `801-1000`. Reject Boolean, zero, and negative chunk sizes. Verify the preparation QC records the VCF header, transcript schema, exact input paths, selected chromosomes, and both index-provenance values.

- [ ] **Step 3: Run focused tests and verify failure**

```bash
python3 -m pytest -q tests/test_carrier_extraction.py -k 'prepare or chunks'
```

Expected: FAIL because the module and functions do not exist.

- [ ] **Step 4: Implement preparation and chunk construction**

Read the transcript header with `open_text()`, resolve it with `TranscriptCarrierSchema.from_header()`, and write the schema JSON. Read VCF samples and contig lengths from the header. Reject missing requested contigs, duplicate chromosomes, and index provenance outside `{"supplied", "generated"}`. Write preparation QC with both paths, both headers, selected chromosomes, contig lengths, sample count, schema, `revel_available: true`, `vcf_index_provenance`, and `transcript_index_provenance`. Use the exact `build_chromosome_chunks()` implementation from Step 2.

- [ ] **Step 5: Write the failing exact gene-match extraction test**

Use the indexed fixture and call:

```python
extract_chromosome_carriers(
    prepared_fixture.vcf_gz,
    prepared_fixture.vat_bgz,
    prepared_fixture.carrier_schema,
    "chr1",
    250,
    tmp_path / "chr1.audit.tsv.gz",
    tmp_path / "chr1.qc.json",
)
```

Read the gzip audit with `csv.DictReader`. Assert:

- `chr1:100:A:C` yields gene-matched rows for both `ENSG000001` and `ENSG000099` without a feature file.
- `ENSG000001` selects `stop_gained`, `HC`, maximum REVEL `0.81`, and classes `lof_hc,lof_hc_or_lc`.
- `ENSG000099` remains a separate gene match.
- The homozygous `chr1:150:G:T` carrier has `sample_alt_allele_count=2` and class `missense`.
- Synonymous and stop-lost audit rows remain present with an empty initial class when applicable.
- No output column contains distance.

- [ ] **Step 6: Run the exact-match test and verify failure**

```bash
python3 -m pytest -q tests/test_carrier_extraction.py -k exact_gene_match
```

Expected: FAIL because `extract_chromosome_carriers()` is absent.

- [ ] **Step 7: Implement chromosome extraction**

Implement `extract_chromosome_carriers(vcf_path: Path, annotation_path: Path, schema_path: Path, chromosome: str, chunk_size_bp: int, audit_output: Path, qc_output: Path) -> None`.

For each non-overlapping chromosome chunk:

1. Load transcript rows from `tabix annotations chromosome:start-end` into `CarrierAnnotationChunkStore`.
2. Stream VCF rows from `tabix vcf chromosome:start-end`.
3. Parse ALT-specific carrier dosage, AC, and AF.
4. Query every gene annotation for the exact allele.
5. Emit one gzip audit row for every carrier and gene annotation.
6. Update only aggregate QC counters and progress logs.

Use the exact audit header from the specification. Format missing values as empty fields. Use `chrom:pos:ref:alt` for `variant_id`. Join class and consequence arrays with commas in deterministic order.

- [ ] **Step 8: Add failure and empty-output cases**

Test unmatched VAT alleles, missing annotation rows, missing REVEL, invalid schema, missing contig length, header-only output, chunk boundaries, multiallelic VCF records, partial genotypes, and tabix subprocess failures. Assert QC contains no sample, gene, transcript, or variant IDs.

- [ ] **Step 9: Run extraction and legacy regression tests**

```bash
python3 -m pytest -q tests/test_carrier_extraction.py tests/test_variants.py tests/test_chromosome_classification.py
```

Expected: PASS.

- [ ] **Step 10: Commit chromosome extraction**

```bash
git add src/rare_variant_enrichment/carrier_extraction.py tests/test_carrier_extraction.py tests/fixtures/transcript_annotations.tsv tests/conftest.py
git commit -m "feat: extract gene-matched variant carriers"
```

---

### Task 5: Gather audits and build the aggregated carrier table

**Files:**
- Create: `src/rare_variant_enrichment/carrier_aggregation.py`
- Create: `tests/test_carrier_aggregation.py`

**Interfaces:**
- Consumes: Per-chromosome audit gzip files and QC JSON files from Task 4.
- Produces: `gather_variant_carriers()` with deterministic final audit, aggregated carrier table, and reconciled gathered QC.

- [ ] **Step 1: Write failing deterministic gather tests**

Create small gzip audit shards in reverse chromosome order. Include one identical duplicate audit record across shards and multiple classes for one record. Call:

```python
gather_variant_carriers(
    [chr2_audit, chr1_audit],
    [chr2_qc, chr1_qc],
    preparation_qc,
    tmp_path / "variant_carrier_audit.tsv.gz",
    tmp_path / "variant_carriers.tsv.gz",
    tmp_path / "variant_carriers.qc.json",
)
```

Assert the final audit header and sort order. Assert the carrier table contains:

```text
sample_id  gene_id  gene_symbol  variant_class  n_variants  variant_ids
S1         ENSG1    GENE1        lof_hc         1           chr1:100:A:C
S1         ENSG1    GENE1        lof_hc_or_lc   1           chr1:100:A:C
```

Add a group with two distinct variants and verify `n_variants=2` plus lexically ordered IDs.

- [ ] **Step 2: Run gather tests and verify import failure**

```bash
python3 -m pytest -q tests/test_carrier_aggregation.py
```

Expected: collection FAILS because the module does not exist.

- [ ] **Step 3: Implement disk-backed gather and aggregation**

Implement `gather_variant_carriers(audit_inputs: Sequence[Path], qc_inputs: Sequence[Path], preparation_qc_path: Path, audit_output: Path, carrier_output: Path, qc_output: Path) -> None`.

Use a temporary SQLite database with an audit primary key of `(sample_id, gene_id, chrom, pos, ref, alt)`. An identical duplicate is ignored and counted. A conflicting duplicate fails. Maintain one non-empty symbol per normalized gene ID and fail on conflicts. Expand `variant_classes` into a second table keyed by sample, gene, class, and variant ID. Iterate ordered SQL results one group at a time to write both gzip outputs without holding the cohort in memory.

- [ ] **Step 4: Add QC reconciliation and validation tests**

Cover mismatched audit/QC input counts, duplicate chromosome QC, malformed JSON, duplicate JSON keys, invalid audit headers, inconsistent symbols, invalid dosage, invalid numeric fields, invalid classes, header-only shards, and a full run with no VAT joins.

Assert gathered QC embeds the preparation QC, includes per-chromosome payloads, totals, unique sample/gene/allele/pair counts, class carrier counts, audit row count, carrier row count, input schemas and paths, both index-provenance values, and the statement `quality_or_frequency_filters_applied: false`.

- [ ] **Step 5: Run gather tests**

```bash
python3 -m pytest -q tests/test_carrier_aggregation.py
```

Expected: PASS.

- [ ] **Step 6: Commit carrier gather**

```bash
git add src/rare_variant_enrichment/carrier_aggregation.py tests/test_carrier_aggregation.py
git commit -m "feat: gather and aggregate variant carriers"
```

---

### Task 6: Expose the extraction stages through the CLI

**Files:**
- Modify: `src/rare_variant_enrichment/cli.py:1-285`
- Modify: `tests/test_cli.py:1-620`

**Interfaces:**
- Consumes: `prepare_carrier_inputs()`, `extract_chromosome_carriers()`, and `gather_variant_carriers()` from Tasks 4 and 5.
- Produces: `prepare-carrier-inputs`, `extract-gene-carriers`, and `gather-gene-carriers` command contracts for WDL tasks.

- [ ] **Step 1: Write failing CLI help and dispatch tests**

Extend the command-list test with all three command names. Add monkeypatched dispatch tests with these exact arguments:

```text
prepare-carrier-inputs
  --vcf variants.vcf.gz
  --annotations transcript.tsv.bgz
  --chromosomes chr1,chr2
  --vcf-index-provenance supplied
  --transcript-index-provenance generated
  --schema-output transcript.schema.json
  --qc-output transcript.prepare.qc.json

extract-gene-carriers
  --vcf variants.vcf.gz
  --annotations transcript.tsv.bgz
  --schema transcript.schema.json
  --chromosome chr1
  --chunk-size-bp 10000000
  --audit-output chr1.audit.tsv.gz
  --qc-output chr1.qc.json

gather-gene-carriers
  --audit-input chr1.audit.tsv.gz
  --qc-input chr1.qc.json
  --preparation-qc transcript.prepare.qc.json
  --audit-output variant_carrier_audit.tsv.gz
  --carrier-output variant_carriers.tsv.gz
  --qc-output variant_carriers.qc.json
```

- [ ] **Step 2: Run CLI tests and verify failure**

```bash
python3 -m pytest -q tests/test_cli.py -k carrier
```

Expected: FAIL because the subcommands are not registered.

- [ ] **Step 3: Add parsers and dispatch branches**

Import the three functions. Add the command names to `COMMANDS`. Use `Path`, `parse_csv_strings`, positive integer validation for chunk size, and `action="append"` for gather inputs. Keep the positional function call order identical to the signatures in Tasks 4 and 5.

- [ ] **Step 4: Run CLI and full parser tests**

```bash
python3 -m pytest -q tests/test_cli.py
rare-variant-enrichment --help
```

Expected: PASS, and help lists all three new commands.

- [ ] **Step 5: Commit CLI wiring**

```bash
git add src/rare_variant_enrichment/cli.py tests/test_cli.py
git commit -m "feat: expose gene carrier extraction commands"
```

---

### Task 7: Add the standalone chromosome-scattered WDL

**Files:**
- Create: `workflows/extract_variant_carriers.wdl`
- Create: `tests/test_carrier_wdl_contract.py`
- Create: `examples/extract_variant_carriers.inputs.json`
- Modify: `.dockstore.yml:1-4`

**Interfaces:**
- Consumes: The Task 6 CLI commands and indexed VCF/transcript inputs.
- Produces: WDL workflow `ExtractVariantCarriers` with final audit, carrier table, gathered QC, per-chromosome QC, schema, and generated-or-validated annotation index.

- [ ] **Step 1: Write a failing miniwdl interface contract**

Load `workflows/extract_variant_carriers.wdl` with the same miniwdl AST helper pattern as `tests/test_wdl_contract.py`. Assert these workflow inputs:

```text
File filtered_vcf
File filtered_vcf_tbi
File transcript_annotations
File? transcript_annotations_tbi
Array[String] chromosomes = [chr1 through chr22]
Int annotation_chunk_size_bp = 10000000
String docker_image = ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main
Int prepare_cpu = 2
Int prepare_memory_gb = 8
Int prepare_disk_gb = 50
Int scatter_cpu = 2
Int scatter_memory_gb = 16
Int scatter_disk_gb = 50
Int gather_cpu = 2
Int gather_memory_gb = 32
Int gather_disk_gb = 100
Int max_retries = 1
Int scatter_preemptible = 2
```

Assert tasks `PrepareCarrierInputs`, `ExtractChromosomeCarriers`, and `GatherVariantCarriers`; one chromosome scatter; no enrichment calls or phenotype/GTF/distance inputs; and these outputs:

```text
File variant_carrier_audit_tsv_gz
File variant_carriers_tsv_gz
File variant_carriers_qc_json
Array[File] chromosome_qc_jsons
File transcript_schema_json
File generated_or_validated_transcript_annotations_tbi
String transcript_index_provenance
```

- [ ] **Step 2: Run the WDL contract and verify missing-file failure**

```bash
python3 -m pytest -q tests/test_carrier_wdl_contract.py
```

Expected: FAIL because the workflow file does not exist.

- [ ] **Step 3: Implement `PrepareCarrierInputs` with explicit logging**

The task must:

1. Log start and localized input sizes.
2. Symlink the transcript BGZ to `transcript_annotations.tsv.bgz`.
3. Symlink a supplied index or generate one with `tabix -f -S 1 -s 1 -b 2 -e 2`.
4. Symlink the VCF and its supplied index.
5. Compare requested contigs with `tabix -l` for both inputs.
6. Run `prepare-carrier-inputs` with VCF provenance `supplied` and the resolved transcript-index provenance.
7. Log index provenance, selected contig count, and completion.

Return the original transcript file, the supplied or generated index, schema JSON, preparation QC, and provenance string.

- [ ] **Step 4: Implement the scattered extraction task with explicit logging**

Localize each indexed input under matching basenames. Run `extract-gene-carriers` for one chromosome. Log the chromosome, chunk size, task start, and output counts read from QC. Do not print audit rows.

- [ ] **Step 5: Implement gather and workflow wiring**

Build repeated CLI arguments from `Array[File]` manifests, pass `PrepareCarrierInputs.preparation_qc_json` to `--preparation-qc`, run `gather-gene-carriers`, and log shard count plus final row counts. Add dynamic disk floors from localized input sizes. Do not call or import the active enrichment workflow.

- [ ] **Step 6: Run WDL syntax and contract checks**

```bash
miniwdl check workflows/extract_variant_carriers.wdl
python3 -m pytest -q tests/test_carrier_wdl_contract.py
```

Expected: PASS.

- [ ] **Step 7: Add example inputs and Dockstore registration**

Create an example JSON with all required file inputs, one chromosome for a small test, the default chunk size, and the published image. Add a second WDL entry to `.dockstore.yml` with primary descriptor `/workflows/extract_variant_carriers.wdl`.

- [ ] **Step 8: Commit the standalone WDL**

```bash
git add workflows/extract_variant_carriers.wdl tests/test_carrier_wdl_contract.py examples/extract_variant_carriers.inputs.json .dockstore.yml
git commit -m "feat: add standalone variant carrier WDL"
```

---

### Task 8: Add end-to-end extraction and GitHub Actions smoke coverage

**Files:**
- Create: `tests/test_carrier_end_to_end.py`
- Create: `tests/test_carrier_wdl_runtime.py`
- Modify: `.github/workflows/python-tests.yml:1-45`
- Modify: `README.md:1-145`

**Interfaces:**
- Consumes: All Python and WDL interfaces from Tasks 1 through 7 and the CI-built `rare-variant-enrichment:test` image.
- Produces: Hand-checked direct-Python and Docker-backed WDL smoke tests, automated only in GitHub Actions for the Docker path.

- [ ] **Step 1: Write the direct end-to-end test**

Run preparation, one-chromosome extraction, and gather on `prepared_fixture`. Assert exact audit rows and aggregated rows for:

- HC stop-gained with maximum REVEL.
- A separate neighboring VAT gene for the same allele.
- Homozygous missense dosage.
- An unclassified synonymous or stop-lost audit row.
- Correct `n_variants` and deterministic variant IDs.
- QC reconciliation between chromosome audit rows, gathered audit rows, and class counts.

- [ ] **Step 2: Run the direct end-to-end test**

```bash
python3 -m pytest -q tests/test_carrier_end_to_end.py
```

Expected: PASS.

- [ ] **Step 3: Write the Docker-backed WDL smoke test**

Follow `_require_wdl_runtime()` from `tests/test_wdl_runtime.py`, but target `workflows/extract_variant_carriers.wdl`. Supply the indexed fixture, `chromosomes=["chr1"]`, `annotation_chunk_size_bp=250`, minimal runtime values, and `docker_image=RARE_VARIANT_ENRICHMENT_TEST_IMAGE`.

Assert all seven public outputs, gzip magic bytes, exact carrier classes, `transcript_index_provenance`, and gathered QC reconciliation. Set the same required-runtime behavior so GitHub Actions fails instead of skipping.

- [ ] **Step 4: Run the WDL smoke locally only when prerequisites already exist**

```bash
python3 -m pytest -q tests/test_carrier_wdl_runtime.py
```

Expected outside CI: SKIP if Docker, miniwdl, or the prebuilt test image is absent. Do not build a local image only to remove the skip.

- [ ] **Step 5: Make GitHub Actions validate both WDLs**

Replace the single-file validation command in `python-tests.yml` with:

```bash
for wdl in workflows/*.wdl; do
  echo "Checking ${wdl}"
  miniwdl check "${wdl}"
done
```

Keep the existing GitHub Actions image build. The CI-built image supplies the runtime for the required WDL smoke tests. The separate `wdl-validation.yml` workflow already checks the `workflows/**/*.wdl` glob and needs no source change.

- [ ] **Step 6: Document standalone extraction**

Add README sections for inputs, exact allele and gene matching, transcript collapse, class rules, both table schemas, QC, example launch, the statement that inputs are prefiltered, and the statement that enrichment runs separately.

- [ ] **Step 7: Run all locally available non-Docker verification**

```bash
python3 -m pytest -q \
  tests/test_carrier_annotations.py \
  tests/test_carrier_annotation_storage.py \
  tests/test_carrier_extraction.py \
  tests/test_carrier_aggregation.py \
  tests/test_carrier_end_to_end.py \
  tests/test_carrier_wdl_contract.py
miniwdl check workflows/extract_variant_carriers.wdl
```

Expected: PASS without a local image build.

- [ ] **Step 8: Commit smoke tests and documentation**

```bash
git add tests/test_carrier_end_to_end.py tests/test_carrier_wdl_runtime.py .github/workflows/python-tests.yml README.md
git commit -m "test: smoke test variant carrier extraction"
```

---

### Task 9: Run final cross-feature verification

**Files:**
- Verify only; no planned source changes.

**Interfaces:**
- Consumes: The complete carrier extraction implementation.
- Produces: Evidence that the new workflow passes and the active LoF enrichment path is unchanged.

- [ ] **Step 1: Compile Python sources and tests**

```bash
python3 -m compileall -q src tests
```

Expected: exit code 0.

- [ ] **Step 2: Run all non-Docker Python tests**

```bash
python3 -m pytest -q
```

Expected: PASS, with Docker-backed runtime tests skipped only when prerequisites are unavailable.

- [ ] **Step 3: Validate every WDL**

```bash
for wdl in workflows/*.wdl; do miniwdl check "${wdl}"; done
```

Expected: both workflow checks PASS.

- [ ] **Step 4: Verify the existing LoF WDL contract explicitly**

```bash
python3 -m pytest -q tests/test_wdl_contract.py tests/test_lof_pc_fixture_end_to_end.py
```

Expected: PASS with no active LoF workflow interface changes.

- [ ] **Step 5: Check formatting and unintended changes**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors. Only intended implementation files differ from the execution base, plus any pre-existing user changes that were preserved.

- [ ] **Step 6: Confirm GitHub Actions smoke status after push**

```bash
gh run list --workflow "Python tests" --limit 5
```

Expected: the run for the implementation branch completes successfully, including the CI-built Docker-backed WDL smoke test.
