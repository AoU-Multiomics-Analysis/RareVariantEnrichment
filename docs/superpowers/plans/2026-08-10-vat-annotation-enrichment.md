# VAT Annotation Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add bounded-memory, gene-matched VAT consequence and LoFTEE enrichment to the chromosome-scattered rare-variant outlier workflow.

**Architecture:** A preparation task validates or creates the generic VAT tabix index and records the resolved schema. Each chromosome task divides the union of maximum-distance TSS windows into non-overlapping 10 Mb chunks, aggregates that chunk's transcript rows in temporary SQLite, streams the matching VCF records, and upserts annotation-aware carrier minima. The gather and statistics stages add annotation family/class dimensions while retaining pooled Fisher testing and global BH correction.

**Tech Stack:** WDL 1.0, Python 3.12+, standard-library `sqlite3`, `tabix`/`bgzip`, miniwdl, pytest, Docker.

## Global Constraints

- Join VCF and VAT alleles exactly on `(chromosome, position, reference, alternate)` in GRCh38; do not silently trim, left-normalize, or change chromosome labels.
- Strip only a terminal numeric Ensembl version suffix from BED and VAT gene IDs before gene matching; preserve the original BED feature ID in outputs.
- Collapse consequences to one most-severe Ensembl term per variant-gene pair; LoFTEE is a separate HC-over-LC family.
- Normalize every valid `gvs_max_af` with `min(AF, 1 - AF)` and retain values at the inclusive default `maximum_gvs_maf <= 0.01`.
- Exclude VAT-absent, missing-frequency, malformed-frequency, inconsistent-frequency, and above-threshold alleles from every family, including baseline.
- Keep Python heap bounded with streamed tabix output and temporary SQLite; do not create chromosome-sized dictionaries of VAT rows, variants, carriers, or phenotype observations.
- Query only the largest TSS window, split into non-overlapping chunks of default size `10000000`; derive smaller distance thresholds from `minimum_distance_bp`.
- Preserve one annotation-aware carrier key per `(sample_id, feature_id, ac_class, annotation_family, annotation_class)` using minimum-distance upserts.
- Publish sample-level carrier records only when `publish_carrier_audit=true`; ordinary QC must contain aggregate counts, not variant or sample identifiers.
- Keep Fisher p-values and BH FDR documented as pooled exploratory screening statistics.
- Use only the Python standard library at runtime; do not add pandas, Hail, cyvcf2, pysam, or an ORM.

---

## File structure

### New files

- `src/rare_variant_enrichment/annotations.py`: VAT schema aliases, gene normalization, consequence parsing/severity, LoFTEE collapse, MAF parsing, and configured annotation classes.
- `src/rare_variant_enrichment/annotation_storage.py`: temporary SQLite aggregation and lookup for one VAT genomic chunk.
- `src/rare_variant_enrichment/vat.py`: VAT header inspection plus supplied/generated generic tabix-index preparation.
- `tests/test_annotations.py`: pure annotation and schema tests.
- `tests/test_annotation_storage.py`: SQLite transcript-collapse and frequency-consistency tests.
- `tests/test_vat.py`: VAT preparation and index tests.
- `tests/test_storage.py`: annotation-aware carrier-store key and ordering tests.
- `tests/fixtures/variant_annotations.tsv`: transcript-granular synthetic VAT source.

### Modified files

- `src/rare_variant_enrichment/variants.py`: non-overlapping query chunks and VAT-aware chromosome classification.
- `src/rare_variant_enrichment/storage.py`: annotation-aware carrier primary key and iteration.
- `src/rare_variant_enrichment/aggregation.py`: six-column carrier parsing, reduction, and zero-VAT-match guard.
- `src/rare_variant_enrichment/statistics.py`: annotation-stratified counts, rows, global BH FDR, and VAT provenance.
- `src/rare_variant_enrichment/cli.py`: `prepare-vat` plus VAT/class/configuration arguments for classification and calculation.
- `src/rare_variant_enrichment/__init__.py`: package version `0.3.0`.
- `workflows/rare_variant_enrichment.wdl`: VAT preparation task, inputs, scatter wiring, calculation wiring, outputs, and workflow version.
- `tests/conftest.py`: bgzipped/indexed VAT fixture paths.
- `tests/test_variants.py`, `tests/test_chromosome_classification.py`: chunk and annotation-aware classification tests.
- `tests/test_aggregation.py`: annotation-aware carrier parsing and gather tests.
- `tests/test_statistics.py`, `tests/test_end_to_end.py`, `tests/test_scale_regression.py`: annotation dimensions, hand-checked results, and heap ceilings.
- `tests/test_cli.py`, `tests/test_wdl_contract.py`, `tests/test_wdl_runtime.py`: public interface, safe command rendering, and Docker-backed execution.
- `examples/rare_variant_enrichment.inputs.json`: required VAT path and new optional settings.
- `README.md`: VAT schema, filtering, output columns, QC, chunking, and revised limitations.

---

### Task 1: Define VAT schema and annotation semantics

**Files:**
- Create: `src/rare_variant_enrichment/annotations.py`
- Create: `tests/test_annotations.py`

**Interfaces:**
- Produces: `VatSchema.from_header(header: Sequence[str]) -> VatSchema`
- Produces: `VatSchema.write_json(path: Path) -> None` and `VatSchema.read_json(path: Path) -> VatSchema`
- Produces: `DEFAULT_CONSEQUENCE_CLASSES: tuple[str, ...]`.
- Produces: `normalize_gene_id(value: str) -> str`
- Produces: `parse_consequence_terms(value: str) -> tuple[str, ...]`
- Produces: `most_severe_consequence(terms: Iterable[str]) -> tuple[str | None, tuple[str, ...]]`
- Produces: `collapse_loftee(values: Iterable[str]) -> str | None`
- Produces: `parse_gvs_max_af(value: str) -> FrequencyValue`
- Produces: `build_annotation_classes(consequence_classes: Sequence[str], loftee_enabled: bool) -> list[AnnotationClass]`
- Consumed by: Tasks 2, 3, 4, 5, and 6.

- [ ] **Step 1: Write schema, gene-ID, consequence, LoFTEE, and MAF tests**

```python
from rare_variant_enrichment.annotations import (
    AnnotationClass,
    VatSchema,
    build_annotation_classes,
    collapse_loftee,
    most_severe_consequence,
    normalize_gene_id,
    parse_consequence_terms,
    parse_gvs_max_af,
)


def test_vat_schema_accepts_official_and_susie_coordinate_names():
    official = VatSchema.from_header([
        "contig", "position", "ref_allele", "alt_allele",
        "gene_id", "consequence", "gvs_max_af", "LoF",
    ])
    susie = VatSchema.from_header([
        "chrom", "pos", "ref", "alt",
        "gene_id", "consequence", "gvs_max_af",
    ])
    assert (official.chromosome, official.position, official.ref, official.alt) == (0, 1, 2, 3)
    assert official.lof == 7
    assert (susie.chromosome, susie.position, susie.ref, susie.alt) == (0, 1, 2, 3)
    assert susie.lof is None


def test_gene_ids_strip_only_terminal_numeric_versions():
    assert normalize_gene_id("ENSG00000123456.17") == "ENSG00000123456"
    assert normalize_gene_id("ENSG00000123456") == "ENSG00000123456"
    assert normalize_gene_id("GENE.A") == "GENE.A"


def test_consequence_arrays_collapse_by_ensembl_severity():
    terms = (
        *parse_consequence_terms("['intron_variant', 'splice_region_variant']"),
        *parse_consequence_terms("missense_variant&stop_gained"),
    )
    selected, unknown = most_severe_consequence(terms)
    assert selected == "stop_gained"
    assert unknown == ()


def test_loftee_and_frequency_normalization_are_deterministic():
    assert collapse_loftee(["LC", "hc", "."]) == "HC"
    converted = parse_gvs_max_af("0.999")
    assert (converted.status, converted.maf, converted.converted) == ("valid", 0.001, True)
    assert parse_gvs_max_af("NA").status == "missing"
    assert parse_gvs_max_af("not-a-number").status == "non_numeric"
    assert parse_gvs_max_af("1.1").status == "out_of_range"


def test_annotation_classes_have_stable_family_order():
    classes = build_annotation_classes(["frameshift_variant", "stop_gained"], True)
    assert classes == [
        AnnotationClass("baseline", "all_rare_variants"),
        AnnotationClass("consequence", "frameshift_variant"),
        AnnotationClass("consequence", "stop_gained"),
        AnnotationClass("loftee", "HC"),
        AnnotationClass("loftee", "LC"),
    ]
```

Also assert that missing/ambiguous aliases, duplicate consequence classes, an empty normalized gene ID, and Boolean/non-finite AF inputs raise exact `ValueError` messages. Assert `build_annotation_classes([], False)` returns the baseline class only.

- [ ] **Step 2: Run the focused tests and confirm they fail**

Run: `.venv/bin/pytest tests/test_annotations.py -q`

Expected: collection fails with `ModuleNotFoundError: No module named 'rare_variant_enrichment.annotations'`.

- [ ] **Step 3: Implement immutable schema and annotation value types**

```python
@dataclass(frozen=True)
class VatSchema:
    header: tuple[str, ...]
    chromosome: int
    position: int
    ref: int
    alt: int
    gene_id: int
    consequence: int
    gvs_max_af: int
    lof: int | None


@dataclass(frozen=True)
class FrequencyValue:
    status: Literal["valid", "missing", "non_numeric", "out_of_range"]
    maf: float | None
    converted: bool


@dataclass(frozen=True, order=True)
class AnnotationClass:
    family: Literal["baseline", "consequence", "loftee"]
    label: str


@dataclass(frozen=True)
class VariantKey:
    chromosome: str
    position: int
    ref: str
    alt: str


@dataclass(frozen=True)
class GeneAnnotation:
    consequence: str | None
    loftee: str | None
```

Resolve aliases from the exact sets `("contig", "chrom")`, `("position", "pos")`, `("ref_allele", "ref")`, and `("alt_allele", "alt")`. Require exact `gene_id`, `consequence`, and `gvs_max_af`; detect exact `LoF` optionally. Serialize indices, header, and `lof_enabled` through JSON with duplicate-key rejection on read.

- [ ] **Step 4: Implement the complete severity order and parsers**

Use this exact ordered tuple, matching current Ensembl severity from most to least severe:

```python
ENSEMBL_CONSEQUENCE_ORDER = (
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
    "sequence_variant",
)
```

Parse bracketed JSON/Python string arrays with `ast.literal_eval`, scalar `&` or comma-delimited values, and the missing tokens `""`, `"."`, `"NA"`, `"NaN"`, `"null"`, and `"[]"`. Deduplicate terms without reordering. Return unknown nonempty terms separately from `most_severe_consequence` so QC can count them.

Normalize LoFTEE by uppercasing stripped values; return HC if present, otherwise LC, otherwise `None`. Parse finite AF values in `[0, 1]`, calculate `min(value, 1.0 - value)`, and mark `converted=True` only for raw values above 0.5.

- [ ] **Step 5: Run annotation tests**

Run: `.venv/bin/pytest tests/test_annotations.py -q`

Expected: all tests pass.

- [ ] **Step 6: Commit the annotation semantic core**

```bash
git add src/rare_variant_enrichment/annotations.py tests/test_annotations.py
git commit -m "feat: define VAT annotation semantics"
```

---

### Task 2: Prepare and index VAT inputs and build non-overlapping query chunks

**Files:**
- Create: `src/rare_variant_enrichment/vat.py`
- Create: `tests/test_vat.py`
- Modify: `src/rare_variant_enrichment/variants.py:34-79`
- Modify: `src/rare_variant_enrichment/cli.py:10-110`
- Modify: `tests/test_variants.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Consumes: `VatSchema` from Task 1.
- Produces: `prepare_vat(vat_path: Path, chromosomes: Sequence[str], schema_output: Path) -> VatSchema`.
- Produces: CLI `prepare-vat --vat --chromosomes --schema-output --loftee-enabled-output`.
- Produces: `QueryChunk(chromosome: str, start: int, end: int)` with one-based inclusive coordinates.
- Produces: `build_query_chunks(features: Sequence[FeatureTss], max_distance: int, chunk_size_bp: int) -> list[QueryChunk]`.
- Consumed by: Task 5 and the WDL in Task 7.

- [ ] **Step 1: Write failing VAT preparation and chunk-boundary tests**

```python
requires_htslib = pytest.mark.skipif(
    shutil.which("tabix") is None or shutil.which("bgzip") is None,
    reason="htslib executables are required",
)


def bgzip_fixture(tmp_path: Path, source_name: str) -> Path:
    source = Path(source_name)
    output = tmp_path / f"{source.name}.bgz"
    with output.open("wb") as handle:
        subprocess.run(["bgzip", "-c", str(source)], stdout=handle, check=True)
    return output


def test_query_chunks_split_merged_windows_without_overlap():
    features = [
        FeatureTss("chr1", 10, "ENSG1.1"),
        FeatureTss("chr1", 25, "ENSG2.2"),
    ]
    chunks = build_query_chunks(features, max_distance=10, chunk_size_bp=10)
    assert chunks == [
        QueryChunk("chr1", 1, 10),
        QueryChunk("chr1", 11, 20),
        QueryChunk("chr1", 21, 30),
        QueryChunk("chr1", 31, 35),
    ]
    assert all(left.end + 1 == right.start for left, right in zip(chunks, chunks[1:]))


@requires_htslib
def test_prepare_vat_generates_and_validates_generic_tabix_index(tmp_path):
    vat = bgzip_fixture(tmp_path, "tests/fixtures/variant_annotations.tsv")
    schema_path = tmp_path / "vat_schema.json"
    schema = prepare_vat(vat, ["chr1"], schema_path)
    assert Path(f"{vat}.tbi").is_file()
    assert schema.lof is not None
    assert VatSchema.read_json(schema_path) == schema
```

Add tests that reject `chunk_size_bp <= 0`, cross-chromosome input to one chunk build, absent requested VAT contigs, missing bgzip support, and a supplied sidecar index that cannot answer `tabix -l`.

- [ ] **Step 2: Run focused tests and confirm the new names are missing**

Run: `.venv/bin/pytest tests/test_vat.py tests/test_variants.py tests/test_cli.py -q`

Expected: failures import `prepare_vat`, `QueryChunk`, or `build_query_chunks`.

- [ ] **Step 3: Implement VAT preparation**

```python
def prepare_vat(
    vat_path: Path,
    chromosomes: Sequence[str],
    schema_output: Path,
) -> VatSchema:
    with open_text(vat_path) as handle:
        header_line = next((line for line in handle if line.strip()), None)
    if header_line is None:
        raise ValueError("VAT TSV is empty")
    schema = VatSchema.from_header(header_line.rstrip("\r\n").split("\t"))
    index_path = Path(f"{vat_path}.tbi")
    if not index_path.is_file():
        subprocess.run([
            "tabix", "-f", "-S", "1",
            "-s", str(schema.chromosome + 1),
            "-b", str(schema.position + 1),
            "-e", str(schema.position + 1),
            str(vat_path),
        ], check=True)
    contigs = set(subprocess.run(
        ["tabix", "-l", str(vat_path)],
        text=True, capture_output=True, check=True,
    ).stdout.splitlines())
    missing = [chromosome for chromosome in chromosomes if chromosome not in contigs]
    if missing:
        raise ValueError("Requested chromosomes are absent from VAT index: " + ", ".join(missing))
    schema.write_json(schema_output)
    return schema
```

Translate `StopIteration`, `OSError`, and `subprocess.CalledProcessError` into input-specific `ValueError` messages while retaining exception chaining.

- [ ] **Step 4: Implement one-based inclusive chunk planning**

```python
@dataclass(frozen=True, order=True)
class QueryChunk:
    chromosome: str
    start: int
    end: int

    @property
    def tabix_region(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}"


def build_query_chunks(features, max_distance, chunk_size_bp):
    _validate_max_distance(max_distance)
    if isinstance(chunk_size_bp, bool) or not isinstance(chunk_size_bp, int) or chunk_size_bp < 1:
        raise ValueError("annotation_chunk_size_bp must be a positive integer")
    windows = sorted(
        (feature.chrom, max(1, feature.tss - max_distance), feature.tss + max_distance)
        for feature in features
    )
    merged = _merge_inclusive_windows(windows)
    return [
        QueryChunk(chrom, chunk_start, min(end, chunk_start + chunk_size_bp - 1))
        for chrom, start, end in merged
        for chunk_start in range(start, end + 1, chunk_size_bp)
    ]
```

Reject feature collections spanning more than one chromosome because classification invokes this function once per chromosome. Keep `merge_tss_windows` for BED output compatibility until Task 5 switches query output to chunk BED rows.

- [ ] **Step 5: Add and dispatch the `prepare-vat` CLI command**

Add required arguments `--vat`, `--chromosomes`, `--schema-output`, and `--loftee-enabled-output`; call `prepare_vat` from `main()`, then write exactly `true\n` or `false\n` according to `schema.lof is not None`. Assert parser dispatch with a monkeypatched function so paths and chromosome arrays are passed unchanged.

- [ ] **Step 6: Run preparation, region, and CLI tests**

Run: `.venv/bin/pytest tests/test_vat.py tests/test_variants.py tests/test_cli.py -q`

Expected: all tests pass, with htslib-dependent tests skipped only when `tabix`/`bgzip` are unavailable.

- [ ] **Step 7: Commit VAT preparation and chunk planning**

```bash
git add src/rare_variant_enrichment/vat.py src/rare_variant_enrichment/variants.py src/rare_variant_enrichment/cli.py tests/test_vat.py tests/test_variants.py tests/test_cli.py
git commit -m "feat: prepare indexed VAT chunks"
```

---

### Task 3: Aggregate transcript annotations in a bounded SQLite chunk store

**Files:**
- Create: `src/rare_variant_enrichment/annotation_storage.py`
- Create: `tests/test_annotation_storage.py`

**Interfaces:**
- Consumes: `VatSchema`, `VariantKey`, `GeneAnnotation`, and Task 1 parsers.
- Produces: `VatChunkStore(directory: Path, schema: VatSchema, maximum_gvs_maf: float, configured_consequences: Sequence[str])`.
- Produces: `VatChunkStore.ingest(fields: Sequence[str]) -> None`.
- Produces: `VatChunkStore.finalize() -> dict[str, int | float | str]`.
- Produces: `VatChunkStore.qualifying_maf(key: VariantKey) -> float | None`.
- Produces: `VatChunkStore.gene_annotation(key: VariantKey, normalized_gene_id: str) -> GeneAnnotation`.
- Consumed by: Task 5.

- [ ] **Step 1: Write failing transcript-collapse and frequency tests**

```python
OFFICIAL_HEADER = [
    "contig", "position", "ref_allele", "alt_allele",
    "gene_id", "consequence", "gvs_max_af", "LoF",
]


def row(
    gene_id: str,
    consequence: str,
    gvs_max_af: str,
    lof: str,
    *,
    position: int = 100,
) -> list[str]:
    return ["chr1", str(position), "A", "C", gene_id, consequence, gvs_max_af, lof]


def test_chunk_store_collapses_transcripts_and_uses_hc_over_lc(tmp_path):
    schema = VatSchema.from_header(OFFICIAL_HEADER)
    with VatChunkStore(tmp_path, schema, 0.01, ["stop_gained", "missense_variant"]) as store:
        store.ingest(row("ENSG000001.1", "['missense_variant']", "0.001", "LC"))
        store.ingest(row("ENSG000001.2", "['stop_gained']", "0.001", "HC"))
        qc = store.finalize()
        key = VariantKey("chr1", 100, "A", "C")
        assert store.qualifying_maf(key) == 0.001
        assert store.gene_annotation(key, "ENSG000001") == GeneAnnotation("stop_gained", "HC")
        assert qc["vat_rows"] == 2
        assert qc["unique_vat_alleles"] == 1


def test_chunk_store_excludes_inconsistent_and_common_frequency(tmp_path):
    schema = VatSchema.from_header(OFFICIAL_HEADER)
    with VatChunkStore(tmp_path, schema, 0.01, ["frameshift_variant"]) as store:
        store.ingest(row("ENSG1", "['frameshift_variant']", "0.001", "HC", position=100))
        store.ingest(row("ENSG1", "['frameshift_variant']", "0.002", "HC", position=100))
        store.ingest(row("ENSG2", "['frameshift_variant']", "0.02", "HC", position=200))
        qc = store.finalize()
        assert store.qualifying_maf(VariantKey("chr1", 100, "A", "C")) is None
        assert store.qualifying_maf(VariantKey("chr1", 200, "A", "C")) is None
        assert qc["inconsistent_frequency_alleles"] == 1
        assert qc["above_maf_threshold_alleles"] == 1
```

Also test missing/non-numeric/out-of-range frequencies, AF above 0.5 conversion, exact 0.01 retention, gene mismatch lookup, unknown terms, pre-collapsed rows, malformed coordinates, and use-before-finalize rejection.

- [ ] **Step 2: Run the storage tests and confirm import failure**

Run: `.venv/bin/pytest tests/test_annotation_storage.py -q`

Expected: `ModuleNotFoundError` for `rare_variant_enrichment.annotation_storage`.

- [ ] **Step 3: Implement the temporary SQLite schema and context manager**

```sql
CREATE TABLE variant_frequency (
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    raw_af REAL,
    maf REAL,
    status TEXT NOT NULL,
    converted INTEGER NOT NULL,
    PRIMARY KEY (chromosome, position, ref, alt)
) WITHOUT ROWID;

CREATE TABLE gene_annotation (
    chromosome TEXT NOT NULL,
    position INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    gene_id TEXT NOT NULL,
    consequence_rank INTEGER,
    consequence TEXT,
    loftee_rank INTEGER,
    loftee TEXT,
    PRIMARY KEY (chromosome, position, ref, alt, gene_id)
) WITHOUT ROWID;

CREATE TABLE unknown_consequence (
    term TEXT PRIMARY KEY,
    row_count INTEGER NOT NULL
) WITHOUT ROWID;
```

Use `journal_mode=OFF`, `synchronous=OFF`, `temp_store=FILE`, and a 4 MiB negative `cache_size`, matching the existing carrier store. Delete the SQLite file in `close()` even when ingestion raises.

- [ ] **Step 4: Implement deterministic upserts and final QC**

For frequency conflicts, preserve the first valid MAF and set `status='inconsistent'` when a later valid MAF differs by more than `1e-12`; any missing or invalid frequency sets the corresponding terminal exclusion status unless the allele is already inconsistent. For gene annotations, retain the smallest Ensembl severity rank and the smallest LoFTEE rank where HC is 0 and LC is 1.

`finalize()` commits, computes unique-allele QC with SQL, records `observed_raw_gvs_max_af`, converted counts, exclusion reasons, unknown-term counts, and prevents later ingestion. `qualifying_maf()` returns a number only for `status='valid'` and `maf <= maximum_gvs_maf`.

- [ ] **Step 5: Run SQLite annotation tests**

Run: `.venv/bin/pytest tests/test_annotation_storage.py -q`

Expected: all tests pass.

- [ ] **Step 6: Commit the chunk store**

```bash
git add src/rare_variant_enrichment/annotation_storage.py tests/test_annotation_storage.py
git commit -m "feat: aggregate VAT annotations on disk"
```

---

### Task 4: Migrate carrier reduction and statistics to annotation-aware keys

**Files:**
- Modify: `src/rare_variant_enrichment/storage.py:7-120`
- Modify: `src/rare_variant_enrichment/aggregation.py:9-80`
- Modify: `src/rare_variant_enrichment/statistics.py:20-375,498-573`
- Modify: `src/rare_variant_enrichment/variants.py:188-310`
- Create: `tests/test_storage.py`
- Modify: `tests/test_aggregation.py`
- Modify: `tests/test_statistics.py`
- Modify: `tests/test_end_to_end.py`
- Modify: `tests/test_scale_regression.py`
- Modify: `tests/test_chromosome_classification.py`

**Interfaces:**
- Consumes: `AnnotationClass` and `build_annotation_classes` from Task 1.
- Produces: six-column carrier TSV `(sample_id, feature_id, ac_class, annotation_family, annotation_class, minimum_distance_bp)`.
- Produces: `MinimumDistanceStore.upsert(sample_id, feature_id, ac_class, annotation_family, annotation_class, minimum_distance_bp) -> None`.
- Produces: `MinimumDistanceStore.iter_feature(feature_id) -> Iterator[tuple[str, str, str, str, int]]`.
- Produces: the following extended calculation interface:

```python
def calculate_enrichment(
    phenotype_bed: Path,
    shared_samples_path: Path,
    carriers_path: Path,
    selected_features_path: Path,
    exact_ac: Sequence[int],
    cumulative_ac_max: Sequence[int],
    z_thresholds: Sequence[float],
    distance_thresholds: Sequence[int],
    tail: str,
    output_tsv: Path,
    output_json: Path,
    *,
    consequence_classes: Sequence[str] = (),
    loftee_enabled: bool = False,
    phenotype_qc_path: Path | None = None,
    chromosome_qc_path: Path | None = None,
    selected_chromosomes: Sequence[str] | None = None,
    container_image: str | None = None,
    workflow_version: str = "unknown",
    max_retries: int = 0,
    index_provenance: str = "unknown",
    vat_index_provenance: str = "unknown",
    maximum_gvs_maf: float = 0.01,
    annotation_chunk_size_bp: int = 10_000_000,
) -> None:
```
- Consumed by: Tasks 5, 6, and 7.

- [ ] **Step 1: Change fixtures to the six-column carrier contract and add failing dimension tests**

```python
carrier_header = (
    "sample_id\tfeature_id\tac_class\tannotation_family\t"
    "annotation_class\tminimum_distance_bp\n"
)


def write_carriers(path: Path, rows: list[tuple[str, str, str, str, str, int]]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(carrier_header)
        for values in rows:
            handle.write("\t".join(map(str, values)) + "\n")
    return path


def test_gather_deduplicates_within_but_not_across_annotation_classes(tmp_path):
    first = write_carriers(tmp_path / "a.tsv", [
        ("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 50),
        ("S1", "ENSG1.1", "AC=1", "consequence", "missense_variant", 20),
    ])
    second = write_carriers(tmp_path / "b.tsv", [
        ("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 10),
    ])
    q1 = tmp_path / "q1.json"
    q2 = tmp_path / "q2.json"
    q1.write_text('{"chromosome":"chr1"}')
    q2.write_text('{"chromosome":"chr2"}')
    output = tmp_path / "gathered.tsv"
    qc_output = tmp_path / "qc.tsv"
    gather_outputs([first, second], [q1, q2], output, qc_output)
    assert output.read_text().splitlines()[1:] == [
        "S1\tENSG1.1\tAC=1\tconsequence\tmissense_variant\t20",
        "S1\tENSG1.1\tAC=1\tconsequence\tstop_gained\t10",
    ]
```

Add a statistics test with baseline, stop-gained, and missense carrier rows. Assert three result strata, distinct hand-checked 2x2 cells, and BH values computed across all three rows together.

- [ ] **Step 2: Run carrier, aggregation, statistics, and end-to-end tests to verify contract failures**

Run: `.venv/bin/pytest tests/test_storage.py tests/test_aggregation.py tests/test_statistics.py tests/test_end_to_end.py tests/test_scale_regression.py tests/test_chromosome_classification.py -q`

Expected: failures report the old four-column schema and missing annotation arguments.

- [ ] **Step 3: Extend `MinimumDistanceStore` and carrier parsing**

Use this primary key:

```sql
PRIMARY KEY (
    feature_id, sample_id, ac_class, annotation_family, annotation_class
)
```

Validate nonempty family/class values, include them in conflict targets and stable ordering, and make `_iter_carrier_file()` yield:

```python
yield sample_id, feature_id, ac_class, annotation_family, annotation_class, distance
```

Update `gather_outputs()` to preserve independent classes. During this migration, existing VCF classification emits only `("baseline", "all_rare_variants")`; Task 5 adds VAT-specific rows.

- [ ] **Step 4: Generalize calculation arrays over configured annotation classes**

```python
annotation_classes = build_annotation_classes(consequence_classes, loftee_enabled)
annotation_indexes = {
    (annotation.family, annotation.label): index
    for index, annotation in enumerate(annotation_classes)
}
carrier_counts = [
    [[0 for _ in distances] for _ in ac_classes]
    for _ in annotation_classes
]
carrier_outlier_counts = [
    [
        [[0 for _ in distances] for _ in ac_classes]
        for _ in annotation_classes
    ]
    for _ in thresholds
]
```

Reject carrier annotation keys not present in `annotation_indexes`. Insert `annotation_family` and `annotation_class` immediately after `tail` in `ENRICHMENT_HEADER`. Loop in stable order `threshold`, annotation class, AC class, distance, then run one `benjamini_hochberg()` call over every emitted p-value.

- [ ] **Step 5: Update all baseline fixtures and bounded-memory ceilings**

Every pre-annotation carrier fixture receives `baseline\tall_rare_variants`. Existing hand-checked cells remain unchanged within that stratum. Keep gather below 4 MiB and calculation below 8 MiB for the existing key counts; print peak MiB values on failure diagnostics.

- [ ] **Step 6: Run the migrated core tests**

Run: `.venv/bin/pytest tests/test_storage.py tests/test_aggregation.py tests/test_statistics.py tests/test_end_to_end.py tests/test_scale_regression.py tests/test_chromosome_classification.py -q`

Expected: all tests pass.

- [ ] **Step 7: Commit the annotation-aware carrier/statistics model**

```bash
git add src/rare_variant_enrichment/storage.py src/rare_variant_enrichment/aggregation.py src/rare_variant_enrichment/statistics.py src/rare_variant_enrichment/variants.py tests/test_storage.py tests/test_aggregation.py tests/test_statistics.py tests/test_end_to_end.py tests/test_scale_regression.py tests/test_chromosome_classification.py
git commit -m "feat: stratify carrier enrichment by annotation"
```

---

### Task 5: Stream VAT and VCF chunks through chromosome classification

**Files:**
- Create: `tests/fixtures/variant_annotations.tsv`
- Modify: `tests/conftest.py:8-70`
- Modify: `src/rare_variant_enrichment/variants.py:149-330`
- Modify: `src/rare_variant_enrichment/cli.py:25-105`
- Modify: `tests/test_chromosome_classification.py`
- Modify: `tests/test_cli.py`
- Modify: `tests/test_end_to_end.py`

**Interfaces:**
- Consumes: `VatSchema`, `VatChunkStore`, `QueryChunk`, and six-column carrier store.
- Produces: `classify_chromosome(vcf_path, vat_path, vat_schema_path, features_path, shared_samples_path, chromosome, exact_ac, cumulative_ac_max, consequence_classes, maximum_gvs_maf, max_distance, annotation_chunk_size_bp, carrier_output, regions_output, qc_output) -> None`.
- Produces: per-chromosome aggregate VAT/VCF/join/frequency/consequence/LoFTEE QC.
- Consumed by: Task 7 WDL and Task 8 end-to-end tests.

- [ ] **Step 1: Add a transcript-granular fixture and failing classifier test**

Create `tests/fixtures/variant_annotations.tsv` with the exact header:

```text
contig	position	ref_allele	alt_allele	gene_id	consequence	gvs_max_af	LoF
```

Include duplicate transcript rows for `chr1:100 A>C` matching `ENSG000001.1`, one missense LC and one stop-gained HC; a neighboring-gene frameshift that must not label the tested gene; a common `chr1:150 G>T`; an AF-above-0.5 convertible allele; and one missing-frequency allele.

```python
@requires_htslib
def test_classification_is_gene_matched_chunked_and_annotation_aware(prepared_fixture, tmp_path):
    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc = tmp_path / "qc.json"
    classify_chromosome(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        prepared_fixture.vat_schema,
        prepared_fixture.features,
        prepared_fixture.samples,
        "chr1", [1, 2], [1, 2],
        ["stop_gained", "frameshift_variant", "missense_variant"],
        0.01, 100, 25,
        carriers, regions, qc,
    )
    rows = carriers.read_text().splitlines()
    assert "S1\tENSG000001.1\tAC=1\tbaseline\tall_rare_variants\t0" in rows
    assert "S1\tENSG000001.1\tAC=1\tconsequence\tstop_gained\t0" in rows
    assert "S1\tENSG000001.1\tAC=1\tloftee\tHC\t0" in rows
    assert not any("frameshift_variant" in row and "ENSG000001.1" in row for row in rows)
    summary = json.loads(qc.read_text())
    assert summary["annotation_chunk_count"] > 1
    assert summary["vcf_tabix_query_count"] == summary["annotation_chunk_count"]
    assert summary["vat_tabix_query_count"] == summary["annotation_chunk_count"]
```

- [ ] **Step 2: Run classifier, CLI, and end-to-end tests to verify failures**

Run: `.venv/bin/pytest tests/test_chromosome_classification.py tests/test_cli.py tests/test_end_to_end.py -q`

Expected: `classify_chromosome()` rejects the new VAT arguments.

- [ ] **Step 3: Extend the prepared fixture with indexed VAT and versioned Ensembl genes**

Add `vat_bgz`, `vat_tbi`, `vat_schema`, and `features` to `PreparedFixture`. Generate the bgzip/index with local htslib or the existing Docker fallback. Change synthetic phenotype/features identifiers to versioned `ENSG000001.1`, `ENSG000002.2`, `ENSG000003.3`, and `ENSG000004.4`; keep z-scores and VCF genotypes unchanged so existing cell counts remain hand-checkable.

- [ ] **Step 4: Replace whole-region VCF streaming with per-chunk dual streaming**

For each `QueryChunk`:

```python
with VatChunkStore(
    carrier_output.parent,
    vat_schema,
    maximum_gvs_maf,
    consequence_classes,
) as annotations:
    _stream_tabix_tsv(vat_path, chunk.tabix_region, annotations.ingest)
    chunk_qc = annotations.finalize()
    for fields in _stream_tabix_vcf(vcf_path, chunk.tabix_region):
        for allele in parse_variant_alleles(fields, sample_ids, shared_samples, qc=qc):
            key = VariantKey(allele.chrom, allele.pos, allele.ref, allele.alt)
            if annotations.qualifying_maf(key) is None:
                continue
            for feature in feature_index.nearby(allele.pos, max_distance):
                gene = normalize_gene_id(feature.feature_id)
                annotation = annotations.gene_annotation(key, gene)
                _upsert_allele_carriers(
                    minimum_distances,
                    allele,
                    feature,
                    ac_classes,
                    annotation,
                    configured_consequence_set,
                )
```

Always upsert baseline after frequency qualification. Upsert consequence only when the selected term is configured. Upsert HC/LC only when present. Continue reducing minimum distance independently for every AC and annotation key.

- [ ] **Step 5: Write chunk BED output and comprehensive aggregate QC**

Write one BED row per non-overlapping query chunk as `chrom`, `start - 1`, `end`. Sum chunk-store QC without storing unknown variant identifiers. Add `annotation_chunk_count`, `vat_tabix_query_count`, `vcf_tabix_query_count`, `vat_joined_alt_alleles`, `vat_unmatched_alt_alleles`, `gene_matched_variant_feature_pairs`, `gene_unmatched_variant_feature_pairs`, and emitted-key counts by family.

If a chromosome has no selected features, emit header-only carrier/region outputs and zero-valued annotation QC without invoking tabix. Preserve strict VCF contig validation.

- [ ] **Step 6: Wire required classify CLI arguments**

Add `--vat`, `--vat-schema`, `--consequence-classes`, `--maximum-gvs-maf`, and `--annotation-chunk-size-bp`. Reuse CSV parsing and validate finite MAF in Python domain code rather than only argparse.

- [ ] **Step 7: Run classifier and miniature Python pipeline tests**

Run: `.venv/bin/pytest tests/test_chromosome_classification.py tests/test_cli.py tests/test_end_to_end.py -q`

Expected: all tests pass, including exact chunk boundaries, gene matching, common/missing exclusion, and HC collapse.

- [ ] **Step 8: Commit chunked VAT-aware classification**

```bash
git add tests/fixtures/variant_annotations.tsv tests/conftest.py src/rare_variant_enrichment/variants.py src/rare_variant_enrichment/cli.py tests/test_chromosome_classification.py tests/test_cli.py tests/test_end_to_end.py
git commit -m "feat: classify annotated variants in chunks"
```

---

### Task 6: Complete annotation configuration, global FDR, and provenance

**Files:**
- Modify: `src/rare_variant_enrichment/statistics.py:20-375,498-573`
- Modify: `src/rare_variant_enrichment/cli.py:40-130`
- Modify: `src/rare_variant_enrichment/aggregation.py:14-38`
- Modify: `tests/test_statistics.py`
- Modify: `tests/test_aggregation.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Consumes: configured consequence classes and VAT `lof_enabled` from Tasks 1 and 2.
- Produces: final result rows across baseline, configured consequence, and optional HC/LC classes, including zero-carrier rows.
- Produces: calculation provenance fields `vat_index`, `maximum_gvs_maf`, `annotation_chunk_size_bp`, `severity_order_version`, `consequence_classes`, and `loftee_enabled`.
- Produces: gathered zero-match guard using summed chromosome QC.

- [ ] **Step 1: Write failing zero-row, global-FDR, provenance, and zero-match tests**

```python
def test_calculation_emits_configured_zero_carrier_annotation_rows(tmp_path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\n"
        "chr1\t99\t100\tENSG1.1\t3\t0\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\n")
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tENSG1.1\n")
    carriers_with_baseline_only = tmp_path / "carriers.tsv"
    carriers_with_baseline_only.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\t"
        "annotation_class\tminimum_distance_bp\n"
        "S1\tENSG1.1\tAC=1\tbaseline\tall_rare_variants\t0\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary = tmp_path / "enrichment.json"
    calculate_enrichment(
        bed, samples, carriers_with_baseline_only, features,
        [1], [], [2.0], [100], "absolute", output, summary,
        consequence_classes=["stop_gained"],
        loftee_enabled=True,
        vat_index_provenance="generated",
        maximum_gvs_maf=0.01,
        annotation_chunk_size_bp=10_000_000,
    )
    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert [(row["annotation_family"], row["annotation_class"], row["n11"]) for row in rows] == [
        ("baseline", "all_rare_variants", "1"),
        ("consequence", "stop_gained", "0"),
        ("loftee", "HC", "0"),
        ("loftee", "LC", "0"),
    ]
    assert json.loads(summary.read_text())["provenance"]["vat_index"] == "generated"
```

Add a gather test whose QC totals have `classified_alt_alleles > 0` and `vat_joined_alt_alleles == 0`; assert `ValueError("No queried VCF ALT alleles matched VAT allele keys")`. Verify an empty-feature chromosome does not trigger that guard when other chromosomes have no eligible queried ALT alleles.

- [ ] **Step 2: Run statistics, aggregation, and CLI tests to verify failures**

Run: `.venv/bin/pytest tests/test_statistics.py tests/test_aggregation.py tests/test_cli.py -q`

Expected: missing provenance parameters or zero-match guard assertions fail.

- [ ] **Step 3: Validate and record annotation calculation configuration**

Require unique configured consequences that exist in the Ensembl order, `0 <= maximum_gvs_maf <= 0.5`, positive integer chunk size, and VAT index provenance in `generated|supplied|unknown`. Include the annotation classes and configuration under `analysis_parameters` and `provenance`; include `severity_order_version = "Ensembl release 116"`.

- [ ] **Step 4: Add calculation CLI arguments and dispatch**

Add `--consequence-classes`, `--loftee-enabled` as `true|false`, `--vat-index-provenance`, `--maximum-gvs-maf`, and `--annotation-chunk-size-bp`. Materialize strings through files in WDL later; CLI tests must prove values with shell metacharacters never become executable command fragments.

- [ ] **Step 5: Add the gathered zero-match safety check**

After validating all chromosome QC records, sum `classified_alt_alleles` and `vat_joined_alt_alleles`. Raise only when the first total is positive and the second is zero. Keep all per-chromosome QC rows available for successful runs.

- [ ] **Step 6: Run focused calculation and gather tests**

Run: `.venv/bin/pytest tests/test_statistics.py tests/test_aggregation.py tests/test_cli.py -q`

Expected: all tests pass and BH FDR remains monotone across all annotation rows.

- [ ] **Step 7: Commit final calculation and provenance behavior**

```bash
git add src/rare_variant_enrichment/statistics.py src/rare_variant_enrichment/aggregation.py src/rare_variant_enrichment/cli.py tests/test_statistics.py tests/test_aggregation.py tests/test_cli.py
git commit -m "feat: report annotation enrichment provenance"
```

---

### Task 7: Wire VAT preparation and annotation dimensions through WDL

**Files:**
- Modify: `workflows/rare_variant_enrichment.wdl:1-575`
- Modify: `tests/test_wdl_contract.py:1-360`

**Interfaces:**
- Consumes: all Python CLI contracts from Tasks 2, 5, and 6.
- Produces: required workflow input `variant_annotation_table` and optional `variant_annotation_table_tbi`.
- Produces: `generated_or_validated_vat_tbi`, `vat_index_provenance`, `vat_schema_json`, and annotation-aware existing outputs.
- Consumed by: Task 8 Docker-backed runtime tests and Task 9 documentation.

- [ ] **Step 1: Extend parsed WDL contract tests before editing WDL**

Assert these exact new public inputs:

```python
"variant_annotation_table": {"type": "File", "default": None},
"variant_annotation_table_tbi": {"type": "File?", "default": None},
"maximum_gvs_maf": {"type": "Float", "default": 0.01},
"annotation_chunk_size_bp": {"type": "Int", "default": 10000000},
"consequence_classes": {
    "type": "Array[String]",
    "default": [
        "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
        "frameshift_variant", "stop_lost", "start_lost",
        "inframe_insertion", "inframe_deletion", "missense_variant",
        "protein_altering_variant", "splice_region_variant",
        "synonymous_variant", "coding_sequence_variant",
    ],
},
```

Update `miniwdl input_template` expectation to require phenotype BED, VCF, and VAT. Assert new task runtime keys, output types, scatter wiring, and safe generated-file boundaries.

- [ ] **Step 2: Run WDL contract tests and confirm failures**

Run: `.venv/bin/pytest tests/test_wdl_contract.py -q`

Expected: public input and task/output contract assertions fail.

- [ ] **Step 3: Add `PrepareVatIndex` WDL task**

The task localizes the VAT as `annotations.tsv.bgz`. If an index is supplied, symlink it to `annotations.tsv.bgz.tbi` and write `supplied`; otherwise write `generated` and let `prepare-vat` create the sidecar. Materialize chromosomes through `write_lines`, invoke:

```bash
rare-variant-enrichment prepare-vat \
  --vat annotations.tsv.bgz \
  --chromosomes "$chromosomes_csv" \
  --schema-output vat_schema.json \
  --loftee-enabled-output loftee_enabled.txt
```

Read `loftee_enabled.txt` with WDL `read_boolean`. Emit the original VAT, selected/generated index, schema, Boolean LoFTEE availability, and provenance.

- [ ] **Step 4: Extend scatter and calculation task inputs**

Localize both VCF and VAT sidecars in `ClassifyChromosome`, pass schema/configuration/chunk inputs to the classifier, and retain one chromosome scatter. Pass consequence classes, LoFTEE availability, both index provenance strings, MAF threshold, and chunk size to `CalculateEnrichment` through materialized files.

Set `workflow_version = "0.3.0"`.

- [ ] **Step 5: Add workflow outputs and provenance**

Add:

```wdl
File generated_or_validated_vat_tbi = PrepareVatIndex.vat_tbi
File vat_schema_json = PrepareVatIndex.vat_schema_json
String vat_index_provenance = PrepareVatIndex.index_provenance
```

Keep `carrier_minimum_distances_tsv` optional and keep chromosome query regions as the per-chromosome chunk BED outputs.

- [ ] **Step 6: Validate command rendering and static WDL**

Run: `.venv/bin/pytest tests/test_wdl_contract.py -q`

Run: `.venv/bin/miniwdl check workflows/rare_variant_enrichment.wdl`

Expected: all contract tests pass and miniwdl reports no errors.

- [ ] **Step 7: Commit WDL wiring**

```bash
git add workflows/rare_variant_enrichment.wdl tests/test_wdl_contract.py
git commit -m "feat: wire VAT annotations through WDL"
```

---

### Task 8: Prove supplied/generated indexes, hand-checked enrichment, and bounded memory end to end

**Files:**
- Modify: `tests/conftest.py`
- Modify: `tests/test_wdl_runtime.py:1-170`
- Modify: `tests/test_end_to_end.py:1-120`
- Modify: `tests/test_scale_regression.py:1-110`

**Interfaces:**
- Consumes: complete Python and WDL workflow from Tasks 1-7.
- Produces: executable acceptance coverage for both index modes and the approved biological semantics.

- [ ] **Step 1: Extend the Docker-backed WDL matrix before running it**

For both parametrized modes:

- supplied VCF index + supplied VAT index + carrier audit enabled;
- generated VCF index + generated VAT index + carrier audit disabled.

Pass a chr1-only workflow input, `maximum_gvs_maf=0.01`, `annotation_chunk_size_bp=25`, and three consequence classes. Assert generated/validated index files exist, provenance strings match, the audit optionality remains correct, and selected chromosomes remain `['chr1']`.

- [ ] **Step 2: Add hand-checked annotation result assertions**

Index rows by:

```python
(
    row["z_threshold"],
    row["annotation_family"],
    row["annotation_class"],
    row["ac_class"],
    row["distance_bp"],
)
```

Assert exact `n11`, `n10`, `n01`, and `n00` for baseline, stop-gained, missense, HC, and LC rows. Assert the common and missing-frequency fixture variants never create carriers. Assert FDR values equal `benjamini_hochberg()` over the complete emitted row list.

- [ ] **Step 3: Add a chunk-store scale regression**

Generate 200,000 transcript rows for one chunk across 25,000 allele-gene keys while tracing Python allocations. Stream rows into `VatChunkStore`, finalize, query representative keys, and assert peak Python heap remains below 12 MiB. Do not read the generated TSV back into a Python list.

- [ ] **Step 4: Run Python acceptance and scale tests**

Run: `.venv/bin/pytest tests/test_end_to_end.py tests/test_scale_regression.py -q`

Expected: all non-Docker tests pass within the existing CI time budget.

- [ ] **Step 5: Run required Docker-backed WDL tests**

Run: `RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME=1 .venv/bin/pytest tests/test_wdl_runtime.py -q`

Expected: both supplied/generated index modes pass; missing Docker/image prerequisites cause an explicit failure rather than a skip.

- [ ] **Step 6: Commit acceptance and scale coverage**

```bash
git add tests/conftest.py tests/test_wdl_runtime.py tests/test_end_to_end.py tests/test_scale_regression.py
git commit -m "test: cover VAT enrichment end to end"
```

---

### Task 9: Document the public interface, bump the release, and run final verification

**Files:**
- Modify: `README.md:11-123`
- Modify: `examples/rare_variant_enrichment.inputs.json`
- Modify: `src/rare_variant_enrichment/__init__.py:1`
- Modify: `pyproject.toml:7-11`

**Interfaces:**
- Consumes: final WDL input/output contract and QC names.
- Produces: user-facing run instructions for VAT annotation enrichment and package/workflow version `0.3.0`.

- [ ] **Step 1: Update the example input JSON**

Add:

```json
"RareVariantEnrichment.variant_annotation_table": "variant_annotations.tsv.bgz",
"RareVariantEnrichment.maximum_gvs_maf": 0.01,
"RareVariantEnrichment.annotation_chunk_size_bp": 10000000
```

Continue omitting both optional index inputs and `chromosomes` so the example demonstrates automatic indexing and the chr1-chr22 default.

- [ ] **Step 2: Rewrite README inputs, outputs, scale, and limitations**

Document official and SuSiE-style VAT aliases, transcript collapse, versioned Ensembl matching, exact allele joins, `gvs_max_af` conversion/filtering, optional LoF, coding consequence defaults, supplied/generated VAT index behavior, coordinate-sort requirement, chunk size, and six-column controlled carrier audit.

Update the enrichment row definition to z threshold × annotation family/class × AC class × distance. State that BH is global across every emitted row and that the baseline also requires valid VAT frequency. Replace the statement that functional/external annotations are not implemented with the remaining non-goals from the approved spec.

- [ ] **Step 3: Bump Python and workflow-facing versions**

Set both `src/rare_variant_enrichment/__init__.py` and `pyproject.toml` to `0.3.0`; WDL was set to `0.3.0` in Task 7. Update provenance assertions to require all three values agree.

- [ ] **Step 4: Run the complete verification suite**

Run: `.venv/bin/pytest -q`

Run: `.venv/bin/miniwdl check workflows/rare_variant_enrichment.wdl`

Run: `.venv/bin/python -m json.tool examples/rare_variant_enrichment.inputs.json`

Run: `git diff --check`

Expected: all tests pass; only explicitly marked Docker prerequisites may skip in the ordinary pytest run; WDL and JSON validation are clean; no whitespace errors are reported.

- [ ] **Step 5: Run the repository review checklist**

Verify manually that:

- no ordinary JSON/TSV QC contains sample IDs, variant IDs, or VAT rows;
- every tabix subprocess receives an argument list rather than a shell string;
- each chunk is non-overlapping and each smaller distance derives from minimum distance;
- HC and LC are mutually exclusive per variant-gene while consequence remains separate;
- all configured zero-carrier classes emit rows and participate in global FDR;
- supplied indexes are validated and generated indexes are workflow outputs; and
- chromosome defaults and explicit chromosome overrides remain unchanged.

- [ ] **Step 6: Commit documentation and release metadata**

```bash
git add README.md examples/rare_variant_enrichment.inputs.json src/rare_variant_enrichment/__init__.py pyproject.toml
git commit -m "docs: release VAT annotation enrichment"
```

- [ ] **Step 7: Record final evidence for review**

Run: `git status --short --branch`

Run: `git log --oneline --decorate -10`

Expected: clean worktree on `codex/vat-annotation-enrichment`, with one focused commit per task and the final verification results available for the PR summary.
