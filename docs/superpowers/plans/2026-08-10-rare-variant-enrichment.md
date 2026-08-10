# Rare Variant Enrichment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and validate a WDL workflow that measures pooled rare-variant enrichment among molecular outliers across z-score, allele-count, and TSS-distance thresholds.

**Architecture:** A WDL 1.0 workflow prepares or validates the VCF tabix index, streams the wide prepare_QTL phenotype BED, and scatters one maximum-window `tabix -R` extraction per chromosome. A dependency-light Python package parses phenotypes and variants, reduces carrier records to minimum TSS distances, gathers chromosome results, and computes pooled 2x2 enrichment statistics without expanding the full phenotype matrix into long format.

**Tech Stack:** WDL 1.0, Python 3.12 standard library, pytest, htslib `bgzip`/`tabix`, miniwdl, Docker/GHCR.

## Global Constraints

- The phenotype input is the prepare_QTL wide BED: `#chr`, `start`, `end`, feature ID, then one scaled-residual z-score column per sample.
- Every BED feature interval must be one base wide, and `TSS = end` converts the BED coordinate to the 1-based TSS used with VCF `POS`.
- Analyze only sample IDs shared by the BED and VCF; report samples unique to either input.
- The rare-variant VCF must be coordinate-sorted and bgzip-compressed. Accept an optional `.tbi`; generate and validate one when absent.
- Query each chromosome once with `tabix -h -R` over merged windows expanded to the largest distance threshold.
- Apply all smaller symmetric TSS-distance thresholds from exact `abs(VCF_POS - TSS)` values without rereading the VCF.
- Evaluate multiallelic VCF records one ALT allele at a time, using the corresponding `INFO/AC` value or genotype-derived AC.
- When `INFO/AC` is absent, derive AC from every VCF sample genotype; restrict only the emitted carrier IDs to the BED/VCF sample intersection.
- Support exact AC classes and cumulative `AC<=N` classes independently.
- Reduce each `(sample_id, feature_id, ac_class)` key to its minimum qualifying distance before enrichment testing.
- Evaluate absolute, positive, or negative outlier tails and multiple z-score thresholds.
- Emit raw 2x2 counts, carrier-rate ratio, continuity-corrected odds ratio, two-sided Fisher exact p-value, and Benjamini-Hochberg FDR.
- Treat pooled Fisher tests as screening statistics because samples and genes recur; state this limitation in output metadata and documentation.
- Use `ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main` as the default runtime image while allowing the WDL caller to override it.
- Do not add functional-annotation filters, ancestry adjustment, watershed inference, or feature-specific tests in this implementation.

## File Map

- `pyproject.toml`: package metadata, console entry point, pytest configuration, and Python version floor.
- `src/rare_variant_enrichment/io.py`: compressed text opening, phenotype BED streaming, VCF header/record parsing, and TSV/JSON writers.
- `src/rare_variant_enrichment/phenotypes.py`: BED validation, sample intersection, TSS extraction, and outlier classification.
- `src/rare_variant_enrichment/variants.py`: AC classes, per-ALT genotype carriers, merged query windows, distance matching, and chromosome classification.
- `src/rare_variant_enrichment/aggregation.py`: minimum-distance carrier gather and chromosome QC gather.
- `src/rare_variant_enrichment/statistics.py`: Fisher exact test, rate/odds ratios, BH correction, and final enrichment table generation.
- `src/rare_variant_enrichment/cli.py`: command-line subcommands used by WDL tasks.
- `workflows/rare_variant_enrichment.wdl`: prepare, scatter, gather, statistics, runtime, and workflow outputs.
- `tests/`: unit, CLI integration, WDL-contract, and miniature end-to-end tests.
- `envs/Dockerfile`: reproducible Python/tabix runtime image.
- `.github/workflows/python-tests.yml`: Python and CLI test execution.
- `.github/workflows/docker-image.yml`: image build triggers for package files and main-branch publishing.
- `.dockstore.yml`: workflow registration.
- `README.md`: input schema, statistics, execution example, output columns, and limitations.

---

### Task 1: Python package, CLI shell, and runtime image

**Files:**
- Create: `pyproject.toml`
- Create: `src/rare_variant_enrichment/__init__.py`
- Create: `src/rare_variant_enrichment/cli.py`
- Create: `tests/test_cli.py`
- Modify: `envs/Dockerfile`
- Create: `.github/workflows/python-tests.yml`

**Interfaces:**
- Consumes: Python 3.12 and repository source files.
- Produces: console command `rare-variant-enrichment` with subcommands `prepare-phenotypes`, `classify-chromosome`, `gather`, and `calculate`.

- [ ] **Step 1: Write the failing CLI contract test**

```python
# tests/test_cli.py
import subprocess
import sys


def test_cli_lists_workflow_subcommands():
    result = subprocess.run(
        [sys.executable, "-m", "rare_variant_enrichment.cli", "--help"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    for command in ("prepare-phenotypes", "classify-chromosome", "gather", "calculate"):
        assert command in result.stdout
```

- [ ] **Step 2: Run the test and verify the package is absent**

Run: `python -m pytest tests/test_cli.py -v`

Expected: FAIL because `rare_variant_enrichment` cannot be imported.

- [ ] **Step 3: Add package metadata and the minimal parser**

```toml
# pyproject.toml
[build-system]
requires = ["setuptools>=75"]
build-backend = "setuptools.build_meta"

[project]
name = "rare-variant-enrichment"
version = "0.1.0"
requires-python = ">=3.12"

[project.scripts]
rare-variant-enrichment = "rare_variant_enrichment.cli:main"

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
```

```python
# src/rare_variant_enrichment/cli.py
import argparse


COMMANDS = ("prepare-phenotypes", "classify-chromosome", "gather", "calculate")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rare-variant-enrichment")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in COMMANDS:
        subparsers.add_parser(command)
    return parser


def main() -> int:
    build_parser().parse_args()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Run the CLI test and verify it passes**

Run: `python -m pytest tests/test_cli.py -v`

Expected: PASS with one test.

- [ ] **Step 5: Define the runtime image and Python CI**

```dockerfile
# envs/Dockerfile
FROM python:3.12-slim-bookworm
RUN apt-get update \
    && apt-get install -y --no-install-recommends tabix ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/rare-variant-enrichment
COPY pyproject.toml ./
COPY src ./src
RUN pip install --no-cache-dir .
ENTRYPOINT []
```

Create `.github/workflows/python-tests.yml` to install the `tabix` apt package, install the Python package with `pip install -e .`, run `python -m pytest -v`, and verify `rare-variant-enrichment --help` on pushes and pull requests touching `src/**`, `tests/**`, or `pyproject.toml`.

- [ ] **Step 6: Commit the package foundation**

```bash
git add pyproject.toml src/rare_variant_enrichment tests/test_cli.py envs/Dockerfile .github/workflows/python-tests.yml
git commit -m "build: scaffold rare variant enrichment package"
```

### Task 2: Stream and validate prepare_QTL phenotype BEDs

**Files:**
- Create: `src/rare_variant_enrichment/io.py`
- Create: `src/rare_variant_enrichment/phenotypes.py`
- Create: `tests/test_phenotypes.py`
- Modify: `src/rare_variant_enrichment/cli.py`

**Interfaces:**
- Consumes: wide BED path, VCF sample-list path, chromosome names, z thresholds, and tail mode.
- Produces: `prepare_phenotypes(phenotype_bed, vcf_samples, chromosomes, z_thresholds, tail, feature_output, sample_output, qc_output) -> None`; TSV columns `chrom`, `tss`, `feature_id`; shared-sample text file; phenotype QC JSON.

- [ ] **Step 1: Write failing BED and outlier-classification tests**

```python
# tests/test_phenotypes.py
import json
from pathlib import Path

import pytest

from rare_variant_enrichment.phenotypes import classify_outlier, prepare_phenotypes


def test_classify_outlier_supports_all_tail_modes():
    assert classify_outlier(-3.0, 2.5, "absolute")
    assert classify_outlier(3.0, 2.5, "positive")
    assert classify_outlier(-3.0, 2.5, "negative")
    assert not classify_outlier(-3.0, 2.5, "positive")


def test_prepare_phenotypes_uses_tss_end_and_shared_samples(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\n"
        "chr1\t99\t100\tGENE1\t3.0\t0.0\tNA\n"
        "chr2\t199\t200\tGENE2\t-4.0\t1.0\t2.5\n"
    )
    vcf_samples = tmp_path / "vcf_samples.txt"
    vcf_samples.write_text("S1\nS3\nS4\n")

    feature_output = tmp_path / "features.tsv"
    sample_output = tmp_path / "shared_samples.txt"
    qc_output = tmp_path / "phenotype_qc.json"
    prepare_phenotypes(
        bed, vcf_samples, ["chr1", "chr2"], [2.0, 3.0], "absolute",
        feature_output, sample_output, qc_output,
    )

    assert feature_output.read_text().splitlines()[1] == "chr1\t100\tGENE1"
    assert sample_output.read_text().splitlines() == ["S1", "S3"]
    qc = json.loads(qc_output.read_text())
    assert qc["shared_sample_count"] == 2
    assert qc["bed_only_samples"] == ["S2"]
    assert qc["vcf_only_samples"] == ["S4"]
    assert qc["non_missing_observations"] == 3


def test_prepare_phenotypes_rejects_non_unit_tss_interval(tmp_path: Path):
    bed = tmp_path / "bad.bed"
    bed.write_text("#chr\tstart\tend\tgene_id\tS1\nchr1\t90\t100\tGENE1\t2.0\n")
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\n")
    with pytest.raises(ValueError, match="one base"):
        prepare_phenotypes(bed, samples, ["chr1"], [2.0], "absolute",
                           tmp_path / "f.tsv", tmp_path / "s.txt", tmp_path / "q.json")
```

- [ ] **Step 2: Run the tests and verify the missing-module failure**

Run: `python -m pytest tests/test_phenotypes.py -v`

Expected: FAIL because `phenotypes.py` does not exist.

- [ ] **Step 3: Implement compressed BED streaming and validation**

In `io.py`, add these exact interfaces:

```python
def open_text(path: Path) -> TextIO:
    raw = path.open("rb")
    if raw.read(2) == b"\x1f\x8b":
        raw.seek(0)
        return TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8")
    raw.seek(0)
    return TextIOWrapper(raw, encoding="utf-8")


def read_nonempty_lines(path: Path) -> list[str]:
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
```

In `phenotypes.py`, add:

```python
def classify_outlier(value: float, threshold: float, tail: str) -> bool:
    if tail == "absolute":
        return abs(value) >= threshold
    if tail == "positive":
        return value >= threshold
    if tail == "negative":
        return value <= -threshold
    raise ValueError(f"Unsupported tail mode: {tail}")
```

Implement `prepare_phenotypes` with parameters `phenotype_bed: Path`, `vcf_samples_path: Path`, `chromosomes: Sequence[str]`, `z_thresholds: Sequence[float]`, `tail: str`, `feature_output: Path`, `sample_output: Path`, and `qc_output: Path`, returning `None`.

The implementation must stream rows, preserve BED sample order in the shared-sample file, reject duplicate feature/sample IDs, validate finite numeric thresholds, accept `NA`, `NaN`, `.`, and empty phenotype values as missing, and count non-missing/outlier observations among shared samples only.

- [ ] **Step 4: Wire the `prepare-phenotypes` CLI arguments**

Add required arguments `--phenotype-bed`, `--vcf-samples`, `--chromosomes`, `--z-thresholds`, `--tail`, `--feature-output`, `--sample-output`, and `--qc-output`. Parse comma-separated lists with shared helpers `parse_csv_strings`, `parse_csv_ints`, and `parse_csv_floats` in `cli.py`.

- [ ] **Step 5: Run focused and full tests**

Run: `python -m pytest tests/test_phenotypes.py tests/test_cli.py -v`

Expected: PASS with all phenotype and CLI tests.

- [ ] **Step 6: Commit phenotype preparation**

```bash
git add src/rare_variant_enrichment/io.py src/rare_variant_enrichment/phenotypes.py src/rare_variant_enrichment/cli.py tests/test_phenotypes.py
git commit -m "feat: validate and stream phenotype BEDs"
```

### Task 3: Allele-count classes, multiallelic carriers, and TSS windows

**Files:**
- Create: `src/rare_variant_enrichment/variants.py`
- Create: `tests/test_variants.py`

**Interfaces:**
- Consumes: exact/cumulative AC parameters, VCF record fields, selected sample indices, chromosome TSS records, and maximum distance.
- Produces: immutable `AcClass`, `VariantAllele`, and `FeatureTss` values; merged BED query regions; per-ALT carrier lists; nearby-feature matches.

- [ ] **Step 1: Write failing AC, multiallelic, and distance tests**

```python
# tests/test_variants.py
from rare_variant_enrichment.variants import (
    FeatureTss,
    build_ac_classes,
    merge_tss_windows,
    nearby_features,
    parse_variant_alleles,
)


def test_build_ac_classes_is_stable_and_deduplicated():
    classes = build_ac_classes([1, 2, 2, 3], [1, 3])
    assert [item.label for item in classes] == ["AC=1", "AC=2", "AC=3", "AC<=1", "AC<=3"]


def test_multiallelic_record_assigns_ac_and_carriers_per_alt():
    fields = "chr1\t100\t.\tA\tC,G\t.\tPASS\tAC=1,2\tGT\t0/1\t0/2\t0/2".split("\t")
    alleles = parse_variant_alleles(fields, ["S1", "S2", "S3"], {"S1", "S2", "S3"})
    assert [(a.alt, a.ac, a.carriers) for a in alleles] == [
        ("C", 1, ("S1",)),
        ("G", 2, ("S2", "S3")),
    ]


def test_genotypes_supply_global_ac_but_only_shared_carriers():
    fields = "chr1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1\t./.\t0/1".split("\t")
    allele = parse_variant_alleles(fields, ["S1", "S2", "S3", "S4"], {"S1", "S2", "S3"})[0]
    assert allele.ac == 4
    assert allele.carriers == ("S1", "S2")


def test_windows_merge_and_distance_boundary_is_inclusive():
    features = [FeatureTss("chr1", 100, "G1"), FeatureTss("chr1", 120, "G2")]
    assert merge_tss_windows(features, 10) == [("chr1", 89, 130)]
    assert [f.feature_id for f in nearby_features(features, 110, 10)] == ["G1", "G2"]
```

- [ ] **Step 2: Run the tests and verify they fail for the absent module**

Run: `python -m pytest tests/test_variants.py -v`

Expected: FAIL because `variants.py` does not exist.

- [ ] **Step 3: Implement the variant data model and AC parsing**

```python
@dataclass(frozen=True)
class AcClass:
    label: str
    kind: Literal["exact", "cumulative"]
    value: int

    def contains(self, ac: int) -> bool:
        return ac == self.value if self.kind == "exact" else ac <= self.value


@dataclass(frozen=True)
class VariantAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    ac: int
    carriers: tuple[str, ...]


@dataclass(frozen=True, order=True)
class FeatureTss:
    chrom: str
    tss: int
    feature_id: str
```

Implement `build_ac_classes`, `parse_variant_alleles`, GT allele-index counting for phased/unphased genotypes, missing-call accounting, and a clear error when `INFO/AC` cardinality does not match the ALT count. Genotype-derived AC must count all VCF samples, while `VariantAllele.carriers` contains only IDs from the supplied shared-sample set.

- [ ] **Step 4: Implement merged BED windows and nearby-feature lookup**

Implement `merge_tss_windows(features: Sequence[FeatureTss], max_distance: int) -> list[tuple[str, int, int]]` and `nearby_features(features: Sequence[FeatureTss], position: int, max_distance: int) -> list[FeatureTss]`.

Use zero-based half-open BED query intervals while preserving the inclusive VCF-distance rule. For TSS 100 and distance 10, the query interval is BED `[89, 110)`, which captures one-based VCF positions 90 through 110.

- [ ] **Step 5: Run variant and regression tests**

Run: `python -m pytest tests/test_variants.py tests/test_phenotypes.py -v`

Expected: PASS.

- [ ] **Step 6: Commit variant primitives**

```bash
git add src/rare_variant_enrichment/variants.py tests/test_variants.py
git commit -m "feat: classify rare alleles and TSS windows"
```

### Task 4: One-pass chromosome tabix classification

**Files:**
- Modify: `src/rare_variant_enrichment/variants.py`
- Modify: `src/rare_variant_enrichment/cli.py`
- Create: `tests/test_chromosome_classification.py`
- Create: `tests/fixtures/rare_variants.vcf`
- Create: `tests/fixtures/features.tsv`
- Create: `tests/fixtures/shared_samples.txt`

**Interfaces:**
- Consumes: adjacent VCF/`.tbi`, compact feature/TSS TSV, shared samples, one chromosome, AC arrays, and maximum distance.
- Produces: `classify_chromosome` with the exact parameters listed in Step 3; merged-region BED, carrier-minimum TSV with columns `sample_id`, `feature_id`, `ac_class`, `minimum_distance_bp`; chromosome QC JSON.

- [ ] **Step 1: Add a failing real-tabix integration test**

```python
# tests/test_chromosome_classification.py
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from rare_variant_enrichment.variants import classify_chromosome


@pytest.mark.skipif(shutil.which("tabix") is None or shutil.which("bgzip") is None,
                    reason="htslib executables are required")
def test_classify_chromosome_extracts_once_and_keeps_minimum_distance(tmp_path: Path):
    plain_vcf = Path("tests/fixtures/rare_variants.vcf")
    vcf = tmp_path / "rare_variants.vcf.gz"
    with vcf.open("wb") as output:
        subprocess.run(["bgzip", "-c", str(plain_vcf)], stdout=output, check=True)
    subprocess.run(["tabix", "-p", "vcf", str(vcf)], check=True)

    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc_path = tmp_path / "qc.json"
    classify_chromosome(
        vcf, Path("tests/fixtures/features.tsv"), Path("tests/fixtures/shared_samples.txt"),
        "chr1", [1, 2, 3], [1, 3], 100,
        carriers, regions, qc_path,
    )

    rows = carriers.read_text().splitlines()
    assert "S1\tGENE1\tAC=1\t0" in rows
    assert "S2\tGENE1\tAC=2\t50" in rows
    qc = json.loads(qc_path.read_text())
    assert qc["tabix_query_count"] == 1
```

The fixture VCF must include: an AC=1 ALT at the GENE1 TSS carried by S1; an AC=2 ALT 50 bp away carried by S2; overlapping GENE1/GENE2 maximum windows; one missing genotype; and one variant outside the maximum window. `features.tsv` must contain `chr1\t100\tGENE1` and `chr1\t120\tGENE2` after its header.

- [ ] **Step 2: Run the integration test and verify the missing function failure**

Run: `python -m pytest tests/test_chromosome_classification.py -v`

Expected: FAIL because `classify_chromosome` is not defined.

- [ ] **Step 3: Implement streaming tabix classification**

Implement `classify_chromosome` with parameters `vcf_path: Path`, `features_path: Path`, `shared_samples_path: Path`, `chromosome: str`, `exact_ac: Sequence[int]`, `cumulative_ac_max: Sequence[int]`, `max_distance: int`, `carrier_output: Path`, `regions_output: Path`, and `qc_output: Path`, returning `None`.

Write merged windows to `regions_output`, launch exactly one `tabix -h -R regions_output vcf_path` subprocess, parse its header and records from stdout, match each ALT to nearby features, and update a dictionary keyed by `(sample_id, feature_id, ac_class)` only when the new absolute distance is smaller. Reject a requested chromosome absent from the VCF; when the chromosome has no BED features, emit header-only carrier/region outputs and QC with `tabix_query_count = 0`. Report extracted records, ALT alleles, missing genotypes, variant-feature pairs, and emitted keys.

- [ ] **Step 4: Wire the `classify-chromosome` CLI**

Add arguments matching the function parameters. The WDL task will create adjacent `variants.vcf.gz` and `variants.vcf.gz.tbi` symlinks before invoking the CLI.

- [ ] **Step 5: Run focused and full Python tests**

Run: `python -m pytest tests/test_chromosome_classification.py tests/test_variants.py -v`

Expected: PASS; the integration test is skipped only outside an htslib-equipped environment.

- [ ] **Step 6: Commit chromosome classification**

```bash
git add src/rare_variant_enrichment/variants.py src/rare_variant_enrichment/cli.py tests/test_chromosome_classification.py tests/fixtures
git commit -m "feat: classify chromosome carriers with tabix"
```

### Task 5: Gather minimum carrier distances and chromosome QC

**Files:**
- Create: `src/rare_variant_enrichment/aggregation.py`
- Create: `tests/test_aggregation.py`
- Modify: `src/rare_variant_enrichment/cli.py`

**Interfaces:**
- Consumes: one or more chromosome carrier TSVs and chromosome QC JSON files.
- Produces: `gather_outputs(carrier_paths, qc_paths, carrier_output, qc_output) -> None`; globally deduplicated carrier-minimum TSV and chromosome QC TSV.

- [ ] **Step 1: Write the failing gather test**

```python
# tests/test_aggregation.py
from pathlib import Path

from rare_variant_enrichment.aggregation import gather_outputs


def test_gather_keeps_global_minimum_distance(tmp_path: Path):
    first = tmp_path / "chr1.tsv"
    second = tmp_path / "chr2.tsv"
    header = "sample_id\tfeature_id\tac_class\tminimum_distance_bp\n"
    first.write_text(header + "S1\tGENE1\tAC=1\t100\nS2\tGENE2\tAC<=3\t20\n")
    second.write_text(header + "S1\tGENE1\tAC=1\t40\n")
    q1 = tmp_path / "chr1.json"
    q2 = tmp_path / "chr2.json"
    q1.write_text('{"chromosome":"chr1","vcf_records_extracted":2}')
    q2.write_text('{"chromosome":"chr2","vcf_records_extracted":1}')

    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "chromosome_qc.tsv"
    gather_outputs([first, second], [q1, q2], carrier_output, qc_output)

    assert "S1\tGENE1\tAC=1\t40" in carrier_output.read_text().splitlines()
    assert qc_output.read_text().splitlines()[1].startswith("chr1\t")
```

- [ ] **Step 2: Run the test and verify the missing-module failure**

Run: `python -m pytest tests/test_aggregation.py -v`

Expected: FAIL because `aggregation.py` does not exist.

- [ ] **Step 3: Implement deterministic gather logic**

Implement `gather_outputs` with parameters `carrier_paths: Sequence[Path]`, `qc_paths: Sequence[Path]`, `carrier_output: Path`, and `qc_output: Path`, returning `None`.

Validate identical carrier headers, retain the minimum integer distance per key, sort output by feature/sample/class, validate that every QC object has one chromosome, and write the union of QC keys as stable TSV columns.

- [ ] **Step 4: Wire the `gather` CLI and run tests**

Use repeatable `--carrier-input` and `--qc-input` arguments plus `--carrier-output` and `--qc-output`.

Run: `python -m pytest tests/test_aggregation.py tests/test_cli.py -v`

Expected: PASS.

- [ ] **Step 5: Commit gather support**

```bash
git add src/rare_variant_enrichment/aggregation.py src/rare_variant_enrichment/cli.py tests/test_aggregation.py
git commit -m "feat: gather carrier distances across chromosomes"
```

### Task 6: Enrichment statistics and final report

**Files:**
- Create: `src/rare_variant_enrichment/statistics.py`
- Create: `tests/test_statistics.py`
- Modify: `src/rare_variant_enrichment/cli.py`

**Interfaces:**
- Consumes: phenotype BED, shared samples, deduplicated carrier minima, exact/cumulative AC definitions, z thresholds, distance thresholds, and tail mode.
- Produces: `fisher_exact_two_sided`, `benjamini_hochberg`, and `calculate_enrichment` with the exact parameters listed in Step 4; enrichment TSV and run-summary JSON.

- [ ] **Step 1: Write failing mathematical tests**

```python
# tests/test_statistics.py
import csv
import json
from pathlib import Path

import pytest

from rare_variant_enrichment.statistics import (
    benjamini_hochberg,
    calculate_enrichment,
    fisher_exact_two_sided,
)


def test_fisher_exact_matches_hand_checked_table():
    assert fisher_exact_two_sided(1, 9, 11, 3) == pytest.approx(0.002759456, rel=1e-6)


def test_bh_is_monotone_in_rank_order():
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03])
    assert adjusted == pytest.approx([0.03, 0.04, 0.04])


def test_calculate_enrichment_counts_distance_specific_carriers(tmp_path: Path):
    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\n"
        "chr1\t99\t100\tGENE1\t3\t0\t-3\t0\n"
    )
    samples = tmp_path / "samples.txt"
    samples.write_text("S1\nS2\nS3\nS4\n")
    carriers = tmp_path / "carriers.tsv"
    carriers.write_text(
        "sample_id\tfeature_id\tac_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\t10\nS2\tGENE1\tAC=1\t100\n"
    )
    output = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"
    calculate_enrichment(bed, samples, carriers, [1], [], [2.0], [50, 100],
                         "absolute", output, summary)

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows[0]["distance_bp"] == "50"
    assert (rows[0]["n11"], rows[0]["n10"], rows[0]["n01"], rows[0]["n00"]) == ("1", "1", "0", "2")
    assert (rows[1]["n11"], rows[1]["n10"], rows[1]["n01"], rows[1]["n00"]) == ("1", "1", "1", "1")
    assert "screening" in json.loads(summary.read_text())["statistical_limitation"].lower()
```

- [ ] **Step 2: Run the tests and verify the missing-module failure**

Run: `python -m pytest tests/test_statistics.py -v`

Expected: FAIL because `statistics.py` does not exist.

- [ ] **Step 3: Implement exact tests and correction without SciPy**

```python
def fisher_exact_two_sided(n11: int, n10: int, n01: int, n00: int) -> float:
    row1 = n11 + n10
    row2 = n01 + n00
    col1 = n11 + n01
    total = row1 + row2

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def probability(x: int) -> float:
        return math.exp(log_choose(col1, x) + log_choose(total - col1, row1 - x)
                        - log_choose(total, row1))

    observed = probability(n11)
    lower = max(0, row1 - (total - col1))
    upper = min(row1, col1)
    return min(1.0, sum(probability(x) for x in range(lower, upper + 1)
                        if probability(x) <= observed * (1.0 + 1e-12)))


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted
```

Compute hypergeometric probabilities from fixed margins using `math.lgamma`, sum tables with probability no greater than the observed table, clamp floating-point results to `[0, 1]`, and apply reverse cumulative minima for BH-adjusted values.

- [ ] **Step 4: Implement streaming enrichment calculation**

Implement `calculate_enrichment` with parameters `phenotype_bed: Path`, `shared_samples_path: Path`, `carriers_path: Path`, `exact_ac: Sequence[int]`, `cumulative_ac_max: Sequence[int]`, `z_thresholds: Sequence[float]`, `distance_thresholds: Sequence[int]`, `tail: str`, `output_tsv: Path`, and `output_json: Path`, returning `None`.

Build the complete class list from `exact_ac` and `cumulative_ac_max` so classes with zero carriers still emit rows. Group carrier minima by feature, stream one BED row at a time, compute total/outlier observations once per z threshold, update carrier and carrier-outlier counts for each applicable AC/distance combination, derive all four cells, calculate ratios and p-values, then apply BH across every emitted row. Use `NA` for undefined uncorrected ratios and report both the uncorrected odds ratio and a 0.5-corrected odds ratio.

Write these TSV columns in order:

```text
z_threshold	tail	ac_class	ac_kind	ac_value	distance_bp	total_observations	outlier_observations	nonoutlier_observations	n11	n10	n01	n00	outlier_carrier_rate	nonoutlier_carrier_rate	carrier_rate_ratio	odds_ratio	odds_ratio_corrected_0_5	fisher_p_value	fisher_fdr_bh
```

- [ ] **Step 5: Wire `calculate` CLI arguments and verify tests**

Run: `python -m pytest tests/test_statistics.py tests/test_phenotypes.py tests/test_aggregation.py -v`

Expected: PASS.

- [ ] **Step 6: Commit enrichment statistics**

```bash
git add src/rare_variant_enrichment/statistics.py src/rare_variant_enrichment/cli.py tests/test_statistics.py
git commit -m "feat: calculate distance-stratified enrichment"
```

### Task 7: WDL scatter/gather orchestration

**Files:**
- Create: `workflows/rare_variant_enrichment.wdl`
- Create: `tests/test_wdl_contract.py`

**Interfaces:**
- Consumes: phenotype BED, rare VCF, optional `.tbi`, chromosomes, z thresholds, exact/cumulative AC arrays, distance thresholds, tail, Docker image, and runtime resources.
- Produces: enrichment TSV/JSON, chromosome QC TSV, deduplicated carrier TSV, prepared VCF index, and per-chromosome merged-window BEDs.

- [ ] **Step 1: Write the failing WDL contract test**

```python
# tests/test_wdl_contract.py
from pathlib import Path


def test_wdl_exposes_required_tasks_inputs_and_outputs():
    text = Path("workflows/rare_variant_enrichment.wdl").read_text()
    for task in ("PrepareVcfIndex", "PreparePhenotypes", "ClassifyChromosome",
                 "GatherCarrierPairs", "CalculateEnrichment"):
        assert f"task {task}" in text
    for name in ("File phenotype_bed", "File rare_variant_vcf", "File? rare_variant_vcf_tbi",
                 "Array[Int] distance_thresholds_bp", "File enrichment_tsv",
                 "File chromosome_qc_tsv"):
        assert name in text
    assert "scatter (chromosome in chromosomes)" in text
```

- [ ] **Step 2: Run the contract test and verify the missing-file failure**

Run: `python -m pytest tests/test_wdl_contract.py -v`

Expected: FAIL because the WDL file does not exist.

- [ ] **Step 3: Implement `PrepareVcfIndex` and `PreparePhenotypes` tasks**

Use WDL version 1.0. Define required workflow inputs `phenotype_bed`, `rare_variant_vcf`, and `chromosomes`; optional `rare_variant_vcf_tbi`; and these defaults:

```wdl
Array[Float] z_thresholds = [2.0, 3.0, 4.0, 5.0]
Array[Int] exact_allele_counts = [1, 2, 3, 4, 5]
Array[Int] cumulative_allele_count_maxima = [1, 2, 3, 5, 10]
Array[Int] distance_thresholds_bp = [1000, 10000, 100000, 1000000]
String outlier_tail = "absolute"
String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
Int prepare_cpu = 2
Int prepare_memory_gb = 8
Int prepare_disk_gb = 50
Int scatter_cpu = 2
Int scatter_memory_gb = 8
Int scatter_disk_gb = 20
Int gather_cpu = 2
Int gather_memory_gb = 16
Int gather_disk_gb = 50
```

`PrepareVcfIndex` must symlink the localized VCF to `variants.vcf.gz`, symlink a supplied index or run `tabix -p vcf`, validate with `tabix -l`, reject requested contigs absent from the index, and write ordered VCF samples from the `#CHROM` header. `PreparePhenotypes` must call the corresponding CLI subcommand and emit `features.tsv`, `shared_samples.txt`, and `phenotype_qc.json`.

- [ ] **Step 4: Implement chromosome scatter and gather tasks**

Compute `Int maximum_distance_bp = max(distance_thresholds_bp)`, scatter `ClassifyChromosome` over the required `chromosomes` input, create adjacent VCF/index symlinks inside every scattered task, and call `classify-chromosome` once. Pass scattered carrier/QC arrays to `GatherCarrierPairs`, using repeatable CLI arguments generated with WDL `sep` expressions.

- [ ] **Step 5: Implement final calculation and workflow outputs**

Call `CalculateEnrichment` with the BED, shared sample list, gathered carrier minima, exact/cumulative AC arrays, z thresholds, all distance thresholds, and tail. Expose these workflow outputs with exact names:

```wdl
output {
    File enrichment_tsv = CalculateEnrichment.enrichment_tsv
    File enrichment_json = CalculateEnrichment.enrichment_json
    File chromosome_qc_tsv = GatherCarrierPairs.chromosome_qc_tsv
    File carrier_minimum_distances_tsv = GatherCarrierPairs.carrier_minimum_distances_tsv
    File generated_or_validated_vcf_tbi = PrepareVcfIndex.vcf_tbi
    Array[File] chromosome_query_regions = ClassifyChromosome.query_regions_bed
}
```

- [ ] **Step 6: Run Python contract and WDL static validation**

Run: `python -m pytest tests/test_wdl_contract.py -v`

Run: `miniwdl check workflows/rare_variant_enrichment.wdl`

Expected: both commands exit 0; miniwdl reports no errors.

- [ ] **Step 7: Commit the workflow**

```bash
git add workflows/rare_variant_enrichment.wdl tests/test_wdl_contract.py
git commit -m "feat: add chromosome-scattered enrichment WDL"
```

### Task 8: End-to-end fixture, repository metadata, and user documentation

**Files:**
- Create: `tests/test_end_to_end.py`
- Create: `examples/rare_variant_enrichment.inputs.json`
- Modify: `README.md`
- Modify: `.dockstore.yml`
- Modify: `.github/workflows/docker-image.yml`

**Interfaces:**
- Consumes: all package commands and the WDL from Tasks 1–7.
- Produces: a runnable miniature CLI test, documented WDL inputs/outputs, Dockstore registration, and build/publish-ready image automation.

- [ ] **Step 1: Write a failing miniature pipeline test**

```python
# tests/test_end_to_end.py
import csv
import shutil
from pathlib import Path

import pytest

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.phenotypes import prepare_phenotypes
from rare_variant_enrichment.statistics import calculate_enrichment
from rare_variant_enrichment.variants import classify_chromosome


@pytest.mark.skipif(shutil.which("tabix") is None or shutil.which("bgzip") is None,
                    reason="htslib executables are required")
def test_miniature_pipeline_emits_every_threshold_combination(prepared_fixture, tmp_path: Path):
    feature_tsv = tmp_path / "features.tsv"
    shared = tmp_path / "shared.txt"
    phenotype_qc = tmp_path / "phenotype_qc.json"
    prepare_phenotypes(prepared_fixture.bed, prepared_fixture.samples, ["chr1"],
                       [2.0, 3.0], "absolute", feature_tsv, shared, phenotype_qc)

    chromosome_carriers = tmp_path / "chr1.carriers.tsv"
    regions = tmp_path / "chr1.regions.bed"
    chromosome_qc = tmp_path / "chr1.qc.json"
    classify_chromosome(prepared_fixture.vcf_gz, feature_tsv, shared, "chr1",
                        [1, 2], [1, 2], 100, chromosome_carriers, regions, chromosome_qc)

    all_carriers = tmp_path / "carriers.tsv"
    all_qc = tmp_path / "chromosome_qc.tsv"
    gather_outputs([chromosome_carriers], [chromosome_qc], all_carriers, all_qc)

    enrichment = tmp_path / "enrichment.tsv"
    summary = tmp_path / "summary.json"
    calculate_enrichment(prepared_fixture.bed, shared, all_carriers,
                         [1, 2], [1, 2], [2.0, 3.0], [10, 100],
                         "absolute", enrichment, summary)
    rows = list(csv.DictReader(enrichment.open(), delimiter="\t"))
    assert len(rows) == 2 * 4 * 2
```

Add a `prepared_fixture` pytest fixture in `tests/conftest.py` that bgzips and indexes the committed plain VCF fixture, creates a four-sample phenotype BED, and writes the VCF sample list.

- [ ] **Step 2: Run the test and verify its first real integration failure**

Run: `python -m pytest tests/test_end_to_end.py -v`

Expected: FAIL until the fixture and all command interfaces are connected exactly as specified.

- [ ] **Step 3: Complete the fixture and make the end-to-end test pass**

Run: `python -m pytest tests/test_end_to_end.py -v`

Expected: PASS in an htslib-equipped environment and exactly 16 enrichment rows.

- [ ] **Step 4: Register and document the workflow**

Set `.dockstore.yml` `primaryDescriptorPath` to `/workflows/rare_variant_enrichment.wdl`. Replace the template README with:

- The biological question and pooled gene–sample analysis unit.
- The prepare_QTL BED schema and one-base TSS convention.
- Required VCF properties and optional index behavior.
- Default/recommended z, exact AC, cumulative AC, and distance arrays.
- A complete `miniwdl run workflows/rare_variant_enrichment.wdl -i examples/rare_variant_enrichment.inputs.json` example.
- Definitions for every enrichment TSV column.
- Sample-overlap and missing-genotype behavior.
- The repeated-observation/Fisher screening limitation.
- Future annotation and watershed extension points.

Create `examples/rare_variant_enrichment.inputs.json` with relative input names `phenotypes.scaled.residualized.bed.gz` and `rare_variants.vcf.gz`, no index field (demonstrating automatic generation), `chr1` through `chr22` and `chrX`, z thresholds `[2.0, 3.0, 4.0, 5.0]`, exact AC `[1, 2, 3, 4, 5]`, cumulative maxima `[1, 2, 3, 5, 10]`, distances `[1000, 10000, 100000, 1000000]`, tail `absolute`, and the default GHCR image.

- [ ] **Step 5: Update image build triggers and publishing behavior**

Extend `.github/workflows/docker-image.yml` path filters to include `pyproject.toml`, `src/**`, and `tests/**`. Build without pushing on pull requests; push the generated GHCR tags on `main` and manual dispatch so the WDL's default runtime image exists.

- [ ] **Step 6: Run complete verification**

Run: `python -m pytest -v`

Run: `miniwdl check workflows/rare_variant_enrichment.wdl`

Run: `docker build -f envs/Dockerfile -t rare-variant-enrichment:test .`

Run: `docker run --rm rare-variant-enrichment:test rare-variant-enrichment --help`

Expected: all tests pass, WDL validation exits 0, the image builds, and the container help lists all four subcommands.

- [ ] **Step 7: Review scope and commit the completed first pass**

Confirm `git diff --check` is clean, no cohort data or credentials are tracked, and the implementation contains no annotation or watershed logic.

```bash
git add README.md .dockstore.yml .github/workflows/docker-image.yml examples tests/test_end_to_end.py tests/conftest.py
git commit -m "docs: document rare variant enrichment workflow"
```
