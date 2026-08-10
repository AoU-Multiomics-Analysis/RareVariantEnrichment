import json
from io import StringIO
import shutil
import subprocess
from pathlib import Path

import pytest

from rare_variant_enrichment.variants import classify_chromosome


class _TabixProcess:
    def __init__(self, output: str):
        self.stdout = StringIO(output)
        self.returncode = 0

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.stdout.close()

    def wait(self) -> int:
        return self.returncode


@pytest.mark.skipif(
    shutil.which("tabix") is None or shutil.which("bgzip") is None,
    reason="htslib executables are required",
)
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
        vcf,
        Path("tests/fixtures/features.tsv"),
        Path("tests/fixtures/shared_samples.txt"),
        "chr1",
        [1, 2, 3],
        [1, 3],
        100,
        carriers,
        regions,
        qc_path,
    )

    rows = carriers.read_text().splitlines()
    assert "S1\tGENE1\tAC=1\t0" in rows
    assert "S2\tGENE1\tAC=2\t50" in rows
    qc = json.loads(qc_path.read_text())
    assert qc["tabix_query_count"] == 1


def test_classify_chromosome_without_features_writes_empty_outputs_without_tabix(
    tmp_path: Path,
):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr2\t100\tGENE2\n")
    samples = tmp_path / "shared_samples.txt"
    samples.write_text("S1\n")
    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc_path = tmp_path / "qc.json"

    classify_chromosome(
        vcf, features, samples, "chr1", [1], [1], 100, carriers, regions, qc_path
    )

    assert carriers.read_text() == "sample_id\tfeature_id\tac_class\tminimum_distance_bp\n"
    assert regions.read_text() == "#chrom\tstart\tend\n"
    assert json.loads(qc_path.read_text()) == {
        "alt_alleles": 0,
        "chromosome": "chr1",
        "emitted_keys": 0,
        "extracted_records": 0,
        "feature_count": 0,
        "merged_region_count": 0,
        "missing_genotypes": 0,
        "tabix_query_count": 0,
        "variant_feature_pairs": 0,
    }


def test_classify_chromosome_rejects_requested_vcf_chromosome_before_tabix(tmp_path: Path):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr2\t100\tGENE2\n")
    samples = tmp_path / "shared_samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="Requested chromosome is absent from VCF: chr2"):
        classify_chromosome(
            vcf,
            features,
            samples,
            "chr2",
            [1],
            [1],
            100,
            tmp_path / "carriers.tsv",
            tmp_path / "regions.bed",
            tmp_path / "qc.json",
        )


def test_classify_chromosome_processes_streamed_record_without_htslib(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
    samples = tmp_path / "shared_samples.txt"
    samples.write_text("S1\n")
    tabix_output = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t100\t.\tA\tC\t.\tPASS\tAC=1\tGT\t0/1\n"
    )
    monkeypatch.setattr(
        "rare_variant_enrichment.variants.subprocess.Popen",
        lambda *_args, **_kwargs: _TabixProcess(tabix_output),
    )
    carriers = tmp_path / "carriers.tsv"
    qc_path = tmp_path / "qc.json"

    classify_chromosome(
        vcf,
        features,
        samples,
        "chr1",
        [1],
        [],
        100,
        carriers,
        tmp_path / "regions.bed",
        qc_path,
    )

    assert carriers.read_text().splitlines() == [
        "sample_id\tfeature_id\tac_class\tminimum_distance_bp",
        "S1\tGENE1\tAC=1\t0",
    ]
    assert json.loads(qc_path.read_text())["extracted_records"] == 1


def test_classify_chromosome_rejects_vcf_without_contig_metadata(tmp_path: Path):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\nchr1\t100\tGENE1\n")
    samples = tmp_path / "shared_samples.txt"
    samples.write_text("S1\n")

    with pytest.raises(ValueError, match="VCF header does not declare contigs"):
        classify_chromosome(
            vcf,
            features,
            samples,
            "chr1",
            [1],
            [],
            100,
            tmp_path / "carriers.tsv",
            tmp_path / "regions.bed",
            tmp_path / "qc.json",
        )
