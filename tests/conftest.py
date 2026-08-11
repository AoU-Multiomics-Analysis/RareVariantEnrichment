from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

import pytest


@dataclass(frozen=True)
class PreparedFixture:
    bed: Path
    samples: Path
    vcf_gz: Path
    vcf_tbi: Path


@pytest.fixture
def prepared_fixture(tmp_path: Path) -> PreparedFixture:
    if shutil.which("bgzip") is None or shutil.which("tabix") is None:
        pytest.skip("htslib executables are required")
    plain_vcf = Path(__file__).parent / "fixtures" / "rare_variants.vcf"
    vcf_gz = tmp_path / "rare_variants.vcf.gz"
    with vcf_gz.open("wb") as output:
        subprocess.run(["bgzip", "-c", str(plain_vcf)], stdout=output, check=True)
    subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)

    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\n"
        "chr1\t99\t100\tGENE1\t3.5\t0\t0\t4\n"
        "chr1\t99\t100\tGENE2\t3.5\t0\t0\t4\n"
        "chr1\t99\t100\tGENE3\t2.5\t0\t0\t4\n"
        "chr1\t99\t100\tGENE4\tNA\t0\t0\t4\n"
    )
    samples = tmp_path / "vcf_samples.txt"
    samples.write_text("S1\nS2\nS3\n")
    return PreparedFixture(bed, samples, vcf_gz, Path(f"{vcf_gz}.tbi"))
