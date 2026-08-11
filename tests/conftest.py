from dataclasses import dataclass
import os
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
    plain_vcf = Path(__file__).parent / "fixtures" / "rare_variants.vcf"
    vcf_gz = tmp_path / "rare_variants.vcf.gz"
    if shutil.which("bgzip") is not None and shutil.which("tabix") is not None:
        with vcf_gz.open("wb") as output:
            subprocess.run(["bgzip", "-c", str(plain_vcf)], stdout=output, check=True)
        subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)
    else:
        _prepare_vcf_with_container(plain_vcf, vcf_gz, tmp_path)

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


def _prepare_vcf_with_container(plain_vcf: Path, vcf_gz: Path, directory: Path) -> None:
    image = os.environ.get(
        "RARE_VARIANT_ENRICHMENT_TEST_IMAGE", "rare-variant-enrichment:test"
    )
    if shutil.which("docker") is None:
        _fixture_unavailable("htslib executables and Docker are unavailable")
    inspected = subprocess.run(
        ["docker", "image", "inspect", image],
        text=True,
        capture_output=True,
        check=False,
    )
    if inspected.returncode != 0:
        _fixture_unavailable(f"Build {image} to prepare the indexed VCF fixture")
    with plain_vcf.open("rb") as input_handle, vcf_gz.open("wb") as output_handle:
        subprocess.run(
            ["docker", "run", "--rm", "-i", image, "bgzip", "-c"],
            stdin=input_handle,
            stdout=output_handle,
            check=True,
        )
    subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{directory.resolve()}:/data",
            "-w",
            "/data",
            image,
            "tabix",
            "-p",
            "vcf",
            vcf_gz.name,
        ],
        check=True,
    )


def _fixture_unavailable(message: str) -> None:
    if os.environ.get("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME") == "1":
        pytest.fail(message)
    pytest.skip(message)
