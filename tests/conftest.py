from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from rare_variant_enrichment.annotations import VatSchema


@dataclass(frozen=True)
class PreparedFixture:
    bed: Path
    features: Path
    samples: Path
    vcf_gz: Path
    vcf_tbi: Path
    vat_bgz: Path
    vat_tbi: Path
    vat_schema: Path


@pytest.fixture
def prepared_fixture(tmp_path: Path) -> PreparedFixture:
    plain_vcf = Path(__file__).parent / "fixtures" / "rare_variants.vcf"
    plain_vat = Path(__file__).parent / "fixtures" / "transcript_annotations.tsv"
    vcf_gz = tmp_path / "rare_variants.vcf.gz"
    vat_bgz = tmp_path / "transcript_annotations.tsv.bgz"
    if shutil.which("bgzip") is not None and shutil.which("tabix") is not None:
        _bgzip(plain_vcf, vcf_gz)
        subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)
        _bgzip(plain_vat, vat_bgz)
        subprocess.run(
            ["tabix", "-f", "-S", "1", "-s", "1", "-b", "2", "-e", "2", str(vat_bgz)],
            check=True,
        )
    else:
        _prepare_tabix_with_container(plain_vcf, vcf_gz, tmp_path, ["-p", "vcf"])
        _prepare_tabix_with_container(
            plain_vat,
            vat_bgz,
            tmp_path,
            ["-f", "-S", "1", "-s", "1", "-b", "2", "-e", "2"],
        )

    vat_schema = tmp_path / "vat_schema.json"
    header = plain_vat.read_text().splitlines()[0].split("\t")
    VatSchema.from_header(header).write_json(vat_schema)

    bed = tmp_path / "phenotypes.bed"
    bed.write_text(
        "#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\n"
        "chr1\t99\t100\tENSG000001.1\t3.5\t0\t0\t4\n"
        "chr1\t149\t150\tENSG000002.2\t3.5\t0\t0\t4\n"
        "chr1\t299\t300\tENSG000003.3\t2.5\t0\t0\t4\n"
        "chr1\t99\t100\tENSG000004.4\tNA\t0\t0\t4\n"
    )
    features = tmp_path / "features.tsv"
    features.write_text(
        "chrom\ttss\tfeature_id\n"
        "chr1\t100\tENSG000001.1\n"
        "chr1\t150\tENSG000002.2\n"
        "chr1\t300\tENSG000003.3\n"
        "chr1\t100\tENSG000004.4\n"
    )
    samples = tmp_path / "vcf_samples.txt"
    samples.write_text("S1\nS2\nS3\n")
    return PreparedFixture(
        bed,
        features,
        samples,
        vcf_gz,
        Path(f"{vcf_gz}.tbi"),
        vat_bgz,
        Path(f"{vat_bgz}.tbi"),
        vat_schema,
    )


def _bgzip(source: Path, output: Path) -> None:
    with output.open("wb") as output_handle:
        subprocess.run(["bgzip", "-c", str(source)], stdout=output_handle, check=True)


def _prepare_tabix_with_container(
    source: Path,
    compressed: Path,
    directory: Path,
    index_arguments: list[str],
) -> None:
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
    with source.open("rb") as input_handle, compressed.open("wb") as output_handle:
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
            *index_arguments,
            compressed.name,
        ],
        check=True,
    )


def _fixture_unavailable(message: str) -> None:
    if os.environ.get("RARE_VARIANT_ENRICHMENT_REQUIRE_WDL_RUNTIME") == "1":
        pytest.fail(message)
    pytest.skip(message)
