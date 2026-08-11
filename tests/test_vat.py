import gzip
from pathlib import Path
import shutil
import subprocess

import pytest

from rare_variant_enrichment.annotations import VatSchema
from rare_variant_enrichment.vat import prepare_vat


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


def write_vat(tmp_path: Path) -> Path:
    vat = tmp_path / "annotations.tsv.gz"
    with gzip.open(vat, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tpos\tref\talt\tgene_id\tconsequence\tgvs_max_af\tLoF\n")
        handle.write("chr1\t10\tA\tG\tENSG1\tmissense_variant\t0.001\tHC\n")
    return vat


@requires_htslib
def test_prepare_vat_generates_and_validates_generic_tabix_index(tmp_path: Path):
    vat = bgzip_fixture(tmp_path, "tests/fixtures/transcript_annotations.tsv")
    schema_path = tmp_path / "vat_schema.json"

    schema = prepare_vat(vat, ["chr1"], schema_path)

    assert Path(f"{vat}.tbi").is_file()
    assert schema.lof is not None
    assert VatSchema.read_json(schema_path) == schema


def test_prepare_vat_rejects_requested_contigs_absent_from_index(tmp_path: Path, monkeypatch):
    vat = write_vat(tmp_path)
    Path(f"{vat}.tbi").touch()
    monkeypatch.setattr(
        "rare_variant_enrichment.vat.subprocess.run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0, "chr1\n"),
    )

    with pytest.raises(ValueError, match=r"^Requested chromosomes are absent from VAT index: chr2$"):
        prepare_vat(vat, ["chr1", "chr2"], tmp_path / "schema.json")


def test_prepare_vat_reports_missing_tabix_or_bgzip_support(tmp_path: Path, monkeypatch):
    vat = write_vat(tmp_path)

    def unavailable(*_args, **_kwargs):
        raise OSError("tabix not found")

    monkeypatch.setattr("rare_variant_enrichment.vat.subprocess.run", unavailable)

    with pytest.raises(ValueError, match=r"^Unable to prepare VAT index"):
        prepare_vat(vat, ["chr1"], tmp_path / "schema.json")


def test_prepare_vat_rejects_supplied_index_that_tabix_cannot_list(tmp_path: Path, monkeypatch):
    vat = write_vat(tmp_path)
    Path(f"{vat}.tbi").touch()

    def invalid_index(command, **_kwargs):
        assert command == ["tabix", "-l", str(vat)]
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr("rare_variant_enrichment.vat.subprocess.run", invalid_index)

    with pytest.raises(ValueError, match=r"^Unable to validate VAT index"):
        prepare_vat(vat, ["chr1"], tmp_path / "schema.json")
