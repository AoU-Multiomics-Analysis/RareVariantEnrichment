import csv
import json
from io import StringIO
from pathlib import Path
import shutil

import pytest

from rare_variant_enrichment.annotations import VatSchema
from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.variants import classify_chromosome, parse_variant_alleles


VAT_HEADER = (
    "chrom\tpos\tref\talt\trsid\tgene_id\tgene_symbol\ttranscript\t"
    "is_canonical_transcript\tconsequence\taa_change\trevel\tLoF\tLoF_filter\t"
    "LoF_flags\tLoF_info\tgvs_max_af\tgvs_max_subpop"
)

requires_htslib = pytest.mark.skipif(
    shutil.which("tabix") is None or shutil.which("bgzip") is None,
    reason="htslib executables are required",
)


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


class _TabixRouter:
    def __init__(self, rows_by_path: dict[Path, list[str]]):
        self.rows_by_path = {str(path): rows for path, rows in rows_by_path.items()}
        self.calls: list[tuple[str, str]] = []

    def __call__(self, command, **_kwargs):
        path = command[1]
        region = command[2]
        self.calls.append((path, region))
        chromosome, interval = region.split(":")
        start_text, end_text = interval.split("-")
        start, end = int(start_text), int(end_text)
        selected = []
        for row in self.rows_by_path[path]:
            fields = row.split("\t")
            if fields[0] == chromosome and start <= int(fields[1]) <= end:
                selected.append(row)
        output = "\n".join(selected)
        return _TabixProcess(output + ("\n" if output else ""))


def _write_inputs(tmp_path: Path, features_text: str, sample_ids: str = "S1\n"):
    vcf = tmp_path / "variants.vcf.gz"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(sample_ids.splitlines())
        + "\n"
    )
    vat = tmp_path / "annotations.tsv.bgz"
    vat.write_text("")
    schema = tmp_path / "vat-schema.json"
    VatSchema.from_header(VAT_HEADER.split("\t")).write_json(schema)
    features = tmp_path / "features.tsv"
    features.write_text("chrom\ttss\tfeature_id\n" + features_text)
    samples = tmp_path / "shared-samples.txt"
    samples.write_text(sample_ids)
    return vcf, vat, schema, features, samples


@requires_htslib
def test_classification_is_gene_matched_chunked_and_annotation_aware(
    prepared_fixture, tmp_path: Path
):
    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc = tmp_path / "qc.json"
    classify_chromosome(
        prepared_fixture.vcf_gz,
        prepared_fixture.vat_bgz,
        prepared_fixture.vat_schema,
        prepared_fixture.features,
        prepared_fixture.samples,
        "chr1",
        [1, 2],
        [1, 2],
        ["stop_gained", "frameshift_variant", "missense_variant"],
        0.01,
        100,
        25,
        carriers,
        regions,
        qc,
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


def test_classification_streams_each_chunk_once_and_reduces_each_annotation_independently(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    vcf, vat, schema, features, samples = _write_inputs(
        tmp_path,
        "chr1\t100\tENSG000001.1\n"
        "chr1\t150\tENSG000002.2\n"
        "chr1\t300\tENSG000003.3\n"
        "chr1\t100\tENSG000004.4\n",
        "S1\nS2\nS3\n",
    )
    vat_rows = Path("tests/fixtures/transcript_annotations.tsv").read_text().splitlines()[1:]
    vat_rows.append(vat_rows[0])
    vat_rows.append(
        "chr1\t120\tA\tG\trs120\tENSG000001.9\tGENE1\tENST000120\ttrue\t"
        "missense_variant\tp.Ala2Val\t0.50\tLC\t.\t.\t.\t0.002\tafr"
    )
    vcf_rows = [
        "chr1\t100\t.\tA\tC\t.\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0",
        "chr1\t100\t.\tA\tG\t.\tPASS\tAC=1\tGT\t0/0\t0/0\t0/1",
        "chr1\t120\t.\tA\tG\t.\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0",
        "chr1\t150\t.\tG\tT\t.\tPASS\tAC=2\tGT\t0/0\t1/1\t0/0",
        "chr1\t180\t.\tC\tG\t.\tPASS\tAC=2\tGT\t0/0\t1/1\t0/0",
        "chr1\t300\t.\tT\tA\t.\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0",
    ]
    router = _TabixRouter({vcf: vcf_rows, vat: vat_rows})
    monkeypatch.setattr("rare_variant_enrichment.variants.subprocess.Popen", router)
    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc_path = tmp_path / "qc.json"

    classify_chromosome(
        vcf,
        vat,
        schema,
        features,
        samples,
        "chr1",
        [1, 2],
        [1, 2],
        ["stop_gained", "frameshift_variant", "missense_variant"],
        0.01,
        100,
        25,
        carriers,
        regions,
        qc_path,
    )

    rows = carriers.read_text().splitlines()
    assert "S1\tENSG000001.1\tAC=1\tbaseline\tall_rare_variants\t0" in rows
    assert "S1\tENSG000001.1\tAC=1\tconsequence\tstop_gained\t0" in rows
    assert "S1\tENSG000001.1\tAC=1\tconsequence\tmissense_variant\t20" in rows
    assert "S1\tENSG000001.1\tAC=1\tloftee\tHC\t0" in rows
    assert "S1\tENSG000001.1\tAC=1\tloftee\tLC\t20" in rows
    assert "S2\tENSG000002.2\tAC=2\tbaseline\tall_rare_variants\t30" in rows
    assert not any("frameshift_variant" in row and "ENSG000001.1" in row for row in rows)
    assert not any(row.startswith("S3\t") for row in rows)
    assert not any("ENSG000003.3" in row for row in rows)

    expected_regions = ["#chrom\tstart\tend"] + [
        f"chr1\t{start - 1}\t{min(400, start + 24)}" for start in range(1, 401, 25)
    ]
    assert regions.read_text().splitlines() == expected_regions
    summary = json.loads(qc_path.read_text())
    assert summary["annotation_chunk_count"] == 16
    assert summary["vat_tabix_query_count"] == 16
    assert summary["vcf_tabix_query_count"] == 16
    assert len(router.calls) == 32
    assert summary["vat_rows"] == 8
    assert summary["duplicate_vat_rows"] == 1
    assert summary["unique_vat_alleles"] == 5
    assert summary["vat_joined_alt_alleles"] == 5
    assert summary["vat_unmatched_alt_alleles"] == 1
    assert summary["above_maf_threshold_alleles"] == 1
    assert summary["missing_frequency_alleles"] == 1
    assert summary["converted_gvs_max_af_values"] == 1
    assert summary["observed_raw_gvs_max_af"] == 0.999
    assert summary["gene_matched_variant_feature_pairs"] == 3
    assert summary["gene_unmatched_variant_feature_pairs"] == 6
    assert summary["baseline_emitted_keys"] == 15
    assert summary["consequence_emitted_keys"] == 6
    assert summary["loftee_emitted_keys"] == 6
    assert summary["emitted_keys"] == 27
    qc_text = qc_path.read_text()
    assert "ENSG" not in qc_text
    assert "rs100" not in qc_text

    gathered_qc = tmp_path / "gathered_qc.tsv"
    gather_outputs(
        [carriers],
        [qc_path],
        tmp_path / "gathered_carriers.tsv",
        gathered_qc,
    )
    with gathered_qc.open() as handle:
        gathered_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert gathered_rows[0]["duplicate_vat_rows"] == "1"


def test_classify_chromosome_without_features_writes_zero_qc_without_tabix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    vcf, vat, schema, features, samples = _write_inputs(
        tmp_path, "chr2\t100\tENSG000002.2\n"
    )
    monkeypatch.setattr(
        "rare_variant_enrichment.variants.subprocess.Popen",
        lambda *_args, **_kwargs: pytest.fail("tabix must not run without selected features"),
    )
    carriers = tmp_path / "carriers.tsv"
    regions = tmp_path / "regions.bed"
    qc_path = tmp_path / "qc.json"

    classify_chromosome(
        vcf, vat, schema, features, samples, "chr1", [1], [1], [], 0.01, 100, 25,
        carriers, regions, qc_path,
    )

    assert carriers.read_text() == (
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\t"
        "minimum_distance_bp\n"
    )
    assert regions.read_text() == "#chrom\tstart\tend\n"
    summary = json.loads(qc_path.read_text())
    assert summary["chromosome"] == "chr1"
    assert summary["feature_count"] == 0
    assert summary["observed_raw_gvs_max_af"] == "none"
    assert all(
        value == 0
        for key, value in summary.items()
        if key not in {"chromosome", "observed_raw_gvs_max_af"}
    )


def test_classify_chromosome_rejects_requested_vcf_chromosome_before_tabix(tmp_path: Path):
    vcf, vat, schema, features, samples = _write_inputs(
        tmp_path, "chr2\t100\tENSG000002.2\n"
    )
    with pytest.raises(ValueError, match="Requested chromosome is absent from VCF: chr2"):
        classify_chromosome(
            vcf, vat, schema, features, samples, "chr2", [1], [1], [], 0.01, 100, 25,
            tmp_path / "carriers.tsv", tmp_path / "regions.bed", tmp_path / "qc.json",
        )


def test_chromosome_qc_preserves_ac_sources_call_states_and_boundary_pairs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    vcf, vat, schema, features, samples = _write_inputs(
        tmp_path, "chr1\t100\tENSG000001.1\n", "S1\nS2\nS3\n"
    )
    vcf_rows = [
        "chr1\t110\t.\tA\tC,G\t.\tPASS\tAC=.,2\tGT\t1/.\t0/2\t0/2"
    ]
    vat_rows = [
        "chr1\t110\tA\tC\t.\tENSG000001.1\tGENE1\tENST1\ttrue\tmissense_variant\t.\t.\t.\t.\t.\t.\t0.001\tglobal",
        "chr1\t110\tA\tG\t.\tENSG000001.1\tGENE1\tENST2\ttrue\tmissense_variant\t.\t.\t.\t.\t.\t.\t0.001\tglobal",
    ]
    router = _TabixRouter({vcf: vcf_rows, vat: vat_rows})
    monkeypatch.setattr("rare_variant_enrichment.variants.subprocess.Popen", router)
    qc_path = tmp_path / "qc.json"

    classify_chromosome(
        vcf, vat, schema, features, samples, "chr1", [1, 2], [], ["missense_variant"],
        0.01, 10, 100, tmp_path / "carriers.tsv", tmp_path / "regions.bed", qc_path,
    )

    qc = json.loads(qc_path.read_text())
    assert qc["info_ac_alt_alleles"] == 1
    assert qc["genotype_ac_fallback_alt_alleles"] == 1
    assert qc["unavailable_ac_alt_alleles"] == 0
    assert qc["partial_genotype_calls"] == 1
    assert qc["fully_missing_genotype_calls"] == 0
    assert qc["boundary_variant_feature_pairs"] == 2


def test_classify_chromosome_rejects_vcf_without_contig_metadata(tmp_path: Path):
    vcf, vat, schema, features, samples = _write_inputs(
        tmp_path, "chr1\t100\tENSG000001.1\n"
    )
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    with pytest.raises(ValueError, match="VCF header does not declare contigs"):
        classify_chromosome(
            vcf, vat, schema, features, samples, "chr1", [1], [], [], 0.01, 100, 25,
            tmp_path / "carriers.tsv", tmp_path / "regions.bed", tmp_path / "qc.json",
        )


def test_committed_vcf_fixture_declares_fields_and_has_concordant_info_ac():
    lines = Path("tests/fixtures/rare_variants.vcf").read_text().splitlines()
    assert any(line.startswith("##INFO=<ID=AC,Number=A,Type=Integer") for line in lines)
    assert any(line.startswith("##FORMAT=<ID=GT,Number=1,Type=String") for line in lines)
    sample_ids = next(line for line in lines if line.startswith("#CHROM")).split("\t")[9:]
    qc: dict[str, int] = {}
    for line in lines:
        if line.startswith("#"):
            continue
        parse_variant_alleles(line.split("\t"), sample_ids, set(sample_ids), qc=qc)
    assert qc["info_genotype_ac_mismatch_alt_alleles"] == 0
