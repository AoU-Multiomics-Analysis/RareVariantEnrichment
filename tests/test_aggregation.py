import json
from pathlib import Path

import pytest

from rare_variant_enrichment.aggregation import gather_outputs


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


def test_gather_keeps_global_minimum_distance(tmp_path: Path):
    first = tmp_path / "chr1.tsv"
    second = tmp_path / "chr2.tsv"
    first.write_text(
        carrier_header
        + "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t100\n"
        + "S2\tGENE2\tAC<=3\tbaseline\tall_rare_variants\t20\n"
    )
    second.write_text(carrier_header + "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t40\n")
    q1 = tmp_path / "chr1.json"
    q2 = tmp_path / "chr2.json"
    q1.write_text('{"chromosome":"chr1","vcf_records_extracted":2}')
    q2.write_text('{"chromosome":"chr2","vcf_records_extracted":1}')

    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "chromosome_qc.tsv"
    gather_outputs([first, second], [q1, q2], carrier_output, qc_output)

    assert "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t40" in carrier_output.read_text().splitlines()
    assert qc_output.read_text().splitlines()[1].startswith("chr1\t")


def test_gather_writes_carriers_and_qc_in_deterministic_order(tmp_path: Path):
    first = tmp_path / "chr2.tsv"
    first.write_text(
        carrier_header
        + "S2\tGENE2\tAC=2\tbaseline\tall_rare_variants\t8\n"
        + "S1\tGENE1\tAC=2\tbaseline\tall_rare_variants\t9\n"
    )
    second = tmp_path / "chr1.tsv"
    second.write_text(
        carrier_header + "S3\tGENE1\tAC=1\tbaseline\tall_rare_variants\t7\n"
    )
    q1 = tmp_path / "chr2.json"
    q1.write_text('{"chromosome":"chr2","z_value":2}')
    q2 = tmp_path / "chr1.json"
    q2.write_text('{"vcf_records_extracted":1,"chromosome":"chr1"}')

    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "chromosome_qc.tsv"
    gather_outputs([first, second], [q1, q2], carrier_output, qc_output)

    assert carrier_output.read_text().splitlines() == [
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp",
        "S1\tGENE1\tAC=2\tbaseline\tall_rare_variants\t9",
        "S3\tGENE1\tAC=1\tbaseline\tall_rare_variants\t7",
        "S2\tGENE2\tAC=2\tbaseline\tall_rare_variants\t8",
    ]
    assert qc_output.read_text().splitlines() == [
        "chromosome\tvcf_records_extracted\tz_value",
        "chr1\t1\t",
        "chr2\t\t2",
    ]


def test_gather_rejects_mismatched_carrier_headers(tmp_path: Path):
    first = tmp_path / "chr1.tsv"
    first.write_text(carrier_header)
    second = tmp_path / "chr2.tsv"
    second.write_text("sample_id\tfeature_id\tac_class\tdistance\n")
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1"}')

    with pytest.raises(ValueError, match="Carrier TSV header does not match"):
        gather_outputs([first, second], [qc, qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")


def test_gather_rejects_malformed_carrier_distance(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(
        carrier_header + "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\tnot-an-int\n"
    )
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1"}')

    with pytest.raises(ValueError, match="minimum_distance_bp must be a non-negative integer"):
        gather_outputs([carrier], [qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")


def test_gather_rejects_mismatched_input_counts(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(carrier_header)
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1"}')

    with pytest.raises(ValueError, match="same number of carrier and QC inputs"):
        gather_outputs([carrier], [], tmp_path / "all.tsv", tmp_path / "qc.tsv")


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ('{"vcf_records_extracted":1}', "QC object must contain one non-empty chromosome"),
        ('{"chromosome":" "}', "QC object must contain one non-empty chromosome"),
        ('{"chromosome":"chr1","chromosome":"chr1"}', "QC JSON contains duplicate key"),
    ],
)
def test_gather_rejects_qc_without_one_valid_chromosome(
    tmp_path: Path, payload: str, message: str
):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(carrier_header)
    qc = tmp_path / "chr1.json"
    qc.write_text(payload)

    with pytest.raises(ValueError, match=message):
        gather_outputs([carrier], [qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")


def test_gather_deduplicates_within_but_not_across_annotation_classes(tmp_path: Path):
    first = write_carriers(
        tmp_path / "a.tsv",
        [
            ("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 50),
            ("S1", "ENSG1.1", "AC=1", "consequence", "missense_variant", 20),
        ],
    )
    second = write_carriers(
        tmp_path / "b.tsv",
        [("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 10)],
    )
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


def test_gather_rejects_run_with_no_vat_allele_key_matches(tmp_path: Path):
    """Changing the all-unmatched join guard must fail this run-level contract."""
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(carrier_header)
    empty_feature_carrier = tmp_path / "chr2.tsv"
    empty_feature_carrier.write_text(carrier_header)
    first_qc = tmp_path / "chr1.json"
    first_qc.write_text(
        '{"chromosome":"chr1","classified_alt_alleles":4,"vat_joined_alt_alleles":0}'
    )
    empty_feature_qc = tmp_path / "chr2.json"
    empty_feature_qc.write_text(
        '{"chromosome":"chr2","classified_alt_alleles":0,"vat_joined_alt_alleles":0}'
    )

    with pytest.raises(ValueError, match="No queried VCF ALT alleles matched VAT allele keys"):
        gather_outputs(
            [carrier, empty_feature_carrier],
            [first_qc, empty_feature_qc],
            tmp_path / "all.tsv",
            tmp_path / "qc.tsv",
        )


def test_gather_allows_empty_feature_chromosomes_when_no_eligible_alt_alleles_exist(
    tmp_path: Path,
):
    """Changing a zero eligible-allele run into a failure must fail this contract."""
    first = tmp_path / "chr1.tsv"
    second = tmp_path / "chr2.tsv"
    first.write_text(carrier_header)
    second.write_text(carrier_header)
    first_qc = tmp_path / "chr1.json"
    second_qc = tmp_path / "chr2.json"
    first_qc.write_text(
        '{"chromosome":"chr1","classified_alt_alleles":0,"vat_joined_alt_alleles":0}'
    )
    second_qc.write_text(
        '{"chromosome":"chr2","classified_alt_alleles":0,"vat_joined_alt_alleles":0}'
    )
    qc_output = tmp_path / "qc.tsv"

    gather_outputs([first, second], [first_qc, second_qc], tmp_path / "all.tsv", qc_output)

    assert qc_output.read_text().splitlines() == [
        "chromosome\tclassified_alt_alleles\tvat_joined_alt_alleles",
        "chr1\t0\t0",
        "chr2\t0\t0",
    ]


@pytest.mark.parametrize(
    ("counter", "value"),
    [
        ("classified_alt_alleles", 0.5),
        ("vat_joined_alt_alleles", 1.0),
        ("classified_alt_alleles", True),
        ("vat_joined_alt_alleles", "1"),
        ("classified_alt_alleles", -1),
    ],
)
def test_gather_rejects_non_integer_zero_match_counters(
    tmp_path: Path, counter: str, value: object
):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(carrier_header)
    qc = tmp_path / "chr1.json"
    payload = {
        "chromosome": "chr1",
        "classified_alt_alleles": 1,
        "vat_joined_alt_alleles": 1,
        counter: value,
    }
    qc.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match=rf"QC {counter} must be a non-negative integer"):
        gather_outputs([carrier], [qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")
