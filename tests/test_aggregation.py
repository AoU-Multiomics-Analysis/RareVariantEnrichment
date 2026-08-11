from pathlib import Path

import pytest

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


def test_gather_writes_carriers_and_qc_in_deterministic_order(tmp_path: Path):
    header = "sample_id\tfeature_id\tac_class\tminimum_distance_bp\n"
    first = tmp_path / "chr2.tsv"
    first.write_text(header + "S2\tGENE2\tAC=2\t8\nS1\tGENE1\tAC=2\t9\n")
    second = tmp_path / "chr1.tsv"
    second.write_text(header + "S3\tGENE1\tAC=1\t7\n")
    q1 = tmp_path / "chr2.json"
    q1.write_text('{"chromosome":"chr2","z_value":2}')
    q2 = tmp_path / "chr1.json"
    q2.write_text('{"vcf_records_extracted":1,"chromosome":"chr1"}')

    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "chromosome_qc.tsv"
    gather_outputs([first, second], [q1, q2], carrier_output, qc_output)

    assert carrier_output.read_text().splitlines() == [
        "sample_id\tfeature_id\tac_class\tminimum_distance_bp",
        "S1\tGENE1\tAC=2\t9",
        "S3\tGENE1\tAC=1\t7",
        "S2\tGENE2\tAC=2\t8",
    ]
    assert qc_output.read_text().splitlines() == [
        "chromosome\tvcf_records_extracted\tz_value",
        "chr1\t1\t",
        "chr2\t\t2",
    ]


def test_gather_rejects_mismatched_carrier_headers(tmp_path: Path):
    first = tmp_path / "chr1.tsv"
    first.write_text("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
    second = tmp_path / "chr2.tsv"
    second.write_text("sample_id\tfeature_id\tac_class\tdistance\n")
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1"}')

    with pytest.raises(ValueError, match="Carrier TSV header does not match"):
        gather_outputs([first, second], [qc, qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")


def test_gather_rejects_malformed_carrier_distance(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(
        "sample_id\tfeature_id\tac_class\tminimum_distance_bp\nS1\tGENE1\tAC=1\tnot-an-int\n"
    )
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1"}')

    with pytest.raises(ValueError, match="minimum_distance_bp must be a non-negative integer"):
        gather_outputs([carrier], [qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")


def test_gather_rejects_mismatched_input_counts(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
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
    carrier.write_text("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
    qc = tmp_path / "chr1.json"
    qc.write_text(payload)

    with pytest.raises(ValueError, match=message):
        gather_outputs([carrier], [qc], tmp_path / "all.tsv", tmp_path / "qc.tsv")
