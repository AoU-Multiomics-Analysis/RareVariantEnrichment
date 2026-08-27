import hashlib
from pathlib import Path

import pytest

from rare_variant_enrichment.artifacts import file_artifact


def test_file_artifact_records_exact_bytes_and_schema(tmp_path: Path):
    source = tmp_path / "audit.tsv.gz"
    source.write_bytes(b"abc")

    assert file_artifact(source, "variant_carrier_audit.tsv.gz", ("a", "b"), 3) == {
        "logical_name": "variant_carrier_audit.tsv.gz",
        "header": ["a", "b"],
        "row_count": 3,
        "size_bytes": 3,
        "sha256": hashlib.sha256(b"abc").hexdigest(),
    }


def test_file_artifact_rejects_negative_row_count(tmp_path: Path):
    source = tmp_path / "audit.tsv.gz"
    source.write_bytes(b"")

    with pytest.raises(ValueError, match="row_count"):
        file_artifact(source, "audit", ("a",), -1)
