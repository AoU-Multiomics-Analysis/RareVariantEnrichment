import csv
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from rare_variant_enrichment.artifacts import file_artifact
from rare_variant_enrichment.carrier_definitions import (
    CARRIER_DEFINITION_HEADER,
    CarrierDefinition,
    build_carrier_definitions,
    read_carrier_definition_config,
)
from rare_variant_enrichment.carrier_extraction import AUDIT_HEADER


def _write_config(tmp_path: Path, payload: object) -> Path:
    path = tmp_path / "definitions.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _audit_row(
    variant_id: str,
    variant_classes: str,
    *,
    sample_id: str = "S1",
    gene_id: str = "ENSG1",
    gene_symbol: str = "GENE1",
    revel: str = "",
) -> dict[str, str]:
    chrom, pos, ref, alt = variant_id.split(":")
    return {
        "sample_id": sample_id,
        "gene_id": gene_id,
        "gene_symbol": gene_symbol,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "variant_id": variant_id,
        "variant_ac": "1",
        "variant_af": "0.001",
        "sample_alt_allele_count": "1",
        "most_severe_consequence": "missense_variant",
        "all_consequences": "missense_variant",
        "unknown_consequences": "",
        "loftee": "",
        "revel": revel,
        "gvs_max_af": "0.001",
        "variant_classes": variant_classes,
    }


def _write_audit_with_qc(
    tmp_path: Path, rows: list[dict[str, str]]
) -> tuple[Path, Path]:
    audit = tmp_path / "variant_carrier_audit.tsv.gz"
    with gzip.open(audit, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    qc = tmp_path / "variant_carriers.qc.json"
    qc.write_text(
        json.dumps(
            {
                "audit_artifact": file_artifact(
                    audit,
                    "variant_carrier_audit.tsv.gz",
                    AUDIT_HEADER,
                    len(rows),
                ),
                "vcf_index_provenance": "supplied",
                "transcript_index_provenance": "generated",
                "quality_or_frequency_filters_applied": False,
            }
        ),
        encoding="utf-8",
    )
    return audit, qc


def _read_carrier_rows(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_definition_config_preserves_order_and_threshold(tmp_path: Path):
    config = _write_config(
        tmp_path,
        {
            "schema_version": 1,
            "definitions": [
                {"name": "lof_hc", "variant_classes": ["lof_hc"]},
                {
                    "name": "missense_revel_ge_0_75",
                    "variant_classes": ["missense"],
                    "minimum_revel": 0.75,
                },
            ],
        },
    )

    parsed = read_carrier_definition_config(config)

    assert parsed.schema_version == 1
    assert parsed.names == ("lof_hc", "missense_revel_ge_0_75")
    assert parsed.definitions[0].minimum_revel is None
    assert parsed.definitions[1].minimum_revel == 0.75


def test_definition_matches_class_union_and_revel_threshold():
    definition = CarrierDefinition(
        name="coding_revel",
        variant_classes=("lof_hc", "missense"),
        minimum_revel=0.75,
    )

    assert definition.matches(frozenset({"missense"}), 0.75)
    assert definition.matches(frozenset({"lof_hc"}), 0.81)
    assert not definition.matches(frozenset({"splice_core"}), 0.9)
    assert not definition.matches(frozenset({"missense"}), 0.74)
    assert not definition.matches(frozenset({"missense"}), None)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ({"schema_version": True, "definitions": []}, "schema_version"),
        ({"schema_version": 2, "definitions": []}, "Unsupported"),
        ({"schema_version": 1, "definitions": []}, "at least one"),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "1bad", "variant_classes": ["missense"]}
                ],
            },
            "name",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "bad", "variant_classes": ["unknown"]}
                ],
            },
            "variant class",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {
                        "name": "bad",
                        "variant_classes": ["missense"],
                        "minimum_revel": 1.1,
                    }
                ],
            },
            "minimum_revel",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {
                        "name": "bad",
                        "variant_classes": ["missense"],
                        "minimum_revel": True,
                    }
                ],
            },
            "minimum_revel",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "bad", "variant_classes": ["missense", "missense"]}
                ],
            },
            "duplicate variant class",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "same", "variant_classes": ["missense"]},
                    {"name": "same", "variant_classes": ["lof_hc"]},
                ],
            },
            "duplicate definition name",
        ),
        (
            {"schema_version": 1, "definitions": [], "extra": 1},
            "unknown top-level key",
        ),
        (
            {
                "schema_version": 1,
                "definitions": [
                    {"name": "bad", "variant_classes": ["missense"], "extra": 1}
                ],
            },
            "unknown definition key",
        ),
    ],
)
def test_definition_config_rejects_invalid_values(
    tmp_path: Path, payload: object, message: str
):
    config = _write_config(tmp_path, payload)

    with pytest.raises(ValueError, match=message):
        read_carrier_definition_config(config)


def test_definition_config_rejects_duplicate_json_keys(tmp_path: Path):
    config = tmp_path / "definitions.json"
    config.write_text(
        '{"schema_version":1,"schema_version":1,"definitions":[]}',
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate JSON key"):
        read_carrier_definition_config(config)


def test_build_carrier_definitions_materializes_overlapping_rules(tmp_path: Path):
    hc = _audit_row("chr1:100:A:C", "lof_hc,lof_hc_or_lc")
    rows = [
        hc,
        _audit_row("chr1:110:A:G", "missense", revel="0.74"),
        _audit_row("chr1:120:C:T", "missense", revel="0.75"),
        _audit_row("chr1:130:G:A", "missense"),
        _audit_row("chr1:140:T:C", "splice_core"),
        _audit_row("chr1:150:A:T", "splice_region"),
        hc.copy(),
    ]
    audit, extraction_qc = _write_audit_with_qc(tmp_path, rows)
    config = _write_config(
        tmp_path,
        {
            "schema_version": 1,
            "definitions": [
                {"name": "lof_hc", "variant_classes": ["lof_hc"]},
                {"name": "missense", "variant_classes": ["missense"]},
                {
                    "name": "missense_revel_ge_0_75",
                    "variant_classes": ["missense"],
                    "minimum_revel": 0.75,
                },
                {
                    "name": "splice_any",
                    "variant_classes": ["splice_core", "splice_region"],
                },
            ],
        },
    )
    output = tmp_path / "carrier_definitions.tsv.gz"
    qc_output = tmp_path / "carrier_definitions.qc.json"

    build_carrier_definitions(
        audit,
        extraction_qc,
        config,
        output,
        qc_output,
        container_image="example.invalid/enrichment@sha256:abc",
    )

    carrier_rows = _read_carrier_rows(output)
    assert list(carrier_rows[0]) == list(CARRIER_DEFINITION_HEADER)
    assert [
        (
            row["sample_id"],
            row["gene_id"],
            row["carrier_definition"],
            row["n_variants"],
            row["variant_ids"],
        )
        for row in carrier_rows
    ] == [
        ("S1", "ENSG1", "lof_hc", "1", "chr1:100:A:C"),
        (
            "S1",
            "ENSG1",
            "missense",
            "3",
            "chr1:110:A:G,chr1:120:C:T,chr1:130:G:A",
        ),
        (
            "S1",
            "ENSG1",
            "missense_revel_ge_0_75",
            "1",
            "chr1:120:C:T",
        ),
        (
            "S1",
            "ENSG1",
            "splice_any",
            "2",
            "chr1:140:T:C,chr1:150:A:T",
        ),
    ]

    manifest = json.loads(qc_output.read_text())
    assert manifest["schema"] == "aou.carrier-definitions-manifest.v1"
    assert manifest["definition_order"] == [
        "lof_hc",
        "missense",
        "missense_revel_ge_0_75",
        "splice_any",
    ]
    assert manifest["audit_counts"] == {
        "input_rows": 7,
        "deduplicated_rows": 6,
        "duplicate_rows": 1,
        "present_revel_rows": 2,
        "missing_revel_rows": 4,
    }
    assert manifest["definition_counts"]["missense"] == {
        "matched_audit_rows": 3,
        "distinct_variants": 3,
        "distinct_sample_gene_pairs": 1,
        "output_rows": 1,
    }
    assert manifest["output_artifact"] == file_artifact(
        output,
        "carrier_definitions.tsv.gz",
        CARRIER_DEFINITION_HEADER,
        4,
    )
    assert manifest["provenance"]["container_image"] == (
        "example.invalid/enrichment@sha256:abc"
    )


def test_build_carrier_definitions_rejects_audit_digest_mismatch(tmp_path: Path):
    audit, extraction_qc = _write_audit_with_qc(
        tmp_path, [_audit_row("chr1:100:A:C", "missense", revel="0.8")]
    )
    config = _write_config(
        tmp_path,
        {
            "schema_version": 1,
            "definitions": [
                {"name": "missense", "variant_classes": ["missense"]}
            ],
        },
    )
    audit.write_bytes(audit.read_bytes() + b"x")

    with pytest.raises(ValueError, match="SHA-256"):
        build_carrier_definitions(
            audit,
            extraction_qc,
            config,
            tmp_path / "output.tsv.gz",
            tmp_path / "output.json",
            container_image="image@sha256:abc",
        )


def test_build_carrier_definitions_is_byte_deterministic(tmp_path: Path):
    audit, extraction_qc = _write_audit_with_qc(
        tmp_path, [_audit_row("chr1:100:A:C", "missense", revel="0.8")]
    )
    config = _write_config(
        tmp_path,
        {
            "schema_version": 1,
            "definitions": [
                {"name": "missense", "variant_classes": ["missense"]}
            ],
        },
    )
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"

    build_carrier_definitions(
        audit,
        extraction_qc,
        config,
        first,
        tmp_path / "first.json",
        container_image="image@sha256:abc",
    )
    build_carrier_definitions(
        audit,
        extraction_qc,
        config,
        second,
        tmp_path / "second.json",
        container_image="image@sha256:abc",
    )

    assert first.read_bytes() == second.read_bytes()


def test_build_carrier_definitions_writes_header_and_zero_counts(tmp_path: Path):
    audit, extraction_qc = _write_audit_with_qc(
        tmp_path, [_audit_row("chr1:100:A:C", "missense", revel="0.7")]
    )
    config = _write_config(
        tmp_path,
        {
            "schema_version": 1,
            "definitions": [
                {"name": "splice_core", "variant_classes": ["splice_core"]}
            ],
        },
    )
    output = tmp_path / "output.tsv.gz"
    qc_output = tmp_path / "output.json"

    build_carrier_definitions(
        audit,
        extraction_qc,
        config,
        output,
        qc_output,
        container_image="image@sha256:abc",
    )

    assert _read_carrier_rows(output) == []
    assert json.loads(qc_output.read_text())["definition_counts"]["splice_core"] == {
        "matched_audit_rows": 0,
        "distinct_variants": 0,
        "distinct_sample_gene_pairs": 0,
        "output_rows": 0,
    }
