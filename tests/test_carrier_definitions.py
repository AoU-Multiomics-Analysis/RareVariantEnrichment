import json
from pathlib import Path

import pytest

from rare_variant_enrichment.carrier_definitions import (
    CarrierDefinition,
    read_carrier_definition_config,
)


def _write_config(tmp_path: Path, payload: object) -> Path:
    path = tmp_path / "definitions.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


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
