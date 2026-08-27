"""Validate and apply named carrier-definition rules."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import Any


SUPPORTED_BASE_CLASSES = frozenset(
    {"lof_hc", "lof_hc_or_lc", "missense", "splice_core", "splice_region"}
)
_NAME_PATTERN = re.compile(r"[A-Za-z][A-Za-z0-9_.-]*\Z")


@dataclass(frozen=True)
class CarrierDefinition:
    """One ordered carrier rule from a schema-versioned configuration."""

    name: str
    variant_classes: tuple[str, ...]
    minimum_revel: float | None = None

    def matches(self, classes: frozenset[str], revel: float | None) -> bool:
        class_match = any(value in classes for value in self.variant_classes)
        threshold_match = self.minimum_revel is None or (
            revel is not None and revel >= self.minimum_revel
        )
        return class_match and threshold_match


@dataclass(frozen=True)
class CarrierDefinitionConfig:
    """Validated carrier-definition configuration."""

    schema_version: int
    definitions: tuple[CarrierDefinition, ...]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.definitions)


def read_carrier_definition_config(path: Path) -> CarrierDefinitionConfig:
    """Read and strictly validate a schema-version 1 definition file."""
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, object_pairs_hook=_reject_duplicate_json_keys)
    except json.JSONDecodeError as error:
        raise ValueError(f"Carrier-definition JSON is invalid: {path}") from error

    if not isinstance(payload, dict):
        raise ValueError("Carrier-definition JSON must contain an object")
    unknown_top_level = set(payload) - {"schema_version", "definitions"}
    if unknown_top_level:
        raise ValueError(
            "Carrier-definition JSON has an unknown top-level key: "
            + ", ".join(sorted(unknown_top_level))
        )
    missing_top_level = {"schema_version", "definitions"} - set(payload)
    if missing_top_level:
        raise ValueError(
            "Carrier-definition JSON is missing a top-level key: "
            + ", ".join(sorted(missing_top_level))
        )

    schema_version = payload["schema_version"]
    if isinstance(schema_version, bool) or not isinstance(schema_version, int):
        raise ValueError("schema_version must be the integer 1")
    if schema_version != 1:
        raise ValueError(f"Unsupported carrier-definition schema_version: {schema_version}")

    raw_definitions = payload["definitions"]
    if not isinstance(raw_definitions, list):
        raise ValueError("definitions must be an array")
    if not raw_definitions:
        raise ValueError("definitions must contain at least one carrier definition")

    definitions: list[CarrierDefinition] = []
    names: set[str] = set()
    for index, raw_definition in enumerate(raw_definitions):
        definition = _parse_definition(raw_definition, index)
        if definition.name in names:
            raise ValueError(f"Carrier-definition JSON has duplicate definition name: {definition.name}")
        names.add(definition.name)
        definitions.append(definition)

    return CarrierDefinitionConfig(
        schema_version=schema_version,
        definitions=tuple(definitions),
    )


def _parse_definition(raw_definition: Any, index: int) -> CarrierDefinition:
    if not isinstance(raw_definition, dict):
        raise ValueError(f"Definition {index} must be an object")
    allowed_keys = {"name", "variant_classes", "minimum_revel"}
    unknown_keys = set(raw_definition) - allowed_keys
    if unknown_keys:
        raise ValueError(
            f"Definition {index} has an unknown definition key: "
            + ", ".join(sorted(unknown_keys))
        )
    missing_keys = {"name", "variant_classes"} - set(raw_definition)
    if missing_keys:
        raise ValueError(
            f"Definition {index} is missing a required key: "
            + ", ".join(sorted(missing_keys))
        )

    name = raw_definition["name"]
    if not isinstance(name, str) or _NAME_PATTERN.fullmatch(name) is None:
        raise ValueError(f"Definition {index} has an invalid name")

    raw_classes = raw_definition["variant_classes"]
    if not isinstance(raw_classes, list) or not raw_classes:
        raise ValueError(f"Definition {name} must have at least one variant class")
    if not all(isinstance(value, str) for value in raw_classes):
        raise ValueError(f"Definition {name} has a non-string variant class")
    if len(raw_classes) != len(set(raw_classes)):
        raise ValueError(f"Definition {name} has a duplicate variant class")
    invalid_classes = sorted(set(raw_classes) - SUPPORTED_BASE_CLASSES)
    if invalid_classes:
        raise ValueError(
            f"Definition {name} has an unsupported variant class: "
            + ", ".join(invalid_classes)
        )

    minimum_revel = raw_definition.get("minimum_revel")
    if minimum_revel is not None:
        if isinstance(minimum_revel, bool) or not isinstance(minimum_revel, (int, float)):
            raise ValueError(f"Definition {name} minimum_revel must be numeric")
        minimum_revel = float(minimum_revel)
        if not math.isfinite(minimum_revel) or not 0.0 <= minimum_revel <= 1.0:
            raise ValueError(f"Definition {name} minimum_revel must be from 0 through 1")

    return CarrierDefinition(
        name=name,
        variant_classes=tuple(raw_classes),
        minimum_revel=minimum_revel,
    )


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"Carrier-definition JSON contains duplicate JSON key: {key}")
        payload[key] = value
    return payload
