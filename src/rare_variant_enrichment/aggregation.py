import json
import math
from pathlib import Path
from typing import Sequence

from rare_variant_enrichment.io import open_text


CARRIER_HEADER = ("sample_id", "feature_id", "ac_class", "minimum_distance_bp")


def gather_outputs(
    carrier_paths: Sequence[Path],
    qc_paths: Sequence[Path],
    carrier_output: Path,
    qc_output: Path,
) -> None:
    """Combine per-chromosome outputs into deterministic carrier and QC tables."""
    if len(carrier_paths) != len(qc_paths):
        raise ValueError("Expected the same number of carrier and QC inputs")
    if not carrier_paths:
        raise ValueError("At least one carrier input and one QC input are required")

    minimum_distances: dict[tuple[str, str, str], int] = {}
    for path in carrier_paths:
        _read_carrier_file(path, minimum_distances)

    qc_records = [_read_qc_file(path) for path in qc_paths]
    chromosomes = [record["chromosome"] for record in qc_records]
    if len(chromosomes) != len(set(chromosomes)):
        raise ValueError("QC inputs must contain unique chromosomes")

    _write_carriers(carrier_output, minimum_distances)
    _write_qc(qc_output, qc_records)


def _read_carrier_file(
    path: Path, minimum_distances: dict[tuple[str, str, str], int]
) -> None:
    with open_text(path) as handle:
        try:
            header = next(handle).rstrip("\r\n").split("\t")
        except StopIteration as error:
            raise ValueError(f"Carrier TSV is empty: {path}") from error
        if tuple(header) != CARRIER_HEADER:
            raise ValueError(f"Carrier TSV header does not match expected header: {path}")

        for line_number, raw_line in enumerate(handle, start=2):
            if not raw_line.strip():
                raise ValueError(f"Carrier TSV line {line_number} is blank: {path}")
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != len(CARRIER_HEADER):
                raise ValueError(f"Carrier TSV line {line_number} must have four columns: {path}")
            sample_id, feature_id, ac_class, distance_text = fields
            if not sample_id or not feature_id or not ac_class:
                raise ValueError(f"Carrier TSV line {line_number} has empty key fields: {path}")
            try:
                distance = int(distance_text)
            except ValueError as error:
                raise ValueError(
                    f"Carrier TSV line {line_number} minimum_distance_bp must be a non-negative integer: {path}"
                ) from error
            if distance < 0:
                raise ValueError(
                    f"Carrier TSV line {line_number} minimum_distance_bp must be a non-negative integer: {path}"
                )

            key = (sample_id, feature_id, ac_class)
            previous_distance = minimum_distances.get(key)
            if previous_distance is None or distance < previous_distance:
                minimum_distances[key] = distance


def _read_qc_file(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_unique_json_object)
    except json.JSONDecodeError as error:
        raise ValueError(f"QC JSON is invalid: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"QC JSON must contain an object: {path}")

    chromosome = payload.get("chromosome")
    if not isinstance(chromosome, str) or not chromosome.strip():
        raise ValueError(f"QC object must contain one non-empty chromosome: {path}")
    for key, value in payload.items():
        _validate_qc_cell(key, value, path)
    return payload


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    payload: dict[str, object] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError(f"QC JSON contains duplicate key: {key}")
        payload[key] = value
    return payload


def _validate_qc_cell(key: str, value: object, path: Path) -> None:
    if not key or any(character in key for character in "\t\r\n"):
        raise ValueError(f"QC JSON key is not valid for TSV output: {path}")
    if isinstance(value, str):
        if any(character in value for character in "\t\r\n"):
            raise ValueError(f"QC JSON value is not valid for TSV output: {path}")
        return
    if isinstance(value, bool) or value is None or isinstance(value, int):
        return
    if isinstance(value, float) and math.isfinite(value):
        return
    raise ValueError(f"QC JSON values must be scalar TSV values: {path}")


def _write_carriers(
    path: Path, minimum_distances: dict[tuple[str, str, str], int]
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(CARRIER_HEADER) + "\n")
        for (sample_id, feature_id, ac_class), distance in sorted(
            minimum_distances.items(), key=lambda item: (item[0][1], item[0][0], item[0][2])
        ):
            handle.write(f"{sample_id}\t{feature_id}\t{ac_class}\t{distance}\n")


def _write_qc(path: Path, qc_records: Sequence[dict[str, object]]) -> None:
    columns = ["chromosome", *sorted({key for record in qc_records for key in record} - {"chromosome"})]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for record in sorted(qc_records, key=lambda item: str(item["chromosome"])):
            handle.write("\t".join(_format_qc_value(record.get(column)) for column in columns) + "\n")


def _format_qc_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)
