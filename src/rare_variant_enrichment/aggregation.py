import json
import math
from pathlib import Path
from typing import Sequence

from rare_variant_enrichment.io import open_text
from rare_variant_enrichment.storage import MinimumDistanceStore


CARRIER_HEADER = (
    "sample_id",
    "feature_id",
    "ac_class",
    "annotation_family",
    "annotation_class",
    "minimum_distance_bp",
)


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

    with MinimumDistanceStore(carrier_output.parent) as minimum_distances:
        for path in carrier_paths:
            for (
                sample_id,
                feature_id,
                ac_class,
                annotation_family,
                annotation_class,
                distance,
            ) in _iter_carrier_file(path):
                minimum_distances.upsert(
                    sample_id,
                    feature_id,
                    ac_class,
                    annotation_family,
                    annotation_class,
                    distance,
                )
        qc_records = [_read_qc_file(path) for path in qc_paths]
        chromosomes = [record["chromosome"] for record in qc_records]
        if len(chromosomes) != len(set(chromosomes)):
            raise ValueError("QC inputs must contain unique chromosomes")
        classified_alt_alleles = _sum_qc_counter(qc_records, "classified_alt_alleles")
        vat_joined_alt_alleles = _sum_qc_counter(qc_records, "vat_joined_alt_alleles")
        if classified_alt_alleles > 0 and vat_joined_alt_alleles == 0:
            raise ValueError("No queried VCF ALT alleles matched VAT allele keys")
        minimum_distances.write_tsv(carrier_output, "feature")
    _write_qc(qc_output, qc_records)


def _iter_carrier_file(path: Path):
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
                raise ValueError(f"Carrier TSV line {line_number} must have six columns: {path}")
            (
                sample_id,
                feature_id,
                ac_class,
                annotation_family,
                annotation_class,
                distance_text,
            ) = fields
            if not all((sample_id, feature_id, ac_class, annotation_family, annotation_class)):
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

            yield sample_id, feature_id, ac_class, annotation_family, annotation_class, distance


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


def _sum_qc_counter(qc_records: Sequence[dict[str, object]], counter: str) -> int:
    total = 0
    for record in qc_records:
        value = record.get(counter, 0)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"QC {counter} must be a non-negative integer")
        total += value
    return total


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
