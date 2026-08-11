"""Preparation and tabix-index validation for variant annotation tables."""

from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Sequence

from rare_variant_enrichment.annotations import VatSchema
from rare_variant_enrichment.io import open_text


def prepare_vat(
    vat_path: Path,
    chromosomes: Sequence[str],
    schema_output: Path,
) -> VatSchema:
    """Validate a bgzipped, coordinate-sorted VAT and emit its schema manifest."""
    try:
        with open_text(vat_path) as handle:
            header_line = next(line for line in handle if line.strip())
    except StopIteration as error:
        raise ValueError("VAT TSV is empty") from error
    except OSError as error:
        raise ValueError(f"Unable to read VAT TSV: {vat_path}") from error

    schema = VatSchema.from_header(header_line.rstrip("\r\n").split("\t"))
    index_path = Path(f"{vat_path}.tbi")
    if not index_path.is_file():
        _create_index(vat_path, schema)
    contigs = _listed_contigs(vat_path)
    missing = [chromosome for chromosome in chromosomes if chromosome not in contigs]
    if missing:
        raise ValueError("Requested chromosomes are absent from VAT index: " + ", ".join(missing))
    schema.write_json(schema_output)
    return schema


def _create_index(vat_path: Path, schema: VatSchema) -> None:
    command = [
        "tabix",
        "-f",
        "-S",
        "1",
        "-s",
        str(schema.chromosome + 1),
        "-b",
        str(schema.position + 1),
        "-e",
        str(schema.position + 1),
        str(vat_path),
    ]
    try:
        subprocess.run(command, check=True)
    except OSError as error:
        raise ValueError(
            f"Unable to prepare VAT index; tabix/bgzip support is required: {vat_path}"
        ) from error
    except subprocess.CalledProcessError as error:
        raise ValueError(
            "Unable to prepare VAT index; VAT must be bgzip-compressed and coordinate-sorted: "
            f"{vat_path}"
        ) from error


def _listed_contigs(vat_path: Path) -> set[str]:
    command = ["tabix", "-l", str(vat_path)]
    try:
        result = subprocess.run(command, text=True, capture_output=True, check=True)
    except OSError as error:
        raise ValueError(
            f"Unable to validate VAT index; tabix/bgzip support is required: {vat_path}"
        ) from error
    except subprocess.CalledProcessError as error:
        raise ValueError(f"Unable to validate VAT index: {vat_path}") from error
    return set(result.stdout.splitlines())
