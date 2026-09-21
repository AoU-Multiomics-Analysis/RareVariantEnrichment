"""Validate ome manifest metadata before WDL declares its paths as File inputs."""

import csv
import logging
from pathlib import Path
from typing import Sequence
from urllib.parse import urlsplit

from rare_variant_enrichment.io import open_text
from rare_variant_enrichment.multiomics import require_local_file, validate_multiomics_inputs


LOGGER = logging.getLogger(__name__)
COLUMNS = ('ome_name', 'phenotype_bed', 'principal_components_tsv', 'additional_covariates_tsv')


def _manifest_path(value: str, column: str, line: int, optional: bool = False) -> str:
    if any(character in value for character in ('\n', '\r', '\t', '\0')):
        raise ValueError(f'Manifest line {line} has a control character in {column}')
    value = value.strip()
    if value in ('', '.'):
        if optional:
            return '.'
        raise ValueError(f'Manifest line {line} requires {column}')
    if value.startswith('gs://'):
        uri = urlsplit(value)
        if not uri.netloc or not uri.path.strip('/') or uri.query or uri.fragment:
            raise ValueError(f'Manifest line {line} requires a GCS bucket and object in {column}')
    elif not Path(value).is_absolute():
        raise ValueError(f'Manifest line {line}: {column} must be a gs:// URI or an absolute local path')
    return value


def prepare_omics_manifest(manifest: Path, thresholds: Sequence[float], output: Path) -> None:
    require_local_file(manifest)
    rows = []
    with open_text(manifest) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if len(set(header)) != len(header):
            raise ValueError('Ome manifest has duplicate header names')
        missing = set(COLUMNS[:3]) - set(header)
        if missing:
            raise ValueError('Ome manifest is missing columns: ' + ', '.join(sorted(missing)))
        unknown = set(header) - set(COLUMNS)
        if unknown:
            raise ValueError('Ome manifest has unknown columns: ' + ', '.join(sorted(unknown)))
        for line, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise ValueError(f'Manifest line {line} has {len(values)} columns; expected {len(header)}')
            row = dict(zip(header, values))
            name = row['ome_name'].strip()
            rows.append([
                name,
                _manifest_path(row['phenotype_bed'], 'phenotype_bed', line),
                _manifest_path(row['principal_components_tsv'], 'principal_components_tsv', line),
                _manifest_path(row.get('additional_covariates_tsv', ''), 'additional_covariates_tsv', line, optional=True),
            ])
    validate_multiomics_inputs([row[0] for row in rows], thresholds)
    # Metadata only: referenced matrices are localized later as explicit WDL Files.
    with output.open('w', encoding='utf-8', newline='') as handle:
        # WDL read_tsv splits literal tabs; it does not decode CSV quoting.
        for row in rows:
            handle.write('\t'.join(row) + '\n')
    LOGGER.info('Validated ome manifest: ome_count=%d', len(rows))
