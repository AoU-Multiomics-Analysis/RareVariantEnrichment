import csv
import subprocess
import sys

import pytest


HEADER = 'ome_name\tphenotype_bed\tprincipal_components_tsv\tlof_carrier_table\tadditional_covariates_tsv\n'


def prepare(tmp_path, text):
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(text)
    return subprocess.run([
        sys.executable, '-m', 'rare_variant_enrichment.cli', 'prepare-omics-manifest',
        '--manifest', str(manifest),
        '--output', str(tmp_path / 'normalized.tsv'),
    ], text=True, capture_output=True)


def test_manifest_keeps_paths_as_metadata_and_optional_fields(tmp_path):
    result = prepare(tmp_path, HEADER +
        'protein\tgs://bucket/protein.bed.gz\tgs://bucket/protein.pcs.tsv\tgs://bucket/protein.lof.tsv\tgs://bucket/age.tsv\n'
        'splicing\tgs://bucket/splicing.bed.gz\tgs://bucket/splicing.pcs.tsv\tgs://bucket/splicing.lof.tsv\t\n')
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'normalized.tsv').read_text().splitlines() == [
        'protein\tgs://bucket/protein.bed.gz\tgs://bucket/protein.pcs.tsv\tgs://bucket/protein.lof.tsv\tgs://bucket/age.tsv',
        'splicing\tgs://bucket/splicing.bed.gz\tgs://bucket/splicing.pcs.tsv\tgs://bucket/splicing.lof.tsv\t.',
    ]


def test_manifest_allows_optional_columns_to_be_omitted_and_headers_reordered(tmp_path):
    result = prepare(tmp_path, 'principal_components_tsv\tlof_carrier_table\tome_name\tphenotype_bed\n'
        'gs://bucket/a.pc\tgs://bucket/shared.lof\ta\tgs://bucket/a.bed\n'
        'gs://bucket/b.pc\tgs://bucket/shared.lof\tb\tgs://bucket/b.bed\n')
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'normalized.tsv').read_text().splitlines()[0] == 'a\tgs://bucket/a.bed\tgs://bucket/a.pc\tgs://bucket/shared.lof\t.'


@pytest.mark.parametrize('text,error', [
    ('ome_name\tphenotype_bed\n', 'principal_components_tsv'),
    ('ome_name\tphenotype_bed\tprincipal_components_tsv\n', 'lof_carrier_table'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\t\t.\n', 'lof_carrier_table'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\t.\t.\n', 'lof_carrier_table'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\trelative.tsv\t.\n', 'absolute'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\tgs://bucket\t.\n', 'object'),
    ('ome_name\tome_name\tphenotype_bed\tprincipal_components_tsv\n', 'duplicate'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\tgs://bucket/lof\t.\n' * 2, 'unique'),
    (HEADER + 'a\trelative.bed\tgs://bucket/pc\tgs://bucket/lof\t.\n', 'absolute'),
    (HEADER + 'a\tgs://bucket\tgs://bucket/pc\tgs://bucket/lof\t.\n', 'object'),
    (HEADER + 'a\tgs://bucket/a\t\tgs://bucket/lof\t.\n', 'principal_components_tsv'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\n', 'columns'),
    (HEADER + 'a\tgs://bucket/a\tgs://bucket/pc\tgs://bucket/lof\t.\n', 'two'),
])
def test_manifest_rejects_invalid_inputs(tmp_path, text, error):
    result = prepare(tmp_path, text)
    assert result.returncode != 0
    assert error.lower() in result.stderr.lower()
    assert not (tmp_path / 'normalized.tsv').exists()



def test_manifest_preserves_quotes_for_wdl_read_tsv(tmp_path):
    manifest = tmp_path / 'manifest.tsv'
    rows = [
        ['rna', '/data/a "quoted".bed', '/data/pc.tsv', '/data/carriers.tsv', '.'],
        ['protein', "/data/b 'quoted'.bed", '/data/pc.tsv', '/data/carriers.tsv', '.'],
    ]
    with manifest.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(HEADER.strip().split('\t'))
        writer.writerows(rows)
    result = prepare(tmp_path, manifest.read_text())
    assert result.returncode == 0, result.stderr
    assert [line.split('\t') for line in (tmp_path / 'normalized.tsv').read_text().splitlines()] == rows
