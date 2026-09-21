import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import pytest


def run_intersections(tmp_path, matrices, thresholds='-2\n-3\n', names=None, haplo=None):
    paths = []
    for index, text in enumerate(matrices):
        path = tmp_path / f'matrix_{index}.tsv.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
        paths.append(str(path))
    (tmp_path / 'paths.txt').write_text('\n'.join(paths) + '\n')
    (tmp_path / 'names.txt').write_text('\n'.join(names or ['rna', 'protein', 'splicing'][:len(paths)]) + '\n')
    (tmp_path / 'thresholds.txt').write_text(thresholds)
    haplo_args = []
    if haplo is not None:
        (tmp_path / 'haplo.tsv').write_text(haplo)
        (tmp_path / 'haplo_paths.txt').write_text(str(tmp_path / 'haplo.tsv') + '\n')
        haplo_args = ['--expression-haplo-file-list', str(tmp_path / 'haplo_paths.txt')]
    return subprocess.run([
        sys.executable, '-m', 'rare_variant_enrichment.cli', 'multiomics-intersections',
        '--matrix-file-list', str(tmp_path / 'paths.txt'),
        '--dataset-id-list', str(tmp_path / 'names.txt'),
        '--threshold-list', str(tmp_path / 'thresholds.txt'),
        '--output-directory', str(tmp_path / 'intersections'),
        '--manifest-output', str(tmp_path / 'manifest.tsv'),
        '--summary-output', str(tmp_path / 'summary.json'),
        *haplo_args,
    ], text=True, capture_output=True)


def read_outputs(tmp_path):
    with (tmp_path / 'manifest.tsv').open() as handle:
        manifest = list(csv.DictReader(handle, delimiter='\t'))
    outputs = {}
    for entry in manifest:
        with gzip.open(tmp_path / 'intersections' / entry['matrix_file'], 'rt') as handle:
            outputs[(entry['datasets'], float(entry['z_threshold']))] = list(csv.reader(handle, delimiter='\t'))
    return manifest, outputs


def test_same_threshold_intersections_align_ids_and_test_only_members(tmp_path):
    result = run_intersections(tmp_path, [
        'gene_id\tS1\tS2\tS3\nG1\t-1\t-4\t-3\nG2\t-4\t-2\tNA\nRNA_ONLY\t-9\t-9\t-9\n',
        'gene_id\tS3\tS1\tS2\tS4\nG2\t-4\t-3\t-1\t-9\nG1\t-3\t-3\t-4\t-9\n',
        'gene_id\tS2\tS1\nG1\t-4\t-3\nG2\t-3\t-4\n',
    ])
    assert result.returncode == 0, result.stderr
    manifest, outputs = read_outputs(tmp_path)
    assert len(manifest) == 8
    assert outputs['protein,splicing', -3] == [
        ['gene_id', 'S1', 'S2'], ['G2', '1', '0'], ['G1', '1', '1'],
    ]
    assert outputs['rna,protein,splicing', -3] == [
        ['gene_id', 'S1', 'S2'], ['G1', '0', '1'], ['G2', '1', '0'],
    ]
    assert outputs['rna,protein', -3] == [
        ['gene_id', 'S1', 'S2', 'S3'], ['G1', '0', '1', '1'], ['G2', '1', '0', 'NA'],
    ]
    assert outputs['rna,splicing', -2] == [
        ['gene_id', 'S1', 'S2'], ['G1', '0', '1'], ['G2', '1', '1'],
    ]
    row = next(row for row in manifest if row['datasets'] == 'rna,protein' and float(row['z_threshold']) == -3)
    assert (row['gene_count'], row['sample_count'], row['outlier_count'], row['nonoutlier_count'], row['missing_count']) == ('2', '3', '3', '2', '1')
    assert json.loads((tmp_path / 'summary.json').read_text())['matrix_count'] == 8


def test_missing_is_unknown_even_when_another_dataset_is_not_an_outlier(tmp_path):
    result = run_intersections(tmp_path, [
        'gene_id\tS1\tS2\nG1\tNA\t4\n', 'gene_id\tS1\tS2\nG1\t4\tNA\n',
    ])
    assert result.returncode == 0, result.stderr
    assert read_outputs(tmp_path)[1]['rna,protein', -2] == [['gene_id', 'S1', 'S2'], ['G1', 'NA', 'NA']]


@pytest.mark.parametrize('second,expected', [
    ('gene_id\tS2\nG1\t-3\n', [['gene_id'], ['G1']]),
    ('gene_id\tS1\nG2\t-3\n', [['gene_id', 'S1']]),
])
def test_empty_overlap_is_reported(tmp_path, second, expected):
    result = run_intersections(tmp_path, ['gene_id\tS1\nG1\t-3\n', second])
    assert result.returncode == 0, result.stderr
    manifest, outputs = read_outputs(tmp_path)
    assert outputs['rna,protein', -2] == expected
    assert all(int(row['outlier_count']) == 0 for row in manifest)


@pytest.mark.parametrize('text,error', [
    ('gene_id\tS1\tS1\nG1\t-3\t-3\n', 'Duplicate sample'),
    ('gene_id\tS1\nG1\t-3\nG1\t-2\n', 'Duplicate gene'),
    ('gene_id\tS1\nG1\tinf\n', 'finite'),
    ('gene_id\tS1\nG1\tbad\n', 'numeric'),
    ('gene_id\tS1\nG1\n', 'columns'),
    ('gene_id\tS1\n\t-3\n', 'empty gene'),
    ('feature\tS1\nG1\t-3\n', 'gene_id'),
])
def test_invalid_matrix_is_rejected(tmp_path, text, error):
    result = run_intersections(tmp_path, [text, 'gene_id\tS1\nG1\t-3\n'])
    assert result.returncode != 0
    assert error.lower() in result.stderr.lower()
    assert not (tmp_path / 'manifest.tsv').exists()


@pytest.mark.parametrize('thresholds,names,error', [
    ('-2\n-2\n', ['rna', 'protein'], 'unique'),
    ('0\n', ['rna', 'protein'], 'negative'),
    ('nan\n', ['rna', 'protein'], 'finite'),
    ('', ['rna', 'protein'], 'threshold'),
    ('-2\n', ['rna', 'rna'], 'unique'),
    ('-2\n', ['rna', '../escape'], 'name'),
    ('-2\n', ['rna'], 'two'),
])
def test_invalid_settings_fail(tmp_path, thresholds, names, error):
    result = run_intersections(tmp_path, ['gene_id\tS1\nG1\t-3\n'] * 2, thresholds, names)
    assert result.returncode != 0
    assert error.lower() in result.stderr.lower()


@pytest.mark.parametrize('cloud', [False, True])
def test_rejects_unreadable_or_unlocalized_matrix(tmp_path, cloud):
    from rare_variant_enrichment.multiomics import build_outlier_intersections
    path = Path('gs://bucket/matrix.tsv.gz') if cloud else tmp_path / 'missing.tsv'
    with pytest.raises(ValueError, match='localization error'):
        build_outlier_intersections([path, path], ['rna', 'protein'], [-3], tmp_path / 'out', tmp_path / 'manifest', tmp_path / 'summary')
    assert not (tmp_path / 'out').exists()


def test_preflight_rejects_duplicate_names_before_analysis(tmp_path):
    (tmp_path / 'names').write_text('rna\nrna\n')
    (tmp_path / 'thresholds').write_text('-3\n')
    result = subprocess.run([
        sys.executable, '-m', 'rare_variant_enrichment.cli', 'validate-multiomics-inputs',
        '--dataset-id-list', str(tmp_path / 'names'), '--threshold-list', str(tmp_path / 'thresholds'),
        '--output', str(tmp_path / 'validated'),
    ], text=True, capture_output=True)
    assert result.returncode != 0
    assert 'unique' in result.stderr
    assert not (tmp_path / 'validated').exists()


def test_real_enrichment_selection_exports_feed_intersections(tmp_path):
    from rare_variant_enrichment.lof_pc import calculate_lof_pc_enrichment, prepare_protein_coding_genes
    from rare_variant_enrichment.pc_selection import analyze_lof_pc_enrichment
    from rare_variant_enrichment.zscore_matrix import export_zscore_matrix
    from rare_variant_enrichment.multiomics import build_outlier_intersections
    fixture = Path(__file__).parent / 'fixtures'
    genes = tmp_path / 'genes.tsv'
    prepare_protein_coding_genes(fixture / 'gene_annotation.gtf', genes, tmp_path / 'genes.qc.json')
    matrix_paths = []
    for index in range(2):
        output = tmp_path / str(index)
        output.mkdir()
        calculate_lof_pc_enrichment(
            fixture / 'lof_pc_phenotypes.bed', fixture / 'lof_carriers.tsv',
            fixture / 'principal_components.tsv', genes, [-0.8], [0],
            output / 'results.tsv', output / 'summary.json', output / 'gene_qc.gz', output / 'analysis_qc.json',
        )
        analyze_lof_pc_enrichment(output / 'results.tsv', output / 'selection.json', output / 'plot.svg', selection_z_thresholds=[-0.8])
        matrix = output / 'zscores.tsv.gz'
        export_zscore_matrix(fixture / 'lof_pc_phenotypes.bed', fixture / 'principal_components.tsv', output / 'selection.json', matrix, output / 'matrix_qc.gz', output / 'matrix_summary.json')
        matrix_paths.append(matrix)
    build_outlier_intersections(matrix_paths, ['rna', 'protein'], [-0.8, -1.4], tmp_path / 'intersections', tmp_path / 'manifest.tsv', tmp_path / 'summary.json')
    manifest, outputs = read_outputs(tmp_path)
    assert [float(row['z_threshold']) for row in manifest] == [-0.8, -1.4]
    assert [int(row['outlier_count']) for row in manifest] == [6, 3]
    assert outputs['rna,protein', -1.4] == [
        ['gene_id', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
        ['ENSG1', '1', '0', '0', '0', '0', '0'],
        ['ENSG2', '0', '0', '0', '0', '0', '1'],
        ['ENSG3', '1', '0', '0', '0', '0', '0'],
    ]


def test_intersections_stream_gene_order_without_sorting_full_vectors(tmp_path, monkeypatch):
    from rare_variant_enrichment import multiomics
    original_connect = multiomics.sqlite3.connect
    plans = []
    def connect(*args, **kwargs):
        connection = original_connect(*args, **kwargs)
        def trace(statement):
            if statement.startswith('SELECT matrix_'):
                plans.extend(row[3] for row in connection.execute('EXPLAIN QUERY PLAN ' + statement))
        connection.set_trace_callback(trace)
        return connection
    monkeypatch.setattr(multiomics.sqlite3, 'connect', connect)
    paths = []
    for index in range(3):
        path = tmp_path / f'{index}.tsv'
        path.write_text('gene_id\tS1\nG2\t-3\nG1\t-1\n')
        paths.append(path)
    multiomics.build_outlier_intersections(paths, ['rna', 'protein', 'splicing'], [-3], tmp_path / 'out', tmp_path / 'index.tsv', tmp_path / 'summary.json')
    assert plans
    assert not any('TEMP B-TREE' in step for step in plans), plans


def test_expression_intersections_require_haplo_and_z_but_other_pairs_ignore_it(tmp_path):
    result = run_intersections(tmp_path, [
        'gene_id\tS1\tS2\tS3\tS4\nG1\t-4\t-4\t-4\t-4\nG2\t-4\t-4\t-4\t-4\n',
        'gene_id\tS4\tS2\tS1\tS3\nG2\t-4\t-4\t-4\t-4\nG1\t-1\t-4\t-4\t-4\n',
        'gene_id\tS1\tS2\tS3\tS4\nG1\t-4\t-4\t-4\t-4\nG2\t-4\t-4\t-4\t-4\n',
    ], names=['protein', 'expression', 'splicing'], haplo=
        'gene_id\tS3\tS1\tS4\tS2\nG1\tNA\t1\t1\t0\nG2\t1\t0\t0\t1\n')
    assert result.returncode == 0, result.stderr
    manifest, outputs = read_outputs(tmp_path)
    for threshold in (-2, -3):
        assert outputs['protein,expression', threshold] == [
            ['gene_id', 'S1', 'S2', 'S3', 'S4'], ['G1', '1', '0', 'NA', '0'], ['G2', '0', '1', '1', '0']]
        assert outputs['protein,expression,splicing', threshold] == outputs['protein,expression', threshold]
        assert outputs['expression,splicing', threshold] == [
            ['gene_id', 'S4', 'S2', 'S1', 'S3'], ['G2', '0', '1', '0', '1'], ['G1', '0', '0', '1', 'NA']]
        assert outputs['protein,splicing', threshold] == [
            ['gene_id', 'S1', 'S2', 'S3', 'S4'], ['G1', '1', '1', '1', '1'], ['G2', '1', '1', '1', '1']]
    for row in manifest:
        assert row['expression_haplo_required'].lower() == str('expression' in row['datasets'].split(',')).lower()
    summary = json.loads((tmp_path / 'summary.json').read_text())
    assert 'haplo' in summary['outlier_rule']
    assert 'haplo' in summary['missing_rule']


def test_expression_requires_haplo_input(tmp_path):
    result = run_intersections(tmp_path, ['gene_id\tS1\nG1\t-3\n'] * 2, names=['expression', 'protein'])
    assert result.returncode != 0
    assert 'requires' in result.stderr and 'haplo' in result.stderr


@pytest.mark.parametrize('haplo,error', [
    ('gene_id\tS1\nG1\t2\n', '0, 1, or NA'),
    ('gene_id\tS1\nG1\t-1\n', '0, 1, or NA'),
    ('gene_id\tS2\nG1\t1\n', 'sample IDs'),
    ('gene_id\tS1\nG2\t1\n', 'gene IDs'),
])
def test_expression_haplo_must_match_expression_ids_and_be_binary(tmp_path, haplo, error):
    result = run_intersections(tmp_path, ['gene_id\tS1\nG1\t-3\n'] * 2, names=['expression', 'protein'], haplo=haplo)
    assert result.returncode != 0
    assert error in result.stderr


def test_haplo_input_without_expression_is_rejected(tmp_path):
    result = run_intersections(tmp_path, ['gene_id\tS1\nG1\t-3\n'] * 2, names=['protein', 'splicing'], haplo='gene_id\tS1\nG1\t1\n')
    assert result.returncode != 0
    assert 'Cannot supply an expression haplo matrix without an expression dataset' in result.stderr


@pytest.mark.parametrize('expression_z,haplo_value', [('4', 'NA'), ('NA', '0')])
def test_missing_required_expression_evidence_remains_na_even_when_other_condition_fails(tmp_path, expression_z, haplo_value):
    result = run_intersections(tmp_path, [
        f'gene_id\tS1\nG1\t{expression_z}\n', 'gene_id\tS1\nG1\t-3\n',
    ], names=['expression', 'protein'], haplo=f'gene_id\tS1\nG1\t{haplo_value}\n')
    assert result.returncode == 0, result.stderr
    assert read_outputs(tmp_path)[1]['expression,protein', -2] == [['gene_id', 'S1'], ['G1', 'NA']]


def test_expression_haplo_rejects_unlocalized_cloud_path(tmp_path):
    from rare_variant_enrichment.multiomics import build_outlier_intersections
    matrix = tmp_path / 'z.tsv'
    matrix.write_text('gene_id\tS1\nG1\t-3\n')
    with pytest.raises(ValueError, match='Input localization error: unresolved cloud URI'):
        build_outlier_intersections([matrix, matrix], ['expression', 'protein'], [-3],
                                    tmp_path / 'out', tmp_path / 'manifest', tmp_path / 'summary',
                                    expression_haplo_path=Path('gs://bucket/haplo.tsv.gz'))
    assert not (tmp_path / 'out').exists()
