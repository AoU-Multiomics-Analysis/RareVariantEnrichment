import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from rare_variant_enrichment.zscore_matrix import export_zscore_matrix


def export(tmp_path, values, *, selected=0, drop=1.0, covariates=False):
    bed = tmp_path / 'logcpm.bed'
    pcs = tmp_path / 'pcs.tsv'
    selection = tmp_path / 'selection.json'
    bed.write_text('#chr\tstart\tend\tgene_id\tS1\tS2\tS3\tS4\tS5\tS6\n' + ''.join(
        'chr1\t0\t1\t' + gene + '\t' + '\t'.join(map(str, row)) + '\n' for gene, row in values
    ))
    pcs.write_text('ID\tPC1\tPC2\n' + ''.join(f'S{i}\t{pc}\t{(-1) ** i}\n' for i, pc in enumerate([-1, -1, 0, 0, 1, 1], 1)))
    selection.write_text(json.dumps({'selection': {'selected_pc_count': selected}}))
    cov = None
    if covariates:
        cov = tmp_path / 'cov.tsv'
        cov.write_text('sample_id\tbatch\nS1\t-2\nS2\t2\nS3\t-0.5\nS4\t0.5\nS5\t0\nS6\t0\n')
    export_zscore_matrix(bed, pcs, selection, tmp_path / 'z.tsv.gz', tmp_path / 'qc.tsv.gz', tmp_path / 'summary.json',
                        additional_covariates_path=cov, haplo_matrix_output=tmp_path / 'haplo.tsv.gz', haplo_logcpm_drop=drop)
    with gzip.open(tmp_path / 'haplo.tsv.gz', 'rt') as handle:
        return list(csv.reader(handle, delimiter='\t'))


def test_haplo_uses_strict_drop_in_logcpm_units_and_keeps_all_genes(tmp_path):
    rows = export(tmp_path, [('ENSG1', [8, 12, 9, 11, 10, 10]), ('ENSG2', [7] * 6)])
    assert rows == [
        ['gene_id', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
        ['ENSG1', '1', '0', '0', '0', '0', '0'],
        ['ENSG2', '0', '0', '0', '0', '0', '0'],
    ]
    summary = json.loads((tmp_path / 'summary.json').read_text())['haplo']
    assert summary['logcpm_drop'] == 1.0
    assert summary['positive_cell_count'] == 1
    assert summary['missing_cell_count'] == 0
    assert summary['additional_covariates_used'] is False


def test_haplo_uses_selected_phenotype_pcs_but_not_fixed_covariates(tmp_path):
    # Expression = 10 + 4*PC1 + batch. PC-only adjustment leaves batch deviations.
    values = [('ENSG1', [4, 8, 9.5, 10.5, 14, 14])]
    rows = export(tmp_path, values, selected=1, covariates=True)
    assert rows[1] == ['ENSG1', '1', '0', '0', '0', '0', '0']
    with gzip.open(tmp_path / 'z.tsv.gz', 'rt') as handle:
        assert list(csv.reader(handle, delimiter='\t'))[1][1:] == ['NA'] * 6


def test_haplo_missing_observations_and_unfit_genes_are_na(tmp_path):
    rows = export(tmp_path, [('ENSG1', ['NA', 6, 10, 10, 12, 12]), ('ENSG2', ['NA'] * 5 + [2])])
    assert rows[1] == ['ENSG1', 'NA', '1', '0', '0', '0', '0']
    assert rows[2] == ['ENSG2'] + ['NA'] * 6
    summary = json.loads((tmp_path / 'summary.json').read_text())['haplo']
    assert summary['missing_cell_count'] == 7
    assert summary['exclusion_counts']['insufficient_dof'] == 1


def test_haplo_drop_is_configurable(tmp_path):
    rows = export(tmp_path, [('ENSG1', [8, 12, 9, 11, 10, 10])], drop=0.5)
    assert rows[1] == ['ENSG1', '1', '0', '1', '0', '0', '0']


@pytest.mark.parametrize('drop', [-1, float('nan'), float('inf')])
def test_haplo_rejects_invalid_drop(tmp_path, drop):
    with pytest.raises(ValueError, match='finite.*non-negative'):
        export(tmp_path, [('ENSG1', [8, 12, 9, 11, 10, 10])], drop=drop)
    assert not (tmp_path / 'haplo.tsv.gz').exists()


def test_haplo_cli_emits_matrix(tmp_path):
    fixtures = Path('tests/fixtures')
    selection = tmp_path / 'selection.json'
    selection.write_text('{"selection": {"selected_pc_count": 0}}')
    result = subprocess.run([
        sys.executable, '-m', 'rare_variant_enrichment.cli', 'export-zscore-matrix',
        '--phenotype-bed', str(fixtures / 'lof_pc_phenotypes.bed'),
        '--principal-components', str(fixtures / 'principal_components.tsv'),
        '--selection-input', str(selection), '--matrix-output', str(tmp_path / 'z.tsv.gz'),
        '--gene-qc-output', str(tmp_path / 'qc.tsv.gz'), '--summary-output', str(tmp_path / 'summary.json'),
        '--haplo-matrix-output', str(tmp_path / 'haplo.tsv.gz'), '--haplo-logcpm-drop', '1',
    ], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    with gzip.open(tmp_path / 'haplo.tsv.gz', 'rt') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    assert rows[1] == ['ENSG1', '1', '1', '0', '0', '0', '0']


def test_haplo_matches_adjusted_logcpm_mean_reference_with_missing_values(tmp_path):
    from rare_variant_enrichment.haplo_matrix import export_haplo_matrix

    rng = np.random.default_rng(32)
    pcs = rng.normal(size=(24, 3))
    y = 20 + pcs @ [4, -2, 1] + rng.normal(size=24) * 2
    incomplete = y.copy()
    incomplete[[2, 7]] = np.nan
    genes = [('ENSG1', y), ('ENSG2', incomplete), ('ENSG3', y + 100)]
    output = tmp_path / 'haplo.tsv.gz'
    export_haplo_matrix(genes, pcs, 2, [str(i) for i in range(24)], output, 1)
    with gzip.open(output, 'rt') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))[1:]
    for row, (_, expression) in zip(rows, genes):
        usable = np.isfinite(expression)
        design = np.column_stack([np.ones(usable.sum()), pcs[usable, :2]])
        beta = np.linalg.lstsq(design, expression[usable], rcond=None)[0]
        adjusted = expression[usable] - pcs[usable, :2] @ beta[1:]
        expected = np.full(24, 'NA', dtype=object)
        expected[usable] = np.where(adjusted < adjusted.mean() - 1, '1', '0')
        assert row[1:] == expected.tolist()
    assert rows[0][1:] == rows[2][1:]


def test_haplo_marks_rank_deficient_design_as_missing(tmp_path):
    from rare_variant_enrichment.haplo_matrix import export_haplo_matrix

    summary = export_haplo_matrix([('ENSG1', np.arange(6.0))], np.ones((6, 1)), 1,
                                  [str(i) for i in range(6)], tmp_path / 'haplo.tsv.gz', 1)
    with gzip.open(tmp_path / 'haplo.tsv.gz', 'rt') as handle:
        assert list(csv.reader(handle, delimiter='\t'))[1] == ['ENSG1'] + ['NA'] * 6
    assert summary['exclusion_counts']['rank_deficiency'] == 1


def test_haplo_zero_drop_does_not_call_constant_gene(tmp_path):
    assert export(tmp_path, [('ENSG1', [7] * 6)], selected=1, drop=0)[1] == ['ENSG1'] + ['0'] * 6


def test_haplo_zero_drop_does_not_call_constant_decimal_gene(tmp_path):
    from rare_variant_enrichment.haplo_matrix import export_haplo_matrix

    output = tmp_path / 'haplo.tsv.gz'
    export_haplo_matrix([('ENSG1', np.full(6, 0.1))], np.random.default_rng(12).normal(size=(6, 2)),
                        2, [f'S{i}' for i in range(6)], output, 0)
    with gzip.open(output, 'rt') as handle:
        assert list(csv.reader(handle, delimiter='\t'))[1] == ['ENSG1'] + ['0'] * 6
