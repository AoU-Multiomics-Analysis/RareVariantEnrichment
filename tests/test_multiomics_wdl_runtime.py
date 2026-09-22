"""The GitHub Actions test job builds the image and requires this smoke test."""

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess

import pytest

from test_wdl_runtime import FIXTURES, TEST_IMAGE, _require_wdl_runtime


@pytest.mark.parametrize("ome_names", [("expression", "protein"), ("protein", "splicing")])
def test_multiomics_wrapper_runs_each_matrix_and_emits_intersections(tmp_path, ome_names):
    miniwdl = _require_wdl_runtime()
    dataset = {
        'phenotype_bed': str((FIXTURES / 'lof_pc_phenotypes.bed').resolve()),
        'principal_components_tsv': str((FIXTURES / 'principal_components.tsv').resolve()),
    }
    carrier_tables = {}
    for index, name in enumerate(ome_names):
        path = tmp_path / f'{name} carriers.tsv'
        content = (FIXTURES / 'lof_carriers.tsv').read_text()
        if index == 1:
            content = content.replace('v1\tHC', 'v1\tLC')
        path.write_text(content)
        carrier_tables[name] = path
    ome_manifest = tmp_path / 'omes.tsv'
    ome_manifest.write_text('ome_name\tphenotype_bed\tprincipal_components_tsv\tlof_carrier_table\n' + ''.join(
        name + '\t' + dataset['phenotype_bed'] + '\t' + dataset['principal_components_tsv'] + '\t' + str(carrier_tables[name]) + '\n'
        for name in ome_names
    ))
    inputs = {
        'ome_manifest': str(ome_manifest),
        'gene_annotation_gtf': str((FIXTURES / 'gene_annotation.gtf').resolve()),
        'negative_z_thresholds': [-0.8], 'selection_z_thresholds': [-0.8],
        'intersection_z_thresholds': [-0.8, -1.4], 'pc_counts': [0],
        'haplo_logcpm_drop': 2.0,
        'docker_image': TEST_IMAGE, 'pc_counts_per_job': 1,
        'prepare_cpu': 1, 'prepare_memory_gb': 1, 'prepare_disk_gb': 1,
        'analysis_cpu': 1, 'analysis_memory_gb': 1, 'analysis_disk_gb': 1,
        'intersection_cpu': 1, 'intersection_memory_gb': 1, 'intersection_disk_gb': 1,
        'max_retries': 0,
    }
    inputs_path = tmp_path / 'inputs.json'
    inputs_path.write_text(json.dumps({'MultiOmicsOutliers.' + key: value for key, value in inputs.items()}))
    output_path = tmp_path / 'outputs.json'
    result = subprocess.run([
        miniwdl, 'run', str(Path('workflows/multiomics_outliers.wdl').resolve()),
        '-i', str(inputs_path), '-d', str(tmp_path / 'run'), '-o', str(output_path), '--no-cache',
    ], text=True, capture_output=True, timeout=180, env={
        **os.environ,
        # Local fixture paths are introduced by the trusted manifest, not by
        # top-level File inputs. Permit miniwdl to read those local fixtures.
        'MINIWDL__FILE_IO__ALLOW_ANY_INPUT': 'true',
    })
    assert result.returncode == 0, result.stderr
    outputs = json.loads(output_path.read_text())['outputs']
    matrices = outputs['MultiOmicsOutliers.matrix_results']
    assert [entry['name'] for entry in matrices] == list(ome_names)
    for index, entry in enumerate(matrices):
        assert Path(entry['lof_carrier_table']).read_text() == carrier_tables[entry['name']].read_text()
        rows = list(csv.DictReader(Path(entry['results_tsv']).open(), delimiter='\t'))
        hc = next(row for row in rows if row['carrier_definition'] == 'HC')
        assert int(hc['carrier_observations']) == (2 if index == 0 else 1)
        selection = json.loads(Path(entry['pc_selection_json']).read_text())
        summary = json.loads(Path(entry['selected_pc_z_scores_summary_json']).read_text())
        assert selection['selection']['selected_pc_count'] == summary['selected_pc_count'] == 0
        assert summary['gene_count'] == 3
        assert Path(entry['phenotype_bed']).is_file()
        if entry['name'] != 'expression':
            assert entry['selected_pc_haplo_calls_tsv_gz'] is None
            assert 'haplo' not in summary
            continue
        with gzip.open(entry['selected_pc_haplo_calls_tsv_gz'], 'rt') as handle:
            haplo = list(csv.reader(handle, delimiter='\t'))
        assert haplo[1] == ['ENSG1', '1', '0', '0', '0', '0', '0']
        assert summary['haplo']['logcpm_drop'] == 2.0
    haplo_paths = outputs['MultiOmicsOutliers.haplo_matrices']
    assert len(haplo_paths) == len(matrices)
    # miniwdl places each output alias under its own output-field directory.
    for path, entry in zip(haplo_paths, matrices):
        if entry['name'] != 'expression':
            assert path is None
            continue
        with gzip.open(path, 'rt') as array_handle, gzip.open(entry['selected_pc_haplo_calls_tsv_gz'], 'rt') as result_handle:
            assert array_handle.read() == result_handle.read()
    manifest = list(csv.DictReader(Path(outputs['MultiOmicsOutliers.intersection_manifest_tsv']).open(), delimiter='\t'))
    assert [float(row['z_threshold']) for row in manifest] == [-0.8, -1.4]
    expected_counts = [3, 3] if 'expression' in ome_names else [6, 3]
    assert [int(row['outlier_count']) for row in manifest] == expected_counts
    assert all(row['expression_haplo_required'].lower() == str('expression' in ome_names).lower() for row in manifest)
    assert all(row['datasets'] == ','.join(ome_names) for row in manifest)
    files = {Path(path).name: path for path in outputs['MultiOmicsOutliers.intersection_matrices']}
    assert len(files) == 2
    for row in manifest:
        with gzip.open(files[row['matrix_file']], 'rt') as handle:
            values = list(csv.reader(handle, delimiter='\t'))
        assert values[0] == ['gene_id', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']
        assert values[1][0] == 'ENSG1'
