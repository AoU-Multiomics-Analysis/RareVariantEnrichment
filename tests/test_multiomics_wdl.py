import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest
import WDL


WORKFLOW = Path('workflows/multiomics_outliers.wdl')


def walk_workflow(node):
    if isinstance(node, WDL.Expr.Apply):
        assert not node.function_name.startswith('write_'), str(node)
    if isinstance(node, WDL.Tree.Call):
        for expression in node.inputs.values():
            walk_workflow(expression)
    else:
        for child in node.children:
            walk_workflow(child)


def test_wrapper_calls_existing_workflow_and_preserves_file_types():
    document = WDL.load(str(WORKFLOW))
    assert document.effective_wdl_version == '1.0'
    workflow = document.workflow
    manifest = next(item for item in workflow.inputs if item.name == 'ome_manifest')
    assert str(manifest.type) == 'File'
    assert 'lof_carrier_table' not in {item.name for item in workflow.inputs}
    assert 'intersection_z_thresholds' not in {item.name for item in workflow.inputs}
    scatter = next(item for item in workflow.body if isinstance(item, WDL.Tree.Scatter))
    assert scatter.variable == 'manifest_row'
    call = next(item for item in scatter.body if isinstance(item, WDL.Tree.Call))
    assert isinstance(call.callee, WDL.Tree.Workflow)
    assert call.callee.name == 'RareVariantEnrichment'
    for field in ('phenotype_bed', 'principal_components_tsv', 'lof_carrier_table', 'additional_covariates_tsv'):
        assert str(call.inputs[field]) == field
    files = [item for item in scatter.body if isinstance(item, WDL.Tree.Decl) and isinstance(item.type, WDL.Type.File)]
    conditional = next(item for item in scatter.body if isinstance(item, WDL.Tree.Conditional))
    files.extend(conditional.body)
    assert [item.name for item in files] == ['phenotype_bed', 'principal_components_tsv', 'lof_carrier_table', 'additional_covariates_tsv']
    row = ['rna', 'gs://bucket/rna.bed', 'gs://bucket/rna.pcs', 'gs://bucket/rna.lof', 'gs://bucket/rna.cov']
    environment = WDL.Env.Bindings().bind('manifest_row', WDL.Value.from_json(WDL.Type.Array(WDL.Type.String()), row))
    for index, declaration in enumerate(files, start=1):
        value = declaration.expr.eval(environment, WDL.StdLib.Base('1.0')).coerce(declaration.type)
        assert isinstance(value, WDL.Value.File)
        assert value.value == row[index]
        localized = WDL.Value.rewrite_paths(value, lambda file: '/localized/' + file.value.rsplit('/', 1)[-1])
        assert isinstance(localized, WDL.Value.File)
        assert localized.value == '/localized/' + row[index].rsplit('/', 1)[-1]
    assert conditional.expr.eval(environment, WDL.StdLib.Base('1.0')).value
    absent = environment.bind('manifest_row', WDL.Value.from_json(WDL.Type.Array(WDL.Type.String()), row[:4] + ['.']))
    assert not conditional.expr.eval(absent, WDL.StdLib.Base('1.0')).value
    outputs = {item.name: str(item.type) for item in workflow.outputs}
    assert outputs['matrix_results'] == 'Array[OmicsResult]'
    assert outputs['z_score_matrices'] == 'Array[File]'
    assert outputs['haplo_matrices'] == 'Array[File?]'
    assert str(call.inputs['ome_name']) == 'dataset_name'
    result = next(item for item in scatter.body if isinstance(item, WDL.Tree.Decl) and item.name == 'result')
    assert str(result.type.members['lof_carrier_table']) == 'File'
    assert str(result.expr.members['lof_carrier_table']) == 'lof_carrier_table'
    assert str(result.type.members['selected_pc_haplo_calls_tsv_gz']) == 'File?'
    assert str(call.inputs['haplo_logcpm_drop']) == 'haplo_logcpm_drop'
    assert outputs['intersection_matrices'] == 'Array[File]'
    intersection = next(item for item in workflow.body if isinstance(item, WDL.Tree.Call) and item.name == 'IntersectMultiOmicsOutliers')
    assert str(intersection.inputs['expression_haplo_matrices']) == 'expression_haplo_matrices'
    haplo_files = next(item for item in workflow.body if isinstance(item, WDL.Tree.Decl) and item.name == 'expression_haplo_matrices')
    assert str(haplo_files.type) == 'Array[File]'
    assert str(haplo_files.expr) == 'select_all(RunMatrix.selected_pc_haplo_calls_tsv_gz)'

    for node in [*workflow.inputs, *workflow.body, *workflow.outputs]:
        walk_workflow(node)


@pytest.mark.parametrize('generated_files_in_cloud', [False, True])
@pytest.mark.parametrize('with_expression', [False, True])
@pytest.mark.parametrize('task_name', ['PrepareOmicsManifest', 'IntersectMultiOmicsOutliers'])
def test_new_tasks_localize_files_before_list_creation_and_run_safely(
    tmp_path, task_name, with_expression, generated_files_in_cloud,
):
    document = WDL.load(str(WORKFLOW))
    task = next(task for task in document.tasks if task.name == task_name)
    local = tmp_path / "local ' \" $(touch INJECTION) `touch INJECTION` $HOME"
    local.mkdir()
    mapping = {}
    for index in range(2):
        path = local / f'matrix{index}.tsv.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write('gene_id\tS1\tS2\nG1\t-3\t1\n')
        mapping[f'gs://bucket/matrix{index}.tsv.gz'] = str(path)
    manifest_path = local / 'omes.tsv'
    manifest_path.write_text('ome_name\tphenotype_bed\tprincipal_components_tsv\tlof_carrier_table\n'
        'rna\tgs://bucket/rna.bed\tgs://bucket/rna.pc\tgs://bucket/rna.lof\n'
        'protein\tgs://bucket/protein.bed\tgs://bucket/protein.pc\tgs://bucket/protein.lof\n')
    matrix_uris = list(mapping)
    mapping['gs://bucket/omes.tsv'] = str(manifest_path)
    haplo_path = local / 'expression_haplo.tsv'
    haplo_path.write_text('gene_id\tS2\tS1\nG1\t1\t0\n')
    mapping['gs://bucket/expression_haplo.tsv'] = str(haplo_path)
    values = {
        'expression_haplo_matrices': ['gs://bucket/expression_haplo.tsv'] if with_expression else [],
        'ome_manifest': 'gs://bucket/omes.tsv',
        'dataset_ids': ['expression' if with_expression else 'rna', 'protein'],
        'z_score_matrices': matrix_uris, 'docker_image': 'unused',
        'cpu': 1, 'memory_gb': 1, 'disk_gb': 1, 'max_retries': 0,
    }
    class StdLib(WDL.StdLib.Base):
        def _virtualize_filename(self, filename):
            if generated_files_in_cloud:
                uri = 'gs://generated-files/' + Path(filename).name
                mapping[uri] = filename
                return uri
            return filename
        def _devirtualize_filename(self, filename):
            return mapping.get(filename, filename)
    stdlib = StdLib('1.0', write_dir=str(tmp_path))
    environment = WDL.Env.Bindings()
    for declaration in task.inputs:
        environment = environment.bind(declaration.name, WDL.Value.from_json(declaration.type, values[declaration.name]))
    for declaration in task.postinputs:
        environment = environment.bind(declaration.name, declaration.expr.eval(environment, stdlib))
    environment = environment.map(lambda item: WDL.Env.Binding(
        item.name, WDL.Value.rewrite_paths(item.value, lambda file: mapping.get(file.value, file.value))
    ))
    if generated_files_in_cloud:
        # Cromwell maps File results after evaluating each command placeholder.
        # A write_lines result wrapped in sub becomes String, so this final
        # localization cannot replace its GCS URI. Model that boundary here.
        parts = []
        for part in task.command.parts:
            if isinstance(part, str):
                parts.append(part)
            else:
                assert not part.options  # These tasks use only scalar placeholders.
                value = part.expr.eval(environment, stdlib)
                localized = WDL.Value.rewrite_paths(
                    value, lambda file: mapping.get(file.value, file.value)
                )
                parts.append(localized.coerce(WDL.Type.String()).value)
        command = ''.join(parts)
    else:
        command = task.command.eval(environment, stdlib).value
    assert 'gs://' not in command
    command = 'rare-variant-enrichment() { ' + shlex.quote(sys.executable) + ' -m rare_variant_enrichment.cli "$@"; }\n' + command
    result = subprocess.run(['bash', '-c', command], cwd=tmp_path, env={**os.environ, 'PYTHONPATH': str(Path('src').resolve())}, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / 'INJECTION').exists()
    if task_name == 'PrepareOmicsManifest':
        assert (tmp_path / 'normalized_omics.tsv').read_text().splitlines() == [
            'rna\tgs://bucket/rna.bed\tgs://bucket/rna.pc\tgs://bucket/rna.lof\t.',
            'protein\tgs://bucket/protein.bed\tgs://bucket/protein.pc\tgs://bucket/protein.lof\t.',
        ]
    else:
        manifest = list(csv.DictReader((tmp_path / 'intersection_manifest.tsv').open(), delimiter='\t'))
        assert len(manifest) == (2 if with_expression else 1)
        for entry in manifest:
            with gzip.open(tmp_path / 'intersections' / entry['matrix_file'], 'rt') as handle:
                assert handle.read().splitlines() == ['gene_id\tS1\tS2', 'G1\t0\t0' if entry['expression_haplo_required'] == 'True' else 'G1\t1\t0']
        assert json.loads((tmp_path / 'intersection_summary.json').read_text())['matrix_count'] == len(manifest)
