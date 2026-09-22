"""Keep the optional covariate File intact until PC-sweep command rendering."""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import pytest
import WDL


@pytest.mark.parametrize('with_covariates', [False, True])
def test_pc_sweep_localizes_optional_covariates_before_command_rendering(tmp_path, with_covariates):
    document = WDL.load('workflows/rare_variant_enrichment.wdl')
    task = next(task for task in document.tasks if task.name == 'CalculateLofPcEnrichment')
    fixtures = Path('tests/fixtures').resolve()
    mapping = {}
    for name, source in {
        'phenotype_bed': 'lof_pc_phenotypes.bed',
        'lof_carrier_table': 'lof_carriers.tsv',
        'principal_components_tsv': 'principal_components.tsv',
        'protein_coding_genes': 'protein_coding_genes.tsv',
        'additional_covariates_tsv': 'genetic_pcs.tsv',
    }.items():
        # Check literal quotes and shell syntax in carrier and covariate paths.
        directory = tmp_path / ("covariates ' \" $(touch INJECTION) `touch INJECTION` $HOME" if name in ('additional_covariates_tsv', 'lof_carrier_table') else 'data')
        directory.mkdir(exist_ok=True)
        local = directory / source
        if name == 'protein_coding_genes':
            local.write_text('gene_id\nENSG1\nENSG2\n')
        else:
            shutil.copyfile(fixtures / source, local)
        mapping['gs://bucket/' + name] = str(local)
    values = {
        **{name: 'gs://bucket/' + name for name in ('phenotype_bed', 'lof_carrier_table', 'principal_components_tsv', 'protein_coding_genes')},
        'additional_covariates_tsv': 'gs://bucket/additional_covariates_tsv' if with_covariates else None,
        'negative_z_thresholds': [-0.8], 'pc_counts': [0, 1], 'pc_grid_mode': 'explicit',
        'docker_image': 'unused', 'cpu': 1, 'memory_gb': 1, 'disk_gb': 1, 'max_retries': 0, 'preemptible': 0,
    }
    class StdLib(WDL.StdLib.Base):
        def _virtualize_filename(self, filename):
            return filename
        def _devirtualize_filename(self, filename):
            return filename
    stdlib = StdLib('1.0', write_dir=str(tmp_path))
    environment = WDL.Env.Bindings()
    for declaration in task.inputs:
        environment = environment.bind(declaration.name, WDL.Value.from_json(declaration.type, values[declaration.name]))
    for declaration in task.postinputs:
        environment = environment.bind(declaration.name, declaration.expr.eval(environment, stdlib).coerce(declaration.type))
    environment = environment.map(lambda binding: WDL.Env.Binding(
        binding.name, WDL.Value.rewrite_paths(binding.value, lambda file: mapping.get(file.value, file.value))
    ))
    command = task.command.eval(environment, stdlib).value
    assert 'gs://' not in command
    assert ('--additional-covariates' in command) == with_covariates
    command = 'rare-variant-enrichment() { ' + shlex.quote(sys.executable) + ' -m rare_variant_enrichment.cli "$@"; }\n' + command
    result = subprocess.run(['bash', '-c', command], cwd=tmp_path, env={**os.environ, 'PYTHONPATH': str(Path('src').resolve())}, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / 'INJECTION').exists()
    qc = json.loads((tmp_path / 'lof_pc_enrichment.analysis_qc.json').read_text())
    assert qc['additional_covariates_supplied'] == with_covariates
    assert qc['additional_covariate_count'] == (2 if with_covariates else 0)
