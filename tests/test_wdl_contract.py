import json
import gzip
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


WORKFLOW = Path("workflows/rare_variant_enrichment.wdl")
MINIWDL = shutil.which("miniwdl") or "/private/tmp/rve-rebuild-venv/bin/miniwdl"


def _miniwdl_python() -> list[str]:
    launcher = Path(MINIWDL).read_text().splitlines()[0].removeprefix("#!")
    command = shlex.split(launcher)
    return command[1:] if Path(command[0]).name == "env" else command


def _inspect_workflow() -> dict[str, object]:
    script = """
import json
import sys
import WDL

document = WDL.load(sys.argv[1])
workflow = document.workflow
stdlib = WDL.StdLib.Base(document.effective_wdl_version)
inputs = {}
for declaration in workflow.inputs:
    default = None
    if declaration.expr is not None:
        default = declaration.expr.eval(WDL.Env.Bindings(), stdlib).json
    inputs[declaration.name] = {"type": str(declaration.type), "default": default}
calls = [item for item in workflow.body if isinstance(item, WDL.Tree.Call)]
scatters = [item for item in workflow.body if isinstance(item, WDL.Tree.Scatter)]
print(json.dumps({
    "inputs": inputs,
    "outputs": {item.name: str(item.type) for item in workflow.outputs},
    "tasks": sorted(item.name for item in document.tasks),
    "calls": [item.name for item in calls],
    "call_inputs": {item.name: {key: str(value) for key, value in item.inputs.items()} for item in calls},
    "scatters": [{
        "variable": item.variable,
        "expr": str(item.expr),
        "calls": [{
            "name": child.name,
            "inputs": {key: str(value) for key, value in child.inputs.items()},
        } for child in item.body if isinstance(child, WDL.Tree.Call)],
    } for item in scatters],
    "runtime_keys": {item.name: sorted(item.runtime) for item in document.tasks},
    "task_inputs": {item.name: {decl.name: str(decl.type) for decl in item.inputs} for item in document.tasks},
    "workflow_body_kinds": [type(item).__name__ for item in workflow.body],
}, sort_keys=True))
"""
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def _workflow_declarations() -> dict[str, str]:
    script = """
import json
import sys
import WDL

document = WDL.load(sys.argv[1])
print(json.dumps({item.name: str(item.expr) for item in document.workflow.body if isinstance(item, WDL.Tree.Decl)}, sort_keys=True))
"""
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_wdl_has_only_the_four_file_public_inputs_and_required_defaults():
    template = subprocess.run(
        [MINIWDL, "input_template", str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert template.returncode == 0, template.stderr
    assert json.loads(template.stdout) == {
        "RareVariantEnrichment.phenotype_bed": "File",
        "RareVariantEnrichment.lof_carrier_table": "File",
        "RareVariantEnrichment.principal_components_tsv": "File",
        "RareVariantEnrichment.gene_annotation_gtf": "File",
    }

    contract = _inspect_workflow()
    assert contract["inputs"] == {
        "phenotype_bed": {"type": "File", "default": None},
        "lof_carrier_table": {"type": "File", "default": None},
        "principal_components_tsv": {"type": "File", "default": None},
        "additional_covariates_tsv": {"type": "File?", "default": None},
        "gene_annotation_gtf": {"type": "File", "default": None},
        "haplo_logcpm_drop": {"type": "Float", "default": 1.0},
        "ome_name": {"type": "String", "default": ""},
        "negative_z_thresholds": {
            "type": "Array[Float]",
            "default": [-2.0, -3.0, -4.0, -5.0, -6.0],
        },
        "selection_z_thresholds": {
            "type": "Array[Float]",
            "default": [-3.0, -4.0, -5.0, -6.0],
        },
        "plateau_fraction": {"type": "Float", "default": 0.95},
        "pc_counts": {"type": "Array[Int]", "default": []},
        "pc_counts_per_job": {"type": "Int", "default": 10},
        "pc_preemptible": {"type": "Int", "default": 2},
        "docker_image": {
            "type": "String",
            "default": "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main",
        },
        "prepare_cpu": {"type": "Int", "default": 2},
        "prepare_memory_gb": {"type": "Int", "default": 32},
        "prepare_disk_gb": {"type": "Int", "default": 500},
        "analysis_cpu": {"type": "Int", "default": 8},
        "analysis_memory_gb": {"type": "Int", "default": 128},
        "analysis_disk_gb": {"type": "Int", "default": 1000},
        "max_retries": {"type": "Int", "default": 1},
    }


def test_wdl_scatter_merge_and_selected_matrix_outputs():
    contract = _inspect_workflow()
    assert contract["tasks"] == [
        "AnalyzeLofPcEnrichment",
        "CalculateLofPcEnrichment",
        "ExportSelectedPcZScores",
        "MergeLofPcEnrichment",
        "PreparePcChunks",
        "PrepareProteinCodingGenes",
    ]
    assert contract["calls"] == [
        "PrepareProteinCodingGenes",
        "PreparePcChunks",
        "MergeLofPcEnrichment",
        "AnalyzeLofPcEnrichment",
        "ExportSelectedPcZScores",
    ]
    assert contract["scatters"] == [{
        "variable": "pc_count_chunk",
        "expr": "pc_count_chunks",
        "calls": [{
            "name": "CalculateLofPcEnrichment",
            "inputs": {
                "phenotype_bed": "phenotype_bed",
                "lof_carrier_table": "lof_carrier_table",
                "principal_components_tsv": "principal_components_tsv",
                "additional_covariates_tsv": "additional_covariates_tsv",
                "protein_coding_genes": "PrepareProteinCodingGenes.protein_coding_genes_tsv",
                "negative_z_thresholds": "negative_z_thresholds",
                "pc_counts": "pc_count_chunk",
                "pc_grid_mode": "pc_grid_mode",
                "docker_image": "docker_image",
                "cpu": "analysis_cpu",
                "memory_gb": "analysis_memory_gb",
                "disk_gb": "dynamic_analysis_disk_gb",
                "max_retries": "max_retries",
                "preemptible": "pc_preemptible",
            },
        }],
    }]
    assert contract["outputs"] == {
        "results_tsv": "File",
        "summary_json": "File",
        "gene_pc_qc_tsv_gz": "File",
        "analysis_qc_json": "File",
        "pc_selection_json": "File",
        "selected_pc_z_scores_tsv_gz": "File",
        "selected_pc_haplo_calls_tsv_gz": "File?",
        "selected_pc_z_scores_gene_qc_tsv_gz": "File",
        "selected_pc_z_scores_summary_json": "File",
        "enrichment_plot_svg": "File",
        "pc_sweep_qc_summary_tsv": "File",
        "pc_sweep_qc_plot_png": "File",
        "protein_coding_genes_tsv": "File",
        "protein_coding_genes_qc_json": "File",
    }
    serialized = WORKFLOW.read_text()
    for forbidden in ("VCF", "vcf", "VAT", "vat", "tabix", "conditional"):
        assert forbidden not in serialized


def test_wdl_wires_chunk_preparation_merge_and_dynamic_disk_floors():
    contract = _inspect_workflow()
    assert contract["call_inputs"] == {
        "PrepareProteinCodingGenes": {
            "gene_annotation_gtf": "gene_annotation_gtf",
            "docker_image": "docker_image",
            "cpu": "prepare_cpu",
            "memory_gb": "prepare_memory_gb",
            "disk_gb": "dynamic_prepare_disk_gb",
            "max_retries": "max_retries",
        },
        "PreparePcChunks": {
            "principal_components_tsv": "principal_components_tsv",
            "pc_counts": "pc_counts",
            "pc_counts_per_job": "pc_counts_per_job",
            "docker_image": "docker_image",
            "cpu": "1",
            "memory_gb": "4",
            "disk_gb": "calculated_pc_chunk_disk_gb",
            "max_retries": "max_retries",
        },
        "MergeLofPcEnrichment": {
            "results_inputs": "CalculateLofPcEnrichment.results_tsv",
            "summary_inputs": "CalculateLofPcEnrichment.summary_json",
            "gene_pc_qc_inputs": "CalculateLofPcEnrichment.gene_pc_qc_tsv_gz",
            "analysis_qc_inputs": "CalculateLofPcEnrichment.analysis_qc_json",
            "docker_image": "docker_image",
            "cpu": "analysis_cpu",
            "memory_gb": "analysis_memory_gb",
            "disk_gb": "dynamic_merge_disk_gb",
            "max_retries": "max_retries",
        },
        "ExportSelectedPcZScores": {
            "phenotype_bed": "phenotype_bed",
            "principal_components_tsv": "principal_components_tsv",
            "additional_covariates_tsv": "additional_covariates_tsv",
            "selection_json": "AnalyzeLofPcEnrichment.selection_json",
            "haplo_logcpm_drop": "haplo_logcpm_drop",
            "ome_name": "ome_name",
            "docker_image": "docker_image",
            "cpu": "analysis_cpu",
            "memory_gb": "analysis_memory_gb",
            "disk_gb": "dynamic_analysis_disk_gb",
            "max_retries": "max_retries",
            "preemptible": "pc_preemptible",
        },
        "AnalyzeLofPcEnrichment": {
            "results_tsv": "MergeLofPcEnrichment.results_tsv",
            "selection_z_thresholds": "selection_z_thresholds",
            "plateau_fraction": "plateau_fraction",
            "docker_image": "docker_image",
            "cpu": "1",
            "memory_gb": "4",
            "disk_gb": "dynamic_merge_disk_gb",
            "max_retries": "max_retries",
        },
    }
    declarations = _workflow_declarations()
    assert declarations["calculated_prepare_disk_gb"] == (
        'ceil((size(gene_annotation_gtf,"GiB") * 2.0 + 20.0))'
    )
    assert declarations["dynamic_prepare_disk_gb"] == (
        "if calculated_prepare_disk_gb > prepare_disk_gb then "
        "calculated_prepare_disk_gb else prepare_disk_gb"
    )
    assert declarations["calculated_pc_chunk_disk_gb"] == (
        'ceil((size(principal_components_tsv,"GiB") * 2.0 + 20.0))'
    )
    assert declarations["pc_grid_mode"] == (
        'if length(pc_counts) == 0 then "adaptive" else "explicit"'
    )
    assert declarations["calculated_analysis_disk_gb"] == (
        'ceil(((size(phenotype_bed,"GiB") + size(lof_carrier_table,"GiB") + '
        'size(principal_components_tsv,"GiB") + size(gene_annotation_gtf,"GiB") + '
        'additional_covariates_size_gb) * 2.0 + 20.0))'
    )
    assert declarations["additional_covariates_size_gb"] == (
        'if defined(additional_covariates_tsv) then '
        'size(select_first([additional_covariates_tsv]),"GiB") else 0.0'
    )
    assert declarations["dynamic_analysis_disk_gb"] == (
        "if calculated_analysis_disk_gb > analysis_disk_gb then "
        "calculated_analysis_disk_gb else analysis_disk_gb"
    )
    assert declarations["calculated_merge_disk_gb"] == (
        'ceil(((size(CalculateLofPcEnrichment.results_tsv,"GiB") + '
        'size(CalculateLofPcEnrichment.summary_json,"GiB") + '
        'size(CalculateLofPcEnrichment.gene_pc_qc_tsv_gz,"GiB") + '
        'size(CalculateLofPcEnrichment.analysis_qc_json,"GiB")) * 2.0 + 20.0))'
    )
    assert declarations["dynamic_merge_disk_gb"] == (
        "if calculated_merge_disk_gb > analysis_disk_gb then "
        "calculated_merge_disk_gb else analysis_disk_gb"
    )


def test_wdl_task_interfaces_and_retries_are_complete():
    contract = _inspect_workflow()
    assert contract["task_inputs"] == {
        "PrepareProteinCodingGenes": {
            "gene_annotation_gtf": "File",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
        },
        "CalculateLofPcEnrichment": {
            "phenotype_bed": "File",
            "lof_carrier_table": "File",
            "principal_components_tsv": "File",
            "additional_covariates_tsv": "File?",
            "protein_coding_genes": "File",
            "negative_z_thresholds": "Array[Float]",
            "pc_counts": "Array[Int]",
            "pc_grid_mode": "String",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
            "preemptible": "Int",
        },
        "PreparePcChunks": {
            "principal_components_tsv": "File",
            "pc_counts": "Array[Int]",
            "pc_counts_per_job": "Int",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
        },
        "MergeLofPcEnrichment": {
            "results_inputs": "Array[File]",
            "summary_inputs": "Array[File]",
            "gene_pc_qc_inputs": "Array[File]",
            "analysis_qc_inputs": "Array[File]",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
        },
        "ExportSelectedPcZScores": {
            "phenotype_bed": "File",
            "principal_components_tsv": "File",
            "additional_covariates_tsv": "File?",
            "selection_json": "File",
            "haplo_logcpm_drop": "Float",
            "ome_name": "String",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
            "preemptible": "Int",
        },
        "AnalyzeLofPcEnrichment": {
            "results_tsv": "File",
            "selection_z_thresholds": "Array[Float]",
            "plateau_fraction": "Float",
            "docker_image": "String",
            "cpu": "Int",
            "memory_gb": "Int",
            "disk_gb": "Int",
            "max_retries": "Int",
        },
    }
    assert all("maxRetries" in keys for keys in contract["runtime_keys"].values())
    assert "preemptible" in contract["runtime_keys"]["CalculateLofPcEnrichment"]


def test_wdl_materializes_array_values_and_quotes_files_before_shell_use():
    dangerous_image = 'example.invalid/image"; touch IMAGE_INJECTION; echo "'
    script = r'''
import json
from pathlib import Path
import sys
import tempfile
import WDL

document = WDL.load(sys.argv[1])
dangerous_threshold = '-2; touch THRESHOLD_INJECTION'
dangerous_pc = '0; touch PC_INJECTION'
dangerous_image = 'example.invalid/image"; touch IMAGE_INJECTION; echo "'

class TestStdLib(WDL.StdLib.Base):
    def _devirtualize_filename(self, filename):
        return filename
    def _virtualize_filename(self, filename):
        return filename

def array(item_type, values):
    return WDL.Value.Array(item_type, values)

values = {
        'gene_annotation_gtf': WDL.Value.File('/tmp/gene annotation.gtf.gz'),
    'phenotype_bed': WDL.Value.File('/tmp/phenotype matrix.bed.gz'),
    'lof_carrier_table': WDL.Value.File('/tmp/lof carriers.tsv'),
    'principal_components_tsv': WDL.Value.File('/tmp/principal components.tsv'),
    'additional_covariates_tsv': WDL.Value.File('/tmp/genetic pcs.tsv'),
    'protein_coding_genes': WDL.Value.File('/tmp/protein coding.tsv'),
    'negative_z_thresholds': array(WDL.Type.Float(), [WDL.Value.Float(-2.0)]),
    'selection_z_thresholds': array(WDL.Type.Float(), [WDL.Value.Float(-3.0), WDL.Value.Float(-4.0)]),
    'plateau_fraction': WDL.Value.Float(0.95),
    'haplo_logcpm_drop': WDL.Value.Float(1.0),
    'ome_name': WDL.Value.String('expression'),
    'pc_counts': array(WDL.Type.Int(), []),
    'pc_grid_mode': WDL.Value.String('adaptive'),
    'pc_counts_per_job': WDL.Value.Int(1),
    'results_inputs': array(WDL.Type.File(), [WDL.Value.File('/tmp/results one.tsv'), WDL.Value.File('/tmp/results two.tsv')]),
    'summary_inputs': array(WDL.Type.File(), [WDL.Value.File('/tmp/summary one.json'), WDL.Value.File('/tmp/summary two.json')]),
    'gene_pc_qc_inputs': array(WDL.Type.File(), [WDL.Value.File('/tmp/gene qc one.tsv.gz'), WDL.Value.File('/tmp/gene qc two.tsv.gz')]),
    'analysis_qc_inputs': array(WDL.Type.File(), [WDL.Value.File('/tmp/analysis qc one.json'), WDL.Value.File('/tmp/analysis qc two.json')]),
    'results_tsv': WDL.Value.File('/tmp/results one.tsv'),
    'selection_json': WDL.Value.File('/tmp/selection.json'),
    'docker_image': WDL.Value.String(dangerous_image),
    'cpu': WDL.Value.Int(1), 'memory_gb': WDL.Value.Int(1), 'disk_gb': WDL.Value.Int(1),
        'max_retries': WDL.Value.Int(1),
        'preemptible': WDL.Value.Int(2),
}
rendered = {}
with tempfile.TemporaryDirectory() as temporary_directory:
    selection_file = Path(temporary_directory) / 'selection.json'
    selection_file.write_text('{"selection": {"selected_pc_count": 20}}')
    values['selection_json'] = WDL.Value.File(str(selection_file))
    for task in document.tasks:
        environment = WDL.Env.Bindings()
        for declaration in task.inputs:
            environment = environment.bind(declaration.name, values[declaration.name])
        task_directory = str(Path(temporary_directory) / task.name)
        Path(task_directory).mkdir()
        stdlib = TestStdLib(document.effective_wdl_version, write_dir=task_directory)
        generated_files = {}
        for declaration in task.postinputs:
            value = declaration.expr.eval(environment, stdlib).coerce(declaration.type)
            environment = environment.bind(declaration.name, value)
            if isinstance(value, WDL.Value.File):
                generated_files[declaration.name] = Path(value.value).read_text().splitlines()
        rendered[task.name] = {'command': task.command.eval(environment, stdlib).value, 'generated_files': generated_files}
print(json.dumps(rendered, sort_keys=True))
'''
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    rendered = json.loads(result.stdout)
    matrix = rendered["ExportSelectedPcZScores"]["command"]
    for suffix in ("z_scores.tsv.gz", "haplo_calls.tsv.gz", "z_scores.gene_qc.tsv.gz", "z_scores.summary.json"):
        assert "expression.20PCs." + suffix in matrix
    assert "selected_pc_" not in matrix
    analysis = rendered["CalculateLofPcEnrichment"]
    assert '"/tmp/phenotype matrix.bed.gz"' in analysis["command"]
    assert "/tmp/lof carriers.tsv" in shlex.split(analysis["command"])
    assert '"/tmp/principal components.tsv"' in analysis["command"]
    assert "'/tmp/genetic pcs.tsv'" in analysis["command"]
    assert '"/tmp/protein coding.tsv"' in analysis["command"]
    assert '--negative-z-thresholds="$negative_z_thresholds_csv"' in analysis["command"]
    assert '--pc-grid-mode "adaptive"' in analysis["command"]
    assert dangerous_image not in analysis["command"]
    assert analysis["generated_files"]["negative_z_thresholds_file"] == ["-2.000000"]
    assert analysis["generated_files"]["pc_counts_file"] == []
    chunks = rendered["PreparePcChunks"]
    assert '"/tmp/principal components.tsv"' in chunks["command"]
    assert '--pc-counts "$pc_counts_csv"' in chunks["command"]
    assert chunks["generated_files"]["pc_counts_file"] == []
    merge = rendered["MergeLofPcEnrichment"]
    assert '--results-input-list' not in merge["command"]
    assert '--summary-input-list' not in merge["command"]
    assert '--gene-pc-qc-input-list' not in merge["command"]
    assert '--analysis-qc-input-list' not in merge["command"]
    assert '--results-input "/tmp/results one.tsv"' in merge["command"]
    assert '--summary-input "/tmp/summary one.json"' in merge["command"]
    assert '--gene-pc-qc-input "/tmp/gene qc one.tsv.gz"' in merge["command"]
    assert '--analysis-qc-input "/tmp/analysis qc one.json"' in merge["command"]
    assert '--results-input "/tmp/results one.tsv" --results-input "/tmp/results two.tsv"' in merge["command"]
    assert '--summary-input "/tmp/summary one.json" --summary-input "/tmp/summary two.json"' in merge["command"]
    assert '--gene-pc-qc-input "/tmp/gene qc one.tsv.gz" --gene-pc-qc-input "/tmp/gene qc two.tsv.gz"' in merge["command"]
    assert '--analysis-qc-input "/tmp/analysis qc one.json" --analysis-qc-input "/tmp/analysis qc two.json"' in merge["command"]
    assert 'printf ' not in merge["command"]
    assert '<<' not in merge["command"]
    assert 'done <' not in merge["command"]
    shell_check = subprocess.run(
        ["bash", "-n"], input=merge["command"], text=True, capture_output=True
    )
    assert shell_check.returncode == 0, shell_check.stderr
    assert merge["generated_files"] == {}
    plot = rendered["AnalyzeLofPcEnrichment"]
    assert 'results-input "/tmp/results one.tsv"' in plot["command"]
    assert '--selection-z-thresholds="$selection_z_thresholds_csv"' in plot["command"]
    assert '--selection-z-thresholds "$selection_z_thresholds_csv"' in plot["command"]
    assert '--plateau-fraction "0.950000"' in plot["command"]
    assert 'Rscript "/opt/rare-variant-enrichment/pc_sweep_qc.R"' in plot["command"]
    assert '--results-input "' in plot["command"]
    assert '--summary-output "pc_sweep_qc_summary.tsv"' in plot["command"]
    assert '--plot-output "pc_sweep_qc_percent_max.png"' in plot["command"]
    assert dangerous_image not in plot["command"]
    assert plot["generated_files"]["selection_z_thresholds_file"] == ["-3.000000", "-4.000000"]


def test_workflow_scope_does_not_write_files():
    script = '''
import sys
import WDL

document = WDL.load(sys.argv[1])
def check(node):
    if isinstance(node, WDL.Expr.Apply):
        assert not node.function_name.startswith("write_"), str(node)
    if isinstance(node, WDL.Tree.Call):
        for expression in node.inputs.values():
            check(expression)
    else:
        for child in node.children:
            check(child)

for node in [*document.workflow.inputs, *document.workflow.body, *document.workflow.outputs]:
    check(node)
'''
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("ome_name", ["expression", "protein", "splicing", "Expression", ""])
@pytest.mark.parametrize("with_covariates", [False, True])
def test_matrix_task_localizes_cloud_files_and_runs_safely(tmp_path, with_covariates, ome_name):
    fixtures = Path(__file__).parent / "fixtures"
    # Shell metacharacters must remain literal parts of each local path.
    localized = tmp_path / "localized ' \" $(touch INJECTION) `touch INJECTION` $HOME"
    localized.mkdir()
    mapping = {}
    for name, source in {
        "phenotype_bed": "lof_pc_phenotypes.bed",
        "principal_components_tsv": "principal_components.tsv",
        "additional_covariates_tsv": "genetic_pcs.tsv",
    }.items():
        path = localized / source
        shutil.copyfile(fixtures / source, path)
        mapping[f"gs://test-bucket/{name}"] = str(path)
    selection_path = localized / "selection.json"
    selection_path.write_text('{"selection": {"selected_pc_count": 0}}')
    mapping["gs://test-bucket/selection_json"] = str(selection_path)
    script = '''
import json
import sys
import WDL

document = WDL.load(sys.argv[1])
task = next(task for task in document.tasks if task.name == "ExportSelectedPcZScores")
mapping = json.loads(sys.argv[2])
environment = WDL.Env.Bindings()
for declaration in task.inputs:
    if declaration.name == "ome_name":
        value = WDL.Value.String(sys.argv[4])
    elif isinstance(declaration.type, WDL.Type.File):
        value = WDL.Value.File("gs://test-bucket/" + declaration.name)
        if declaration.name == "additional_covariates_tsv" and sys.argv[3] == "false":
            value = WDL.Value.Null()
    elif isinstance(declaration.type, WDL.Type.Float):
        value = WDL.Value.Float(1.0)
    elif isinstance(declaration.type, WDL.Type.Int):
        value = WDL.Value.Int(1)
    else:
        value = WDL.Value.String("test-image")
    environment = environment.bind(declaration.name, value)
class LocalizedStdLib(WDL.StdLib.Base):
    def _devirtualize_filename(self, filename):
        return mapping.get(filename, filename)
stdlib = LocalizedStdLib(document.effective_wdl_version)
for declaration in task.postinputs:
    environment = environment.bind(declaration.name, declaration.expr.eval(environment, stdlib))
# Simulate the runner's localization step after input evaluation.
environment = environment.map(
    lambda binding: WDL.Env.Binding(
        binding.name, WDL.Value.rewrite_paths(binding.value, lambda f: mapping[f.value])
    )
)
print(task.command.eval(environment, stdlib).value)
'''
    rendered = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW), json.dumps(mapping),
         "true" if with_covariates else "false", ome_name],
        text=True, capture_output=True, check=False,
    )
    assert rendered.returncode == 0, rendered.stderr
    assert "gs://" not in rendered.stdout
    if not with_covariates:
        assert "--additional-covariates" not in rendered.stdout
    command = (
        "rare-variant-enrichment() { " + shlex.quote(sys.executable)
        + ' -m rare_variant_enrichment.cli "$@"; }\n' + rendered.stdout
    )
    result = subprocess.run(
        ["bash", "-c", command], cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(Path("src").resolve())},
        text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "INJECTION").exists()
    prefix = (ome_name or "phenotype") + ".0PCs"
    with gzip.open(tmp_path / f"{prefix}.z_scores.tsv.gz", "rt") as handle:
        rows = handle.read().splitlines()
    assert len(rows) == 4  # All three genes, including the noncoding gene.
    assert rows[0].split("\t") == ["gene_id", "S1", "S2", "S3", "S4", "S5"] + (
        [] if with_covariates else ["S6"]
    )

    haplo_path = tmp_path / f"{prefix}.haplo_calls.tsv.gz"
    summary = json.loads((tmp_path / f"{prefix}.z_scores.summary.json").read_text())
    if ome_name != "expression":
        assert not haplo_path.exists()
        assert "haplo" not in summary
        assert "--haplo-matrix-output" not in rendered.stdout
        return
    assert "haplo" in summary
    assert summary["haplo"]["additional_covariates_used"] == with_covariates
    assert summary["haplo"]["additional_covariate_count"] == (2 if with_covariates else 0)
    with gzip.open(haplo_path, "rt") as handle:
        haplo_rows = handle.read().splitlines()
    assert haplo_rows[0] == rows[0]
    assert len(haplo_rows) == len(rows)
    expected_haplo = ["0"] * (5 if with_covariates else 6)
    assert haplo_rows[1].split("\t")[1:] == expected_haplo
