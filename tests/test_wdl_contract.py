import json
import shlex
import shutil
import subprocess
from pathlib import Path


WORKFLOW = Path("workflows/rare_variant_enrichment.wdl")
MINIWDL = shutil.which("miniwdl") or "/opt/homebrew/bin/miniwdl"


def _miniwdl_python() -> list[str]:
    launcher = Path(MINIWDL).read_text().splitlines()[0].removeprefix("#!")
    command = shlex.split(launcher)
    if Path(command[0]).name == "env":
        command = command[1:]
    return command


def _parsed_workflow_inputs() -> dict[str, dict[str, object]]:
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
print(json.dumps(inputs, sort_keys=True))
"""
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def _rendered_task_shell_boundaries() -> dict[str, dict[str, object]]:
    script = r'''
import json
from pathlib import Path
import sys
import tempfile

import WDL

document = WDL.load(sys.argv[1])
dangerous_chromosome = 'chr1"; touch CHROMOSOME_INJECTION; echo "'
dangerous_tail = 'absolute"; touch TAIL_INJECTION; echo "'
dangerous_carrier = '/tmp/carrier input; touch CARRIER_INJECTION'
dangerous_qc = '/tmp/qc input; touch QC_INJECTION'

def array(item_type, values):
    return WDL.Value.Array(item_type, values)

class TestStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

values = {
    "rare_variant_vcf": WDL.Value.File("/tmp/rare variants.vcf.gz"),
    "rare_variant_vcf_tbi": WDL.Value.File("/tmp/rare variants.vcf.gz.tbi"),
    "chromosomes": array(WDL.Type.String(), [WDL.Value.String(dangerous_chromosome)]),
    "docker_image": WDL.Value.String("example.invalid/image:latest"),
    "cpu": WDL.Value.Int(1),
    "memory_gb": WDL.Value.Int(1),
    "disk_gb": WDL.Value.Int(1),
    "phenotype_bed": WDL.Value.File("/tmp/phenotypes.bed"),
    "vcf_samples": WDL.Value.File("/tmp/vcf_samples.txt"),
    "z_thresholds": array(WDL.Type.Float(), [WDL.Value.Float(2.0)]),
    "outlier_tail": WDL.Value.String(dangerous_tail),
    "distance_thresholds_bp": array(WDL.Type.Int(), [WDL.Value.Int(1000)]),
    "features_tsv": WDL.Value.File("/tmp/features.tsv"),
    "shared_samples": WDL.Value.File("/tmp/shared_samples.txt"),
    "chromosome": WDL.Value.String(dangerous_chromosome),
    "exact_allele_counts": array(WDL.Type.Int(), []),
    "cumulative_allele_count_maxima": array(WDL.Type.Int(), [WDL.Value.Int(1)]),
    "maximum_distance_bp": WDL.Value.Int(1000),
    "carrier_pairs": array(WDL.Type.File(), [WDL.Value.File(dangerous_carrier)]),
    "chromosome_qc": array(WDL.Type.File(), [WDL.Value.File(dangerous_qc)]),
    "carrier_minimum_distances_tsv": WDL.Value.File("/tmp/carriers.tsv"),
}

rendered = {}
with tempfile.TemporaryDirectory() as temporary_directory:
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
        rendered[task.name] = {
            "command": task.command.eval(environment, stdlib).value,
            "generated_files": generated_files,
        }

print(json.dumps(rendered, sort_keys=True))
'''
    result = subprocess.run(
        [*_miniwdl_python(), "-c", script, str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_wdl_public_input_and_default_contract():
    template_result = subprocess.run(
        [MINIWDL, "input_template", str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert template_result.returncode == 0, template_result.stderr
    assert json.loads(template_result.stdout) == {
        "RareVariantEnrichment.chromosomes": ["String"],
        "RareVariantEnrichment.phenotype_bed": "File",
        "RareVariantEnrichment.rare_variant_vcf": "File",
    }

    assert _parsed_workflow_inputs() == {
        "phenotype_bed": {"type": "File", "default": None},
        "rare_variant_vcf": {"type": "File", "default": None},
        "rare_variant_vcf_tbi": {"type": "File?", "default": None},
        "chromosomes": {"type": "Array[String]", "default": None},
        "z_thresholds": {"type": "Array[Float]", "default": [2.0, 3.0, 4.0, 5.0]},
        "exact_allele_counts": {"type": "Array[Int]", "default": [1, 2, 3, 4, 5]},
        "cumulative_allele_count_maxima": {
            "type": "Array[Int]",
            "default": [1, 2, 3, 5, 10],
        },
        "distance_thresholds_bp": {
            "type": "Array[Int]",
            "default": [1000, 10000, 100000, 1000000],
        },
        "outlier_tail": {"type": "String", "default": "absolute"},
        "docker_image": {
            "type": "String",
            "default": "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main",
        },
        "prepare_cpu": {"type": "Int", "default": 2},
        "prepare_memory_gb": {"type": "Int", "default": 8},
        "prepare_disk_gb": {"type": "Int", "default": 50},
        "scatter_cpu": {"type": "Int", "default": 2},
        "scatter_memory_gb": {"type": "Int", "default": 8},
        "scatter_disk_gb": {"type": "Int", "default": 20},
        "gather_cpu": {"type": "Int", "default": 2},
        "gather_memory_gb": {"type": "Int", "default": 16},
        "gather_disk_gb": {"type": "Int", "default": 50},
    }


def test_wdl_materializes_shell_values_before_rendering_commands():
    rendered = _rendered_task_shell_boundaries()
    dangerous_chromosome = 'chr1"; touch CHROMOSOME_INJECTION; echo "'
    dangerous_tail = 'absolute"; touch TAIL_INJECTION; echo "'
    dangerous_carrier = "/tmp/carrier input; touch CARRIER_INJECTION"
    dangerous_qc = "/tmp/qc input; touch QC_INJECTION"

    payloads_by_task = {
        "PrepareVcfIndex": [dangerous_chromosome],
        "PreparePhenotypes": [dangerous_chromosome, dangerous_tail],
        "ClassifyChromosome": [dangerous_chromosome],
        "GatherCarrierPairs": [dangerous_carrier, dangerous_qc],
        "CalculateEnrichment": [dangerous_tail],
    }
    expected_files_by_task = {
        "PrepareVcfIndex": [[dangerous_chromosome]],
        "PreparePhenotypes": [[dangerous_chromosome], ["2.000000"], [dangerous_tail]],
        "DetermineMaximumDistance": [["1000"]],
        "ClassifyChromosome": [[dangerous_chromosome], [], ["1"]],
        "GatherCarrierPairs": [[dangerous_carrier], [dangerous_qc]],
        "CalculateEnrichment": [[], ["1"], ["2.000000"], ["1000"], [dangerous_tail]],
    }

    for task_name, payloads in payloads_by_task.items():
        command = rendered[task_name]["command"]
        for payload in payloads:
            assert payload not in command
    for task_name, expected_files in expected_files_by_task.items():
        generated_files = list(rendered[task_name]["generated_files"].values())
        for expected_file in expected_files:
            assert expected_file in generated_files
