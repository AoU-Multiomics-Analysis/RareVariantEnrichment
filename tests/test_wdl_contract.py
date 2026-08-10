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
