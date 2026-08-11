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


def _parsed_outputs_and_runtime_keys() -> dict[str, object]:
    script = """
import json
import sys

import WDL

document = WDL.load(sys.argv[1])
print(json.dumps({
    "outputs": {declaration.name: str(declaration.type) for declaration in document.workflow.outputs},
    "runtime_keys": {task.name: sorted(task.runtime) for task in document.tasks},
    "task_inputs": {
        task.name: {declaration.name: str(declaration.type) for declaration in task.inputs}
        for task in document.tasks
    },
    "task_outputs": {
        task.name: {declaration.name: str(declaration.type) for declaration in task.outputs}
        for task in document.tasks
    },
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


def _parsed_workflow_wiring() -> dict[str, str | None]:
    script = """
import json
import sys

import WDL

document = WDL.load(sys.argv[1])
workflow = document.workflow

def identifier(expression):
    if not isinstance(expression, WDL.Expr.Get):
        return None
    if not isinstance(expression.expr, WDL.Expr.Ident):
        return None
    return expression.expr.name

calls = {
    item.name: item
    for item in workflow.body
    if isinstance(item, WDL.Tree.Call)
}
scatter = next(
    item for item in workflow.body if isinstance(item, WDL.Tree.Scatter)
)
scatters = [
    item for item in workflow.body if isinstance(item, WDL.Tree.Scatter)
]
classify = next(
    item for item in scatter.body if isinstance(item, WDL.Tree.Call)
)
calculate = calls["CalculateEnrichment"]
print(json.dumps({
    "scatter_count": len(scatters),
    "prepare_vcf_chromosomes": identifier(
        calls["PrepareVcfIndex"].inputs["chromosomes"]
    ),
    "prepare_vat_chromosomes": identifier(
        calls["PrepareVatIndex"].inputs["chromosomes"]
    ),
    "prepare_phenotypes_chromosomes": identifier(
        calls["PreparePhenotypes"].inputs["chromosomes"]
    ),
    "scatter_variable": scatter.variable,
    "scatter_collection": identifier(scatter.expr),
    "classify_chromosome": identifier(classify.inputs["chromosome"]),
    "classify_vat": identifier(classify.inputs["variant_annotation_table"]),
    "classify_vat_tbi": identifier(classify.inputs["variant_annotation_table_tbi"]),
    "classify_vat_schema": identifier(classify.inputs["vat_schema_json"]),
    "classify_consequence_classes": identifier(classify.inputs["consequence_classes"]),
    "classify_maximum_gvs_maf": identifier(classify.inputs["maximum_gvs_maf"]),
    "classify_annotation_chunk_size_bp": identifier(
        classify.inputs["annotation_chunk_size_bp"]
    ),
    "calculate_features_tsv": identifier(calculate.inputs["features_tsv"]),
    "calculate_selected_chromosomes": identifier(
        calculate.inputs["selected_chromosomes"]
    ),
    "calculate_consequence_classes": identifier(
        calculate.inputs["consequence_classes"]
    ),
    "calculate_vat_schema": identifier(calculate.inputs["vat_schema_json"]),
    "calculate_loftee_enabled": identifier(calculate.inputs["loftee_enabled"]),
    "calculate_vcf_index_provenance": identifier(
        calculate.inputs["vcf_index_provenance"]
    ),
    "calculate_vat_index_provenance": identifier(
        calculate.inputs["vat_index_provenance"]
    ),
    "calculate_maximum_gvs_maf": identifier(calculate.inputs["maximum_gvs_maf"]),
    "calculate_annotation_chunk_size_bp": identifier(
        calculate.inputs["annotation_chunk_size_bp"]
    ),
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


def _rendered_task_shell_boundaries(
    consequence_classes: list[str] | None = None,
) -> dict[str, dict[str, object]]:
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
dangerous_image = 'example.invalid/image"; touch IMAGE_INJECTION; echo "'
dangerous_consequence = 'stop_gained; touch CONSEQUENCE_INJECTION'
dangerous_vcf_provenance = 'supplied; touch VCF_PROVENANCE_INJECTION'
dangerous_vat_provenance = 'generated; touch VAT_PROVENANCE_INJECTION'
configured_consequences = (
    [dangerous_consequence] if len(sys.argv) == 2 else json.loads(sys.argv[2])
)

def array(item_type, values):
    return WDL.Value.Array(item_type, values)

class TestStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

values = {
    "rare_variant_vcf": WDL.Value.File("/tmp/rare variants.vcf.gz"),
    "rare_variant_vcf_tbi": WDL.Value.File("/tmp/rare variants.vcf.gz.tbi"),
    "variant_annotation_table": WDL.Value.File("/tmp/annotation table.tsv.bgz"),
    "variant_annotation_table_tbi": WDL.Value.File("/tmp/annotation table.tsv.bgz.tbi"),
    "vat_schema_json": WDL.Value.File("/tmp/vat_schema.json"),
    "chromosomes": array(WDL.Type.String(), [WDL.Value.String(dangerous_chromosome)]),
    "docker_image": WDL.Value.String(dangerous_image),
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
    "consequence_classes": array(
        WDL.Type.String(), [WDL.Value.String(value) for value in configured_consequences]
    ),
    "maximum_gvs_maf": WDL.Value.Float(0.01),
    "annotation_chunk_size_bp": WDL.Value.Int(500),
    "maximum_distance_bp": WDL.Value.Int(1000),
    "carrier_pairs": array(WDL.Type.File(), [WDL.Value.File(dangerous_carrier)]),
    "chromosome_qc": array(WDL.Type.File(), [WDL.Value.File(dangerous_qc)]),
    "carrier_minimum_distances_tsv": WDL.Value.File("/tmp/carriers.tsv"),
    "source_carrier_minimum_distances_tsv": WDL.Value.File("/tmp/carriers.tsv"),
    "phenotype_qc_json": WDL.Value.File("/tmp/phenotype_qc.json"),
    "chromosome_qc_tsv": WDL.Value.File("/tmp/chromosome_qc.tsv"),
    "selected_chromosomes": array(
        WDL.Type.String(), [WDL.Value.String(dangerous_chromosome)]
    ),
    "loftee_enabled": WDL.Value.Boolean(True),
    "vcf_index_provenance": WDL.Value.String(dangerous_vcf_provenance),
    "vat_index_provenance": WDL.Value.String(dangerous_vat_provenance),
    "workflow_version": WDL.Value.String("0.3.0"),
    "max_retries": WDL.Value.Int(2),
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
    command = [*_miniwdl_python(), "-c", script, str(WORKFLOW)]
    if consequence_classes is not None:
        command.append(json.dumps(consequence_classes))
    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_wdl_public_input_and_default_contract():
    expected_autosomes = [f"chr{chromosome}" for chromosome in range(1, 23)]
    template_result = subprocess.run(
        [MINIWDL, "input_template", str(WORKFLOW)],
        text=True,
        capture_output=True,
        check=False,
    )
    assert template_result.returncode == 0, template_result.stderr
    assert json.loads(template_result.stdout) == {
        "RareVariantEnrichment.phenotype_bed": "File",
        "RareVariantEnrichment.rare_variant_vcf": "File",
        "RareVariantEnrichment.variant_annotation_table": "File",
    }

    inputs = _parsed_workflow_inputs()
    assert inputs == {
        "phenotype_bed": {"type": "File", "default": None},
        "rare_variant_vcf": {"type": "File", "default": None},
        "rare_variant_vcf_tbi": {"type": "File?", "default": None},
        "variant_annotation_table": {"type": "File", "default": None},
        "variant_annotation_table_tbi": {"type": "File?", "default": None},
        "chromosomes": {
            "type": "Array[String]",
            "default": expected_autosomes,
        },
        "z_thresholds": {"type": "Array[Float]", "default": [2.0, 3.0, 4.0, 5.0]},
        "exact_allele_counts": {"type": "Array[Int]", "default": [1, 2, 3, 4, 5]},
        "cumulative_allele_count_maxima": {
            "type": "Array[Int]",
            "default": [1, 2, 3, 5, 10],
        },
        "consequence_classes": {
            "type": "Array[String]",
            "default": [
                "splice_acceptor_variant",
                "splice_donor_variant",
                "stop_gained",
                "frameshift_variant",
                "stop_lost",
                "start_lost",
                "inframe_insertion",
                "inframe_deletion",
                "missense_variant",
                "protein_altering_variant",
                "splice_region_variant",
                "synonymous_variant",
                "coding_sequence_variant",
            ],
        },
        "maximum_gvs_maf": {"type": "Float", "default": 0.01},
        "annotation_chunk_size_bp": {"type": "Int", "default": 10000000},
        "distance_thresholds_bp": {
            "type": "Array[Int]",
            "default": [10000, 100000],
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
        "max_retries": {"type": "Int", "default": 1},
        "publish_carrier_audit": {"type": "Boolean", "default": False},
    }
    assert inputs["exact_allele_counts"]["default"] == [1, 2, 3, 4, 5]
    assert inputs["cumulative_allele_count_maxima"]["default"] == [1, 2, 3, 5, 10]
    assert inputs["distance_thresholds_bp"]["default"] == [10_000, 100_000]


def test_example_uses_default_autosomes():
    inputs = json.loads(
        (Path("examples") / "rare_variant_enrichment.inputs.json").read_text()
    )
    assert "RareVariantEnrichment.chromosomes" not in inputs


def test_wdl_workflow_wires_chromosomes_and_prepared_features():
    wiring = _parsed_workflow_wiring()
    assert wiring == {
        "scatter_count": 1,
        "prepare_vcf_chromosomes": "chromosomes",
        "prepare_vat_chromosomes": "chromosomes",
        "prepare_phenotypes_chromosomes": "chromosomes",
        "scatter_variable": "chromosome",
        "scatter_collection": "chromosomes",
        "classify_chromosome": "chromosome",
        "classify_vat": "PrepareVatIndex.vat",
        "classify_vat_tbi": "PrepareVatIndex.vat_tbi",
        "classify_vat_schema": "PrepareVatIndex.vat_schema_json",
        "classify_consequence_classes": "consequence_classes",
        "classify_maximum_gvs_maf": "maximum_gvs_maf",
        "classify_annotation_chunk_size_bp": "annotation_chunk_size_bp",
        "calculate_features_tsv": "PreparePhenotypes.features_tsv",
        "calculate_selected_chromosomes": "chromosomes",
        "calculate_consequence_classes": "consequence_classes",
        "calculate_vat_schema": "PrepareVatIndex.vat_schema_json",
        "calculate_loftee_enabled": "PrepareVatIndex.loftee_enabled",
        "calculate_vcf_index_provenance": "PrepareVcfIndex.index_provenance",
        "calculate_vat_index_provenance": "PrepareVatIndex.index_provenance",
        "calculate_maximum_gvs_maf": "maximum_gvs_maf",
        "calculate_annotation_chunk_size_bp": "annotation_chunk_size_bp",
    }
    assert wiring["classify_chromosome"] == wiring["scatter_variable"]


def test_wdl_materializes_shell_values_before_rendering_commands():
    rendered = _rendered_task_shell_boundaries()
    dangerous_chromosome = 'chr1"; touch CHROMOSOME_INJECTION; echo "'
    dangerous_tail = 'absolute"; touch TAIL_INJECTION; echo "'
    dangerous_carrier = "/tmp/carrier input; touch CARRIER_INJECTION"
    dangerous_qc = "/tmp/qc input; touch QC_INJECTION"
    dangerous_image = 'example.invalid/image"; touch IMAGE_INJECTION; echo "'
    dangerous_consequence = "stop_gained; touch CONSEQUENCE_INJECTION"
    dangerous_vcf_provenance = "supplied; touch VCF_PROVENANCE_INJECTION"
    dangerous_vat_provenance = "generated; touch VAT_PROVENANCE_INJECTION"

    payloads_by_task = {
        "PrepareVcfIndex": [dangerous_chromosome],
        "PrepareVatIndex": [dangerous_chromosome],
        "PreparePhenotypes": [dangerous_chromosome, dangerous_tail],
        "ClassifyChromosome": [dangerous_chromosome, dangerous_consequence],
        "GatherCarrierPairs": [dangerous_carrier, dangerous_qc],
        "CalculateEnrichment": [
            dangerous_tail,
            dangerous_chromosome,
            dangerous_image,
            dangerous_consequence,
            dangerous_vcf_provenance,
            dangerous_vat_provenance,
        ],
    }
    expected_files_by_task = {
        "PrepareVcfIndex": [[dangerous_chromosome]],
        "PrepareVatIndex": [[dangerous_chromosome]],
        "PreparePhenotypes": [[dangerous_chromosome], ["2.000000"], [dangerous_tail]],
        "DetermineMaximumDistance": [["1000"]],
        "ClassifyChromosome": [
            [dangerous_chromosome],
            [],
            ["1"],
            [dangerous_consequence],
        ],
        "GatherCarrierPairs": [[dangerous_carrier], [dangerous_qc]],
        "CalculateEnrichment": [
            [],
            ["1"],
            ["2.000000"],
            ["1000"],
            [dangerous_tail],
            [dangerous_chromosome],
            [dangerous_image],
            [dangerous_vcf_provenance],
            [dangerous_vat_provenance],
            [dangerous_consequence],
            ["true"],
            ["0.010000"],
            ["500"],
            ["0.3.0"],
        ],
    }

    for task_name, payloads in payloads_by_task.items():
        command = rendered[task_name]["command"]
        for payload in payloads:
            assert payload not in command
    for task_name, expected_files in expected_files_by_task.items():
        generated_files = list(rendered[task_name]["generated_files"].values())
        for expected_file in expected_files:
            assert expected_file in generated_files

    prepare_vat_command = rendered["PrepareVatIndex"]["command"]
    assert 'ln -s "/tmp/annotation table.tsv.bgz" annotations.tsv.bgz' in prepare_vat_command
    assert (
        'ln -s "/tmp/annotation table.tsv.bgz.tbi" annotations.tsv.bgz.tbi'
        in prepare_vat_command
    )
    classify_command = rendered["ClassifyChromosome"]["command"]
    assert '--vat annotations.tsv.bgz' in classify_command
    assert 'ln -s "/tmp/annotation table.tsv.bgz.tbi" annotations.tsv.bgz.tbi' in classify_command
    calculate_command = rendered["CalculateEnrichment"]["command"]
    assert '--vat-schema "/tmp/vat_schema.json"' in calculate_command


def test_wdl_materializes_empty_consequence_arrays_for_both_consumers():
    rendered = _rendered_task_shell_boundaries([])

    for task_name in ("ClassifyChromosome", "CalculateEnrichment"):
        assert rendered[task_name]["generated_files"]["consequence_classes_file"] == []
        assert '--consequence-classes "$consequence_classes_csv"' in rendered[task_name][
            "command"
        ]


def test_wdl_exposes_qc_provenance_optional_audit_and_required_retry_runtime():
    contract = _parsed_outputs_and_runtime_keys()
    assert contract["outputs"] == {
        "carrier_minimum_distances_tsv": "File?",
        "chromosome_qc_tsv": "File",
        "chromosome_query_regions": "Array[File]",
        "enrichment_json": "File",
        "enrichment_tsv": "File",
        "generated_or_validated_vat_tbi": "File",
        "generated_or_validated_vcf_tbi": "File",
        "phenotype_qc_json": "File",
        "vat_index_provenance": "String",
        "vat_schema_json": "File",
        "vcf_index_provenance": "String",
    }
    assert contract["task_inputs"]["PrepareVatIndex"] == {
        "variant_annotation_table": "File",
        "variant_annotation_table_tbi": "File?",
        "chromosomes": "Array[String]",
        "docker_image": "String",
        "cpu": "Int",
        "memory_gb": "Int",
        "disk_gb": "Int",
        "max_retries": "Int",
    }
    assert contract["task_outputs"]["PrepareVatIndex"] == {
        "vat": "File",
        "vat_tbi": "File",
        "vat_schema_json": "File",
        "loftee_enabled": "Boolean",
        "index_provenance": "String",
    }
    assert contract["task_inputs"]["CalculateEnrichment"]["vat_schema_json"] == "File"
    assert all(
        "maxRetries" in runtime_keys
        for runtime_keys in contract["runtime_keys"].values()
    )
