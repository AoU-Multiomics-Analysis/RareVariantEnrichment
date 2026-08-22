from pathlib import Path

import WDL


WORKFLOW = Path("workflows/extract_variant_carriers.wdl")


def test_carrier_wdl_has_standalone_gene_matched_interface():
    document = WDL.load(str(WORKFLOW))
    workflow = document.workflow
    assert workflow is not None

    inputs = {declaration.name: str(declaration.type) for declaration in workflow.inputs}
    assert inputs == {
        "filtered_vcf": "File",
        "filtered_vcf_tbi": "File",
        "transcript_annotations": "File",
        "transcript_annotations_tbi": "File?",
        "chromosomes": "Array[String]",
        "annotation_chunk_size_bp": "Int",
        "docker_image": "String",
        "prepare_cpu": "Int",
        "prepare_memory_gb": "Int",
        "prepare_disk_gb": "Int",
        "scatter_cpu": "Int",
        "scatter_memory_gb": "Int",
        "scatter_disk_gb": "Int",
        "gather_cpu": "Int",
        "gather_memory_gb": "Int",
        "gather_disk_gb": "Int",
        "max_retries": "Int",
        "scatter_preemptible": "Int",
    }
    assert sorted(task.name for task in document.tasks) == [
        "ExtractChromosomeCarriers", "GatherVariantCarriers", "PrepareCarrierInputs"
    ]
    assert {output.name: str(output.type) for output in workflow.outputs} == {
        "variant_carrier_audit_tsv_gz": "File",
        "variant_carriers_tsv_gz": "File",
        "variant_carriers_qc_json": "File",
        "chromosome_qc_jsons": "Array[File]",
        "transcript_schema_json": "File",
        "generated_or_validated_transcript_annotations_tbi": "File",
        "transcript_index_provenance": "String",
    }
    source = WORKFLOW.read_text()
    assert "scatter (chromosome in chromosomes)" in source
    for forbidden in ("phenotype", "gtf", "distance", "calculate-enrichment"):
        assert forbidden not in source.casefold()


def test_carrier_wdl_commands_have_operational_logging():
    document = WDL.load(str(WORKFLOW))
    for task in document.tasks:
        command = str(task.command)
        assert "Starting" in command
        assert "Completed" in command
        assert "count" in command.casefold()
