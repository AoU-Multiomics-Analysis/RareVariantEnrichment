from pathlib import Path

import WDL


WORKFLOW = Path("workflows/carrier_enrichment.wdl")


def test_carrier_enrichment_wdl_has_exact_public_contract():
    document = WDL.load(str(WORKFLOW))
    workflow = document.workflow
    assert workflow is not None

    assert {declaration.name: str(declaration.type) for declaration in workflow.inputs} == {
        "variant_carrier_audit_tsv_gz": "File",
        "variant_carriers_qc_json": "File",
        "carrier_definitions_json": "File",
        "phenotype_bed": "File",
        "principal_components_tsv": "File",
        "additional_covariates_tsv": "File?",
        "gene_annotation_gtf": "File",
        "negative_z_thresholds": "Array[Float]",
        "selection_z_thresholds": "Array[Float]",
        "plateau_fraction": "Float",
        "pc_counts": "Array[Int]",
        "pc_counts_per_job": "Int",
        "pc_selection_carrier_definitions": "Array[String]",
        "pc_preemptible": "Int",
        "docker_image": "String",
        "prepare_cpu": "Int",
        "prepare_memory_gb": "Int",
        "prepare_disk_gb": "Int",
        "analysis_cpu": "Int",
        "analysis_memory_gb": "Int",
        "analysis_disk_gb": "Int",
        "max_retries": "Int",
    }
    assert sorted(task.name for task in document.tasks) == [
        "AnalyzeCarrierPcEnrichment",
        "BuildCarrierDefinitions",
        "CalculateCarrierPcEnrichment",
        "MergeCarrierPcEnrichment",
        "PreparePcChunks",
        "PrepareProteinCodingGenes",
    ]
    assert {output.name: str(output.type) for output in workflow.outputs} == {
        "carrier_definitions_tsv_gz": "File",
        "carrier_definitions_qc_json": "File",
        "results_tsv": "File",
        "summary_json": "File",
        "gene_pc_qc_tsv_gz": "File",
        "analysis_qc_json": "File",
        "pc_selection_json": "File",
        "enrichment_plot_svg": "File",
        "pc_sweep_qc_summary_tsv": "File",
        "pc_sweep_qc_plot_png": "File",
        "protein_coding_genes_tsv": "File",
        "protein_coding_genes_qc_json": "File",
    }


def test_carrier_enrichment_wdl_scatter_and_commands_are_generic_and_logged():
    document = WDL.load(str(WORKFLOW))
    source = WORKFLOW.read_text()

    assert "scatter (pc_count_chunk in pc_count_chunks)" in source
    assert "scatter (carrier" not in source
    for command in (
        "build-carrier-definitions",
        "carrier-pc-enrichment",
        "merge-carrier-pc-enrichment",
        "analyze-carrier-pc-enrichment",
    ):
        assert command in source
    assert '--container-image "~{docker_image}"' in source
    assert '--carrier-definitions "$carrier_definitions_csv"' in source
    assert 'Rscript "/opt/rare-variant-enrichment/pc_sweep_qc.R"' in source

    for task in document.tasks:
        command = str(task.command)
        assert "Starting" in command
        assert "Completed" in command
        assert "count" in command.casefold()

    calculate = next(
        task for task in document.tasks if task.name == "CalculateCarrierPcEnrichment"
    )
    analyze = next(
        task for task in document.tasks if task.name == "AnalyzeCarrierPcEnrichment"
    )
    assert "write_lines(negative_z_thresholds)" in source
    assert "write_lines(pc_counts)" in source
    assert "write_lines(carrier_definitions)" in source
    assert (
        'size(BuildCarrierDefinitions.carrier_definitions_tsv_gz, "GiB")'
        in source
    )
    assert (
        'size(BuildCarrierDefinitions.carrier_definitions_qc_json, "GiB")'
        in source
    )
    assert "preemptible" in calculate.runtime
    assert "preemptible" not in analyze.runtime


def test_carrier_enrichment_wdl_is_registered_with_dockstore():
    dockstore = Path(".dockstore.yml").read_text()
    assert "primaryDescriptorPath: /workflows/carrier_enrichment.wdl" in dockstore
