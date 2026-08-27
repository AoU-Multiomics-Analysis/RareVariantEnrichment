version 1.0

task BuildCarrierDefinitions {
    input {
        File variant_carrier_audit_tsv_gz
        File variant_carriers_qc_json
        File carrier_definitions_json
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail

        echo "Starting carrier-definition construction"
        audit_count=$(gzip -cd "~{variant_carrier_audit_tsv_gz}" | tail -n +2 | wc -l | tr -d ' ')
        echo "Input audit row count: $audit_count"

        rare-variant-enrichment build-carrier-definitions \
            --audit "~{variant_carrier_audit_tsv_gz}" \
            --extraction-qc "~{variant_carriers_qc_json}" \
            --definitions "~{carrier_definitions_json}" \
            --container-image "~{docker_image}" \
            --output "carrier_definitions.tsv.gz" \
            --qc-output "carrier_definitions.qc.json"

        carrier_count=$(gzip -cd carrier_definitions.tsv.gz | tail -n +2 | wc -l | tr -d ' ')
        definition_count=$(python -c 'import json,sys; print(len(json.load(open(sys.argv[1], encoding="utf-8"))["definition_order"]))' carrier_definitions.qc.json)
        definition_names=$(python -c 'import json,sys; print(",".join(json.load(open(sys.argv[1], encoding="utf-8"))["definition_order"]))' carrier_definitions.qc.json)
        echo "Selected definition names: $definition_names"
        echo "Completed carrier-definition construction; definition count: $definition_count; carrier row count: $carrier_count"
    >>>

    output {
        File carrier_definitions_tsv_gz = "carrier_definitions.tsv.gz"
        File carrier_definitions_qc_json = "carrier_definitions.qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task PrepareProteinCodingGenes {
    input {
        File gene_annotation_gtf
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail

        echo "Starting protein-coding gene preparation; input byte count: $(wc -c < "~{gene_annotation_gtf}")"
        rare-variant-enrichment prepare-protein-coding-genes \
            --gtf "~{gene_annotation_gtf}" \
            --genes-output "protein_coding_genes.tsv" \
            --qc-output "protein_coding_genes.qc.json"
        gene_count=$(tail -n +2 protein_coding_genes.tsv | wc -l | tr -d ' ')
        echo "Completed protein-coding gene preparation; gene count: $gene_count"
    >>>

    output {
        File protein_coding_genes_tsv = "protein_coding_genes.tsv"
        File protein_coding_genes_qc_json = "protein_coding_genes.qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task PreparePcChunks {
    input {
        File principal_components_tsv
        Array[Int] pc_counts
        Int pc_counts_per_job
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File pc_counts_file = write_lines(pc_counts)

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        echo "Starting PC chunk preparation; requested PC count: ~{length(pc_counts)}"
        pc_count_values=()
        while IFS= read -r value; do
            pc_count_values+=("$value")
        done < "~{pc_counts_file}"
        pc_counts_csv="$(join_by_comma "${pc_count_values[@]}")"
        echo "PC chunk values: ${pc_counts_csv:-adaptive grid}"

        rare-variant-enrichment pc-chunks \
            --principal-components "~{principal_components_tsv}" \
            --pc-counts "$pc_counts_csv" \
            --pc-counts-per-job "~{pc_counts_per_job}" \
            --output "pc_chunks.json"

        chunk_count=$(python -c 'import json,sys; print(len(json.load(open(sys.argv[1], encoding="utf-8"))))' pc_chunks.json)
        echo "Completed PC chunk preparation; chunk count: $chunk_count"
    >>>

    output {
        File pc_chunks_json = "pc_chunks.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task CalculateCarrierPcEnrichment {
    input {
        File phenotype_bed
        File carrier_definitions_tsv_gz
        File carrier_definitions_qc_json
        File principal_components_tsv
        File? additional_covariates_tsv
        File protein_coding_genes
        Array[Float] negative_z_thresholds
        Array[Int] pc_counts
        String pc_grid_mode
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
        Int preemptible
    }

    File negative_z_thresholds_file = write_lines(negative_z_thresholds)
    File pc_counts_file = write_lines(pc_counts)
    String additional_covariates_argument = if defined(additional_covariates_tsv) then "--additional-covariates \"~{select_first([additional_covariates_tsv])}\"" else ""

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        negative_z_threshold_values=()
        while IFS= read -r value; do
            negative_z_threshold_values+=("$value")
        done < "~{negative_z_thresholds_file}"
        if ((${#negative_z_threshold_values[@]} == 0)); then
            echo "negative_z_thresholds must contain at least one negative threshold" >&2
            exit 2
        fi
        pc_count_values=()
        while IFS= read -r value; do
            pc_count_values+=("$value")
        done < "~{pc_counts_file}"
        negative_z_thresholds_csv="$(join_by_comma "${negative_z_threshold_values[@]}")"
        pc_counts_csv="$(join_by_comma "${pc_count_values[@]}")"
        definition_count=$(python -c 'import json,sys; print(len(json.load(open(sys.argv[1], encoding="utf-8"))["definition_order"]))' "~{carrier_definitions_qc_json}")
        definition_names=$(python -c 'import json,sys; print(",".join(json.load(open(sys.argv[1], encoding="utf-8"))["definition_order"]))' "~{carrier_definitions_qc_json}")
        echo "Starting carrier/PC enrichment; definition count: $definition_count; PC count: ${#pc_count_values[@]}; threshold count: ${#negative_z_threshold_values[@]}"
        echo "Selected definition names: $definition_names"
        echo "PC chunk values: $pc_counts_csv"

        rare-variant-enrichment carrier-pc-enrichment \
            --phenotype-bed "~{phenotype_bed}" \
            --carrier-table "~{carrier_definitions_tsv_gz}" \
            --carrier-manifest "~{carrier_definitions_qc_json}" \
            --principal-components "~{principal_components_tsv}" \
            ~{additional_covariates_argument} \
            --protein-coding-genes "~{protein_coding_genes}" \
            --negative-z-thresholds="$negative_z_thresholds_csv" \
            --pc-counts "$pc_counts_csv" \
            --pc-grid-mode "~{pc_grid_mode}" \
            --results-output "carrier_pc_enrichment.tsv" \
            --summary-output "carrier_pc_enrichment.summary.json" \
            --gene-pc-qc-output "carrier_pc_enrichment.gene_pc_qc.tsv.gz" \
            --analysis-qc-output "carrier_pc_enrichment.analysis_qc.json"

        result_count=$(tail -n +2 carrier_pc_enrichment.tsv | wc -l | tr -d ' ')
        echo "Completed carrier/PC enrichment; result row count: $result_count"
    >>>

    output {
        File results_tsv = "carrier_pc_enrichment.tsv"
        File summary_json = "carrier_pc_enrichment.summary.json"
        File gene_pc_qc_tsv_gz = "carrier_pc_enrichment.gene_pc_qc.tsv.gz"
        File analysis_qc_json = "carrier_pc_enrichment.analysis_qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
        preemptible: preemptible
    }
}

task MergeCarrierPcEnrichment {
    input {
        Array[File] results_inputs
        Array[File] summary_inputs
        Array[File] gene_pc_qc_inputs
        Array[File] analysis_qc_inputs
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail

        echo "Starting carrier/PC merge; shard count: ~{length(results_inputs)}"
        rare-variant-enrichment merge-carrier-pc-enrichment \
            --results-input "~{sep="\" --results-input \"" results_inputs}" \
            --summary-input "~{sep="\" --summary-input \"" summary_inputs}" \
            --gene-pc-qc-input "~{sep="\" --gene-pc-qc-input \"" gene_pc_qc_inputs}" \
            --analysis-qc-input "~{sep="\" --analysis-qc-input \"" analysis_qc_inputs}" \
            --results-output "carrier_pc_enrichment.tsv" \
            --summary-output "carrier_pc_enrichment.summary.json" \
            --gene-pc-qc-output "carrier_pc_enrichment.gene_pc_qc.tsv.gz" \
            --analysis-qc-output "carrier_pc_enrichment.analysis_qc.json"

        result_count=$(tail -n +2 carrier_pc_enrichment.tsv | wc -l | tr -d ' ')
        echo "Completed carrier/PC merge; result row count: $result_count"
    >>>

    output {
        File results_tsv = "carrier_pc_enrichment.tsv"
        File summary_json = "carrier_pc_enrichment.summary.json"
        File gene_pc_qc_tsv_gz = "carrier_pc_enrichment.gene_pc_qc.tsv.gz"
        File analysis_qc_json = "carrier_pc_enrichment.analysis_qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task AnalyzeCarrierPcEnrichment {
    input {
        File results_tsv
        File carrier_definitions_qc_json
        Array[String] carrier_definitions
        Array[Float] selection_z_thresholds
        Float plateau_fraction
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File carrier_definitions_file = write_lines(carrier_definitions)
    File selection_z_thresholds_file = write_lines(selection_z_thresholds)

    command <<<
        set -euo pipefail

        join_by_comma() {
            local IFS=,
            printf '%s' "$*"
        }

        selection_z_threshold_values=()
        while IFS= read -r value; do
            selection_z_threshold_values+=("$value")
        done < "~{selection_z_thresholds_file}"
        carrier_definition_values=()
        while IFS= read -r value; do
            carrier_definition_values+=("$value")
        done < "~{carrier_definitions_file}"
        selection_z_thresholds_csv="$(join_by_comma "${selection_z_threshold_values[@]}")"
        if ((${#carrier_definition_values[@]} == 0)); then
            carrier_definitions_csv=$(python -c 'import json,sys; print(",".join(json.load(open(sys.argv[1], encoding="utf-8"))["definition_order"]))' "~{carrier_definitions_qc_json}")
        else
            carrier_definitions_csv="$(join_by_comma "${carrier_definition_values[@]}")"
        fi
        definition_count=$(python -c 'import sys; print(len(sys.argv[1].split(",")))' "$carrier_definitions_csv")
        echo "Starting carrier/PC analysis; selected definition count: $definition_count; threshold count: ${#selection_z_threshold_values[@]}"
        echo "Selected definition names: $carrier_definitions_csv"

        rare-variant-enrichment analyze-carrier-pc-enrichment \
            --results-input "~{results_tsv}" \
            --selection-output "carrier_pc_selection.json" \
            --plot-output "carrier_pc_enrichment.svg" \
            --carrier-definitions "$carrier_definitions_csv" \
            --selection-z-thresholds="$selection_z_thresholds_csv" \
            --plateau-fraction "~{plateau_fraction}"

        Rscript "/opt/rare-variant-enrichment/pc_sweep_qc.R" \
            --results-input "~{results_tsv}" \
            --selection-z-thresholds "$selection_z_thresholds_csv" \
            --carrier-definitions "$carrier_definitions_csv" \
            --summary-output "pc_sweep_qc_summary.tsv" \
            --plot-output "pc_sweep_qc_percent_max.png"

        summary_count=$(tail -n +2 pc_sweep_qc_summary.tsv | wc -l | tr -d ' ')
        echo "Completed carrier/PC analysis; QC summary row count: $summary_count"
    >>>

    output {
        File selection_json = "carrier_pc_selection.json"
        File plot_svg = "carrier_pc_enrichment.svg"
        File pc_sweep_qc_summary_tsv = "pc_sweep_qc_summary.tsv"
        File pc_sweep_qc_plot_png = "pc_sweep_qc_percent_max.png"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

workflow CarrierEnrichment {
    input {
        File variant_carrier_audit_tsv_gz
        File variant_carriers_qc_json
        File carrier_definitions_json
        File phenotype_bed
        File principal_components_tsv
        File? additional_covariates_tsv
        File gene_annotation_gtf
        Array[Float] negative_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
        Array[Float] selection_z_thresholds = [-3.0, -4.0, -5.0, -6.0]
        Float plateau_fraction = 0.95
        Array[Int] pc_counts = []
        Int pc_counts_per_job = 10
        Array[String] pc_selection_carrier_definitions = []
        Int pc_preemptible = 2
        String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
        Int prepare_cpu = 2
        Int prepare_memory_gb = 32
        Int prepare_disk_gb = 500
        Int analysis_cpu = 8
        Int analysis_memory_gb = 128
        Int analysis_disk_gb = 1000
        Int max_retries = 1
    }

    Int calculated_materialization_disk_gb = ceil(((size(variant_carrier_audit_tsv_gz, "GiB") + size(variant_carriers_qc_json, "GiB") + size(carrier_definitions_json, "GiB")) * 2.0 + 20.0))
    Int dynamic_materialization_disk_gb = if calculated_materialization_disk_gb > prepare_disk_gb then calculated_materialization_disk_gb else prepare_disk_gb
    Int calculated_prepare_disk_gb = ceil((size(gene_annotation_gtf, "GiB") * 2.0 + 20.0))
    Int dynamic_prepare_disk_gb = if calculated_prepare_disk_gb > prepare_disk_gb then calculated_prepare_disk_gb else prepare_disk_gb
    Float additional_covariates_size_gb = if defined(additional_covariates_tsv) then size(select_first([additional_covariates_tsv]), "GiB") else 0.0
    Int calculated_pc_chunk_disk_gb = ceil((size(principal_components_tsv, "GiB") * 2.0 + 20.0))
    String pc_grid_mode = if length(pc_counts) == 0 then "adaptive" else "explicit"

    call BuildCarrierDefinitions {
        input:
            variant_carrier_audit_tsv_gz = variant_carrier_audit_tsv_gz,
            variant_carriers_qc_json = variant_carriers_qc_json,
            carrier_definitions_json = carrier_definitions_json,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = dynamic_materialization_disk_gb,
            max_retries = max_retries
    }

    call PrepareProteinCodingGenes {
        input:
            gene_annotation_gtf = gene_annotation_gtf,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = dynamic_prepare_disk_gb,
            max_retries = max_retries
    }

    call PreparePcChunks {
        input:
            principal_components_tsv = principal_components_tsv,
            pc_counts = pc_counts,
            pc_counts_per_job = pc_counts_per_job,
            docker_image = docker_image,
            cpu = 1,
            memory_gb = 4,
            disk_gb = calculated_pc_chunk_disk_gb,
            max_retries = max_retries
    }

    Array[Array[Int]] pc_count_chunks = read_json(PreparePcChunks.pc_chunks_json)
    Int calculated_analysis_disk_gb = ceil(((size(phenotype_bed, "GiB") + size(BuildCarrierDefinitions.carrier_definitions_tsv_gz, "GiB") + size(BuildCarrierDefinitions.carrier_definitions_qc_json, "GiB") + size(principal_components_tsv, "GiB") + size(gene_annotation_gtf, "GiB") + additional_covariates_size_gb) * 2.0 + 20.0))
    Int dynamic_analysis_disk_gb = if calculated_analysis_disk_gb > analysis_disk_gb then calculated_analysis_disk_gb else analysis_disk_gb

    scatter (pc_count_chunk in pc_count_chunks) {
        call CalculateCarrierPcEnrichment {
            input:
                phenotype_bed = phenotype_bed,
                carrier_definitions_tsv_gz = BuildCarrierDefinitions.carrier_definitions_tsv_gz,
                carrier_definitions_qc_json = BuildCarrierDefinitions.carrier_definitions_qc_json,
                principal_components_tsv = principal_components_tsv,
                additional_covariates_tsv = additional_covariates_tsv,
                protein_coding_genes = PrepareProteinCodingGenes.protein_coding_genes_tsv,
                negative_z_thresholds = negative_z_thresholds,
                pc_counts = pc_count_chunk,
                pc_grid_mode = pc_grid_mode,
                docker_image = docker_image,
                cpu = analysis_cpu,
                memory_gb = analysis_memory_gb,
                disk_gb = dynamic_analysis_disk_gb,
                max_retries = max_retries,
                preemptible = pc_preemptible
        }
    }

    Int calculated_merge_disk_gb = ceil(((size(CalculateCarrierPcEnrichment.results_tsv, "GiB") + size(CalculateCarrierPcEnrichment.summary_json, "GiB") + size(CalculateCarrierPcEnrichment.gene_pc_qc_tsv_gz, "GiB") + size(CalculateCarrierPcEnrichment.analysis_qc_json, "GiB")) * 2.0 + 20.0))
    Int dynamic_merge_disk_gb = if calculated_merge_disk_gb > analysis_disk_gb then calculated_merge_disk_gb else analysis_disk_gb

    call MergeCarrierPcEnrichment {
        input:
            results_inputs = CalculateCarrierPcEnrichment.results_tsv,
            summary_inputs = CalculateCarrierPcEnrichment.summary_json,
            gene_pc_qc_inputs = CalculateCarrierPcEnrichment.gene_pc_qc_tsv_gz,
            analysis_qc_inputs = CalculateCarrierPcEnrichment.analysis_qc_json,
            docker_image = docker_image,
            cpu = analysis_cpu,
            memory_gb = analysis_memory_gb,
            disk_gb = dynamic_merge_disk_gb,
            max_retries = max_retries
    }

    call AnalyzeCarrierPcEnrichment {
        input:
            results_tsv = MergeCarrierPcEnrichment.results_tsv,
            carrier_definitions_qc_json = BuildCarrierDefinitions.carrier_definitions_qc_json,
            carrier_definitions = pc_selection_carrier_definitions,
            selection_z_thresholds = selection_z_thresholds,
            plateau_fraction = plateau_fraction,
            docker_image = docker_image,
            cpu = 1,
            memory_gb = 4,
            disk_gb = dynamic_merge_disk_gb,
            max_retries = max_retries
    }

    output {
        File carrier_definitions_tsv_gz = BuildCarrierDefinitions.carrier_definitions_tsv_gz
        File carrier_definitions_qc_json = BuildCarrierDefinitions.carrier_definitions_qc_json
        File results_tsv = MergeCarrierPcEnrichment.results_tsv
        File summary_json = MergeCarrierPcEnrichment.summary_json
        File gene_pc_qc_tsv_gz = MergeCarrierPcEnrichment.gene_pc_qc_tsv_gz
        File analysis_qc_json = MergeCarrierPcEnrichment.analysis_qc_json
        File pc_selection_json = AnalyzeCarrierPcEnrichment.selection_json
        File enrichment_plot_svg = AnalyzeCarrierPcEnrichment.plot_svg
        File pc_sweep_qc_summary_tsv = AnalyzeCarrierPcEnrichment.pc_sweep_qc_summary_tsv
        File pc_sweep_qc_plot_png = AnalyzeCarrierPcEnrichment.pc_sweep_qc_plot_png
        File protein_coding_genes_tsv = PrepareProteinCodingGenes.protein_coding_genes_tsv
        File protein_coding_genes_qc_json = PrepareProteinCodingGenes.protein_coding_genes_qc_json
    }
}
