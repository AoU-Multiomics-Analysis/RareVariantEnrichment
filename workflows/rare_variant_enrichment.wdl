version 1.0

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

        rare-variant-enrichment prepare-protein-coding-genes \
            --gtf "~{gene_annotation_gtf}" \
            --genes-output "protein_coding_genes.tsv" \
            --qc-output "protein_coding_genes.qc.json"
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

task CalculateLofPcEnrichment {
    input {
        File phenotype_bed
        File lof_carrier_table
        File principal_components_tsv
        File protein_coding_genes
        Array[Float] negative_z_thresholds
        Array[Int] pc_counts
        String pc_grid_mode
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    File negative_z_thresholds_file = write_lines(negative_z_thresholds)
    File pc_counts_file = write_lines(pc_counts)

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

        rare-variant-enrichment lof-pc-enrichment \
            --phenotype-bed "~{phenotype_bed}" \
            --lof-carriers "~{lof_carrier_table}" \
            --principal-components "~{principal_components_tsv}" \
            --protein-coding-genes "~{protein_coding_genes}" \
            --negative-z-thresholds="$negative_z_thresholds_csv" \
            --pc-counts "$pc_counts_csv" \
            --pc-grid-mode "~{pc_grid_mode}" \
            --results-output "lof_pc_enrichment.tsv" \
            --summary-output "lof_pc_enrichment.summary.json" \
            --gene-pc-qc-output "lof_pc_enrichment.gene_pc_qc.tsv.gz" \
            --analysis-qc-output "lof_pc_enrichment.analysis_qc.json"
    >>>

    output {
        File results_tsv = "lof_pc_enrichment.tsv"
        File summary_json = "lof_pc_enrichment.summary.json"
        File gene_pc_qc_tsv_gz = "lof_pc_enrichment.gene_pc_qc.tsv.gz"
        File analysis_qc_json = "lof_pc_enrichment.analysis_qc.json"
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

        pc_count_values=()
        while IFS= read -r value; do
            pc_count_values+=("$value")
        done < "~{pc_counts_file}"
        pc_counts_csv="$(join_by_comma "${pc_count_values[@]}")"

        rare-variant-enrichment pc-chunks \
            --principal-components "~{principal_components_tsv}" \
            --pc-counts "$pc_counts_csv" \
            --pc-counts-per-job "~{pc_counts_per_job}" \
            --output "pc_chunks.json"
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

task MergeLofPcEnrichment {
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

        results_args=()
        while IFS= read -r path; do results_args+=(--results-input "$path"); done <<'EOF'
        ~{sep="\n" results_inputs}
        EOF
        summary_args=()
        while IFS= read -r path; do summary_args+=(--summary-input "$path"); done <<'EOF'
        ~{sep="\n" summary_inputs}
        EOF
        gene_pc_qc_args=()
        while IFS= read -r path; do gene_pc_qc_args+=(--gene-pc-qc-input "$path"); done <<'EOF'
        ~{sep="\n" gene_pc_qc_inputs}
        EOF
        analysis_qc_args=()
        while IFS= read -r path; do analysis_qc_args+=(--analysis-qc-input "$path"); done <<'EOF'
        ~{sep="\n" analysis_qc_inputs}
        EOF

        rare-variant-enrichment merge-lof-pc-enrichment \
            "${results_args[@]}" \
            "${summary_args[@]}" \
            "${gene_pc_qc_args[@]}" \
            "${analysis_qc_args[@]}" \
            --results-output "lof_pc_enrichment.tsv" \
            --summary-output "lof_pc_enrichment.summary.json" \
            --gene-pc-qc-output "lof_pc_enrichment.gene_pc_qc.tsv.gz" \
            --analysis-qc-output "lof_pc_enrichment.analysis_qc.json"
    >>>

    output {
        File results_tsv = "lof_pc_enrichment.tsv"
        File summary_json = "lof_pc_enrichment.summary.json"
        File gene_pc_qc_tsv_gz = "lof_pc_enrichment.gene_pc_qc.tsv.gz"
        File analysis_qc_json = "lof_pc_enrichment.analysis_qc.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

task AnalyzeLofPcEnrichment {
    input {
        File results_tsv
        Array[Float] selection_z_thresholds
        Float plateau_fraction
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

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
        selection_z_thresholds_csv="$(join_by_comma "${selection_z_threshold_values[@]}")"

        rare-variant-enrichment analyze-lof-pc-enrichment \
            --results-input "~{results_tsv}" \
            --selection-output "lof_pc_selection.json" \
            --plot-output "lof_pc_enrichment.svg" \
            --selection-z-thresholds "$selection_z_thresholds_csv" \
            --plateau-fraction "~{plateau_fraction}"

        Rscript "/opt/rare-variant-enrichment/pc_sweep_qc.R" \
            --results-input "~{results_tsv}" \
            --selection-z-thresholds "$selection_z_thresholds_csv" \
            --summary-output "pc_sweep_qc_summary.tsv" \
            --plot-output "pc_sweep_qc_percent_max.png"
    >>>

    output {
        File selection_json = "lof_pc_selection.json"
        File plot_svg = "lof_pc_enrichment.svg"
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

workflow RareVariantEnrichment {
    input {
        File phenotype_bed
        File lof_carrier_table
        File principal_components_tsv
        File gene_annotation_gtf
        Array[Float] negative_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
        Array[Float] selection_z_thresholds = [-3.0, -4.0, -5.0, -6.0]
        Float plateau_fraction = 0.95
        Array[Int] pc_counts = []
        Int pc_counts_per_job = 10
        String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
        Int prepare_cpu = 2
        Int prepare_memory_gb = 32
        Int prepare_disk_gb = 500
        Int analysis_cpu = 8
        Int analysis_memory_gb = 128
        Int analysis_disk_gb = 1000
        Int max_retries = 1
    }

    Int calculated_prepare_disk_gb = ceil((size(gene_annotation_gtf, "GiB") * 2.0 + 20.0))
    Int dynamic_prepare_disk_gb = if calculated_prepare_disk_gb > prepare_disk_gb then calculated_prepare_disk_gb else prepare_disk_gb
    Int calculated_analysis_disk_gb = ceil(((size(phenotype_bed, "GiB") + size(lof_carrier_table, "GiB") + size(principal_components_tsv, "GiB") + size(gene_annotation_gtf, "GiB")) * 2.0 + 20.0))
    Int dynamic_analysis_disk_gb = if calculated_analysis_disk_gb > analysis_disk_gb then calculated_analysis_disk_gb else analysis_disk_gb
    Int calculated_pc_chunk_disk_gb = ceil((size(principal_components_tsv, "GiB") * 2.0 + 20.0))
    String pc_grid_mode = if length(pc_counts) == 0 then "adaptive" else "explicit"

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

    scatter (pc_count_chunk in pc_count_chunks) {
        call CalculateLofPcEnrichment {
            input:
                phenotype_bed = phenotype_bed,
                lof_carrier_table = lof_carrier_table,
                principal_components_tsv = principal_components_tsv,
                protein_coding_genes = PrepareProteinCodingGenes.protein_coding_genes_tsv,
                negative_z_thresholds = negative_z_thresholds,
                pc_counts = pc_count_chunk,
                pc_grid_mode = pc_grid_mode,
                docker_image = docker_image,
                cpu = analysis_cpu,
                memory_gb = analysis_memory_gb,
                disk_gb = dynamic_analysis_disk_gb,
                max_retries = max_retries
        }
    }

    Int calculated_merge_disk_gb = ceil(((size(CalculateLofPcEnrichment.results_tsv, "GiB") + size(CalculateLofPcEnrichment.summary_json, "GiB") + size(CalculateLofPcEnrichment.gene_pc_qc_tsv_gz, "GiB") + size(CalculateLofPcEnrichment.analysis_qc_json, "GiB")) * 2.0 + 20.0))
    Int dynamic_merge_disk_gb = if calculated_merge_disk_gb > analysis_disk_gb then calculated_merge_disk_gb else analysis_disk_gb

    call MergeLofPcEnrichment {
        input:
            results_inputs = CalculateLofPcEnrichment.results_tsv,
            summary_inputs = CalculateLofPcEnrichment.summary_json,
            gene_pc_qc_inputs = CalculateLofPcEnrichment.gene_pc_qc_tsv_gz,
            analysis_qc_inputs = CalculateLofPcEnrichment.analysis_qc_json,
            docker_image = docker_image,
            cpu = analysis_cpu,
            memory_gb = analysis_memory_gb,
            disk_gb = dynamic_merge_disk_gb,
            max_retries = max_retries
    }

    call AnalyzeLofPcEnrichment {
        input:
            results_tsv = MergeLofPcEnrichment.results_tsv,
            selection_z_thresholds = selection_z_thresholds,
            plateau_fraction = plateau_fraction,
            docker_image = docker_image,
            cpu = 1,
            memory_gb = 4,
            disk_gb = dynamic_merge_disk_gb,
            max_retries = max_retries
    }

    output {
        File results_tsv = MergeLofPcEnrichment.results_tsv
        File summary_json = MergeLofPcEnrichment.summary_json
        File gene_pc_qc_tsv_gz = MergeLofPcEnrichment.gene_pc_qc_tsv_gz
        File analysis_qc_json = MergeLofPcEnrichment.analysis_qc_json
        File pc_selection_json = AnalyzeLofPcEnrichment.selection_json
        File enrichment_plot_svg = AnalyzeLofPcEnrichment.plot_svg
        File pc_sweep_qc_summary_tsv = AnalyzeLofPcEnrichment.pc_sweep_qc_summary_tsv
        File pc_sweep_qc_plot_png = AnalyzeLofPcEnrichment.pc_sweep_qc_plot_png
        File protein_coding_genes_tsv = PrepareProteinCodingGenes.protein_coding_genes_tsv
        File protein_coding_genes_qc_json = PrepareProteinCodingGenes.protein_coding_genes_qc_json
    }
}
