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

workflow RareVariantEnrichment {
    input {
        File phenotype_bed
        File lof_carrier_table
        File principal_components_tsv
        File gene_annotation_gtf
        Array[Float] negative_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
        Array[Int] pc_counts = []
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

    call PrepareProteinCodingGenes {
        input:
            gene_annotation_gtf = gene_annotation_gtf,
            docker_image = docker_image,
            cpu = prepare_cpu,
            memory_gb = prepare_memory_gb,
            disk_gb = dynamic_prepare_disk_gb,
            max_retries = max_retries
    }

    call CalculateLofPcEnrichment {
        input:
            phenotype_bed = phenotype_bed,
            lof_carrier_table = lof_carrier_table,
            principal_components_tsv = principal_components_tsv,
            protein_coding_genes = PrepareProteinCodingGenes.protein_coding_genes_tsv,
            negative_z_thresholds = negative_z_thresholds,
            pc_counts = pc_counts,
            docker_image = docker_image,
            cpu = analysis_cpu,
            memory_gb = analysis_memory_gb,
            disk_gb = dynamic_analysis_disk_gb,
            max_retries = max_retries
    }

    output {
        File results_tsv = CalculateLofPcEnrichment.results_tsv
        File summary_json = CalculateLofPcEnrichment.summary_json
        File gene_pc_qc_tsv_gz = CalculateLofPcEnrichment.gene_pc_qc_tsv_gz
        File analysis_qc_json = CalculateLofPcEnrichment.analysis_qc_json
        File protein_coding_genes_tsv = PrepareProteinCodingGenes.protein_coding_genes_tsv
        File protein_coding_genes_qc_json = PrepareProteinCodingGenes.protein_coding_genes_qc_json
    }
}
