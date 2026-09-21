version 1.0

import "rare_variant_enrichment.wdl" as enrichment

struct OmicsResult {
    String name
    File phenotype_bed
    File principal_components_tsv
    File? additional_covariates_tsv
    File selected_pc_z_scores_tsv_gz
    File? selected_pc_haplo_calls_tsv_gz
    File selected_pc_z_scores_gene_qc_tsv_gz
    File selected_pc_z_scores_summary_json
    File results_tsv
    File summary_json
    File gene_pc_qc_tsv_gz
    File analysis_qc_json
    File pc_selection_json
    File enrichment_plot_svg
    File pc_sweep_qc_summary_tsv
    File pc_sweep_qc_plot_png
    File protein_coding_genes_tsv
    File protein_coding_genes_qc_json
}

task PrepareOmicsManifest {
    input {
        File ome_manifest
        Array[Float] thresholds
        String docker_image
        Int max_retries
    }

    command <<<
        set -euo pipefail
        echo "Starting multi-omics input validation" >&2
        rare-variant-enrichment prepare-omics-manifest \
            --manifest '~{sub(ome_manifest, "'", "'\"'\"'")}' \
            --threshold-list '~{sub(write_lines(thresholds), "'", "'\"'\"'")}' \
            --output "normalized_omics.tsv"
        echo "Completed multi-omics input validation" >&2
    >>>

    output {
        File normalized_manifest_tsv = "normalized_omics.tsv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        maxRetries: max_retries
    }
}

task IntersectMultiOmicsOutliers {
    input {
        Array[File] z_score_matrices
        Array[File] expression_haplo_matrices
        Array[String] dataset_ids
        Array[Float] thresholds
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
    }

    command <<<
        set -euo pipefail
        echo "Starting multi-omics outlier intersections" >&2
        # Build the file list here, after the matrix files have been localized.
        rare-variant-enrichment multiomics-intersections \
            --matrix-file-list '~{sub(write_lines(z_score_matrices), "'", "'\"'\"'")}' \
            --expression-haplo-file-list '~{sub(write_lines(expression_haplo_matrices), "'", "'\"'\"'")}' \
            --dataset-id-list '~{sub(write_lines(dataset_ids), "'", "'\"'\"'")}' \
            --threshold-list '~{sub(write_lines(thresholds), "'", "'\"'\"'")}' \
            --output-directory "intersections" \
            --manifest-output "intersection_manifest.tsv" \
            --summary-output "intersection_summary.json"
        echo "Completed multi-omics outlier intersections" >&2
    >>>

    output {
        Array[File] matrices = glob("intersections/*.tsv.gz")
        File manifest_tsv = "intersection_manifest.tsv"
        File summary_json = "intersection_summary.json"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        maxRetries: max_retries
    }
}

workflow MultiOmicsOutliers {
    input {
        File ome_manifest
        File lof_carrier_table
        File gene_annotation_gtf
        Float haplo_logcpm_drop = 1.0
        Array[Float] intersection_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
        Array[Float] negative_z_thresholds = [-2.0, -3.0, -4.0, -5.0, -6.0]
        Array[Float] selection_z_thresholds = [-3.0, -4.0, -5.0, -6.0]
        Float plateau_fraction = 0.95
        Array[Int] pc_counts = []
        Int pc_counts_per_job = 10
        Int pc_preemptible = 2
        String docker_image = "ghcr.io/aou-multiomics-analysis/rarevariantenrichment:main"
        Int prepare_cpu = 2
        Int prepare_memory_gb = 32
        Int prepare_disk_gb = 500
        Int analysis_cpu = 8
        Int analysis_memory_gb = 128
        Int analysis_disk_gb = 1000
        Int intersection_cpu = 1
        Int intersection_memory_gb = 8
        Int intersection_disk_gb = 1000
        Int max_retries = 1
    }

    call PrepareOmicsManifest {
        input:
            ome_manifest = ome_manifest,
            thresholds = intersection_z_thresholds,
            docker_image = docker_image,
            max_retries = max_retries
    }

    # Manifest paths start as metadata. Declare Files before downstream localization.
    Array[Array[String]] manifest_rows = read_tsv(PrepareOmicsManifest.normalized_manifest_tsv)
    scatter (manifest_row in manifest_rows) {
        String dataset_name = manifest_row[0]
        File phenotype_bed = manifest_row[1]
        File principal_components_tsv = manifest_row[2]
        if (manifest_row[3] != ".") {
            File additional_covariates_tsv = manifest_row[3]
        }

        call enrichment.RareVariantEnrichment as RunMatrix {
            input:
                phenotype_bed = phenotype_bed,
                principal_components_tsv = principal_components_tsv,
                additional_covariates_tsv = additional_covariates_tsv,
                lof_carrier_table = lof_carrier_table,
                gene_annotation_gtf = gene_annotation_gtf,
                haplo_logcpm_drop = haplo_logcpm_drop,
                ome_name = dataset_name,
                negative_z_thresholds = negative_z_thresholds,
                selection_z_thresholds = selection_z_thresholds,
                plateau_fraction = plateau_fraction,
                pc_counts = pc_counts,
                pc_counts_per_job = pc_counts_per_job,
                pc_preemptible = pc_preemptible,
                docker_image = docker_image,
                prepare_cpu = prepare_cpu,
                prepare_memory_gb = prepare_memory_gb,
                prepare_disk_gb = prepare_disk_gb,
                analysis_cpu = analysis_cpu,
                analysis_memory_gb = analysis_memory_gb,
                analysis_disk_gb = analysis_disk_gb,
                max_retries = max_retries
        }

        OmicsResult result = object {
            name: dataset_name,
            phenotype_bed: phenotype_bed,
            principal_components_tsv: principal_components_tsv,
            additional_covariates_tsv: additional_covariates_tsv,
            selected_pc_z_scores_tsv_gz: RunMatrix.selected_pc_z_scores_tsv_gz,
            selected_pc_haplo_calls_tsv_gz: RunMatrix.selected_pc_haplo_calls_tsv_gz,
            selected_pc_z_scores_gene_qc_tsv_gz: RunMatrix.selected_pc_z_scores_gene_qc_tsv_gz,
            selected_pc_z_scores_summary_json: RunMatrix.selected_pc_z_scores_summary_json,
            results_tsv: RunMatrix.results_tsv,
            summary_json: RunMatrix.summary_json,
            gene_pc_qc_tsv_gz: RunMatrix.gene_pc_qc_tsv_gz,
            analysis_qc_json: RunMatrix.analysis_qc_json,
            pc_selection_json: RunMatrix.pc_selection_json,
            enrichment_plot_svg: RunMatrix.enrichment_plot_svg,
            pc_sweep_qc_summary_tsv: RunMatrix.pc_sweep_qc_summary_tsv,
            pc_sweep_qc_plot_png: RunMatrix.pc_sweep_qc_plot_png,
            protein_coding_genes_tsv: RunMatrix.protein_coding_genes_tsv,
            protein_coding_genes_qc_json: RunMatrix.protein_coding_genes_qc_json
        }
    }

    Array[File] expression_haplo_matrices = select_all(RunMatrix.selected_pc_haplo_calls_tsv_gz)
    Int calculated_intersection_disk_gb = ceil((size(RunMatrix.selected_pc_z_scores_tsv_gz, "GiB") + size(expression_haplo_matrices, "GiB")) * 8.0 + 20.0)
    Int dynamic_intersection_disk_gb = if calculated_intersection_disk_gb > intersection_disk_gb then calculated_intersection_disk_gb else intersection_disk_gb

    call IntersectMultiOmicsOutliers {
        input:
            z_score_matrices = RunMatrix.selected_pc_z_scores_tsv_gz,
            expression_haplo_matrices = expression_haplo_matrices,
            dataset_ids = dataset_name,
            thresholds = intersection_z_thresholds,
            docker_image = docker_image,
            cpu = intersection_cpu,
            memory_gb = intersection_memory_gb,
            disk_gb = dynamic_intersection_disk_gb,
            max_retries = max_retries
    }

    output {
        Array[OmicsResult] matrix_results = result
        Array[String] dataset_ids = dataset_name
        Array[File] z_score_matrices = RunMatrix.selected_pc_z_scores_tsv_gz
        Array[File?] haplo_matrices = RunMatrix.selected_pc_haplo_calls_tsv_gz
        Array[File] intersection_matrices = IntersectMultiOmicsOutliers.matrices
        File intersection_manifest_tsv = IntersectMultiOmicsOutliers.manifest_tsv
        File intersection_summary_json = IntersectMultiOmicsOutliers.summary_json
    }
}
