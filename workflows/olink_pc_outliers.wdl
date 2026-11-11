version 1.0

struct OlinkPcResult {
    String npx_type
    Int pcs_removed
    File z_scores
    File protein_outliers
    File sample_outliers
    File undefined_residual_assays
    File summary
}

workflow OlinkPcOutliers {
    input {
        File npx_file
        Array[Int] pc_counts
        # Override when the workspace needs a different image registry or tag.
        String docker_image = "us-central1-docker.pkg.dev/wb-perky-tomato-1843/olink-images/rarevariantenrichment-olink:main"
        Array[String] excluded_assays = ["GBP1", "MAP2K1"]
        Boolean strip_sample_suffix = true
        Float connectivity_cutoff = -3.0
        Float protein_cutoff = 3.0
        Int block_size = 500
        Int prepare_cpu = 8
        Int prepare_memory_gb = 128
        Int prepare_disk_gb = 500
        Int adjust_cpu = 8
        Int adjust_memory_gb = 64
        Int adjust_disk_gb = 200
        Int max_retries = 1
        Int preemptible = 0
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    call ValidateOlinkPcGrid {
        input:
            pc_counts = pc_counts,
            docker_image = docker_image,
            max_retries = max_retries,
            zones = zones
    }

    scatter (metric in ["NPX", "PCNormalizedNPX"]) {
        call PrepareOlink {
            input:
                npx_file = npx_file,
                metric = metric,
                max_pcs = ValidateOlinkPcGrid.max_pcs,
                excluded_assays = excluded_assays,
                strip_sample_suffix = strip_sample_suffix,
                connectivity_cutoff = connectivity_cutoff,
                block_size = block_size,
                docker_image = docker_image,
                cpu = prepare_cpu,
                memory_gb = prepare_memory_gb,
                disk_gb = prepare_disk_gb,
                max_retries = max_retries,
                preemptible = preemptible,
                zones = zones
        }

        scatter (pc_count in pc_counts) {
            call AdjustOlink {
                input:
                    prepared_file = PrepareOlink.prepared_file,
                    metric = metric,
                    n_pcs = pc_count,
                    protein_cutoff = protein_cutoff,
                    docker_image = docker_image,
                    cpu = adjust_cpu,
                    memory_gb = adjust_memory_gb,
                    disk_gb = adjust_disk_gb,
                    max_retries = max_retries,
                    preemptible = preemptible,
                    zones = zones
            }

            OlinkPcResult result = object {
                npx_type: metric,
                pcs_removed: pc_count,
                z_scores: AdjustOlink.z_scores,
                protein_outliers: AdjustOlink.protein_outliers,
                sample_outliers: AdjustOlink.sample_outliers,
                undefined_residual_assays: AdjustOlink.undefined_residual_assays,
                summary: AdjustOlink.summary
            }
        }
    }

    call CollectOlinkSummaries {
        input:
            summary_files = flatten(AdjustOlink.summary),
            docker_image = docker_image,
            max_retries = max_retries,
            zones = zones
    }

    output {
        Array[OlinkPcResult] results = flatten(result)
        Array[File] z_score_matrices = flatten(AdjustOlink.z_scores)
        Array[File] protein_outlier_files = flatten(AdjustOlink.protein_outliers)
        Array[File] sample_outlier_files = flatten(AdjustOlink.sample_outliers)
        Array[File] undefined_residual_assay_files = flatten(AdjustOlink.undefined_residual_assays)
        Array[File] run_summaries = flatten(AdjustOlink.summary)
        Array[File] sample_qc_files = PrepareOlink.sample_qc
        Array[File] excluded_assay_files = PrepareOlink.excluded_assays_file
        Array[File] preparation_summaries = PrepareOlink.summary
        Array[File] session_info_files = PrepareOlink.session_info
        File pc_plan = ValidateOlinkPcGrid.pc_plan
        File summary = CollectOlinkSummaries.summary
    }

    meta {
        description: "Generate Olink residual Z-score matrices and outlier calls for several PC counts in All of Us Workbench or Terra."
    }
}

task ValidateOlinkPcGrid {
    input {
        Array[Int] pc_counts
        String docker_image
        Int max_retries
        String zones
    }

    command <<<
        set -euo pipefail
        echo "[validate] Starting Olink PC-grid validation" >&2
        Rscript /opt/olink/olink_outliers.R validate \
            --pc-counts '~{sep="," pc_counts}' \
            --output-dir validation
        echo "[validate] Completed Olink PC-grid validation" >&2
    >>>

    output {
        Int max_pcs = read_int("validation/max_pcs.txt")
        File pc_plan = "validation/pc_plan.tsv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 SSD"
        bootDiskSizeGb: 20
        maxRetries: max_retries
        zones: zones
    }
}

task PrepareOlink {
    input {
        File npx_file
        String metric
        Int max_pcs
        Array[String] excluded_assays
        Boolean strip_sample_suffix
        Float connectivity_cutoff
        Int block_size
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
        Int preemptible
        String zones
    }

    # Account for input decompression and intermediate dense matrices.
    Int calculated_disk_gb = ceil(size(npx_file, "GiB") * 12.0 + 50.0)
    Int effective_disk_gb = if calculated_disk_gb > disk_gb then calculated_disk_gb else disk_gb

    command <<<
        set -euo pipefail
        echo "[prepare] Starting Olink sample QC and shared PCA" >&2
        export OMP_NUM_THREADS='~{cpu}'
        export OPENBLAS_NUM_THREADS='~{cpu}'
        export MKL_NUM_THREADS='~{cpu}'
        Rscript /opt/olink/olink_outliers.R prepare \
            --input '~{sub(npx_file, "'", "'\"'\"'")}' \
            --metric '~{sub(metric, "'", "'\"'\"'")}' \
            --max-pcs '~{max_pcs}' \
            --excluded-assays-file '~{write_lines(excluded_assays)}' \
            --strip-sample-suffix '~{true="true" false="false" strip_sample_suffix}' \
            --connectivity-cutoff '~{connectivity_cutoff}' \
            --block-size '~{block_size}' \
            --threads '~{cpu}' \
            --output-dir prepared
        echo "[prepare] Completed Olink sample QC and shared PCA" >&2
    >>>

    output {
        File prepared_file = "prepared/prepared.rds"
        File sample_qc = "prepared/sample_qc." + metric + ".tsv"
        File excluded_assays_file = "prepared/excluded_assays." + metric + ".tsv"
        File summary = "prepared/prepare_summary." + metric + ".tsv"
        File session_info = "prepared/sessionInfo." + metric + ".txt"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{effective_disk_gb} SSD"
        bootDiskSizeGb: 20
        maxRetries: max_retries
        preemptible: preemptible
        zones: zones
    }
}

task AdjustOlink {
    input {
        File prepared_file
        String metric
        Int n_pcs
        Float protein_cutoff
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_gb
        Int max_retries
        Int preemptible
        String zones
    }

    Int calculated_disk_gb = ceil(size(prepared_file, "GiB") * 6.0 + 20.0)
    Int effective_disk_gb = if calculated_disk_gb > disk_gb then calculated_disk_gb else disk_gb

    command <<<
        set -euo pipefail
        echo "[adjust] Starting Olink PC adjustment and outlier calls" >&2
        export OMP_NUM_THREADS='~{cpu}'
        export OPENBLAS_NUM_THREADS='~{cpu}'
        export MKL_NUM_THREADS='~{cpu}'
        Rscript /opt/olink/olink_outliers.R adjust \
            --prepared-file '~{sub(prepared_file, "'", "'\"'\"'")}' \
            --metric '~{sub(metric, "'", "'\"'\"'")}' \
            --n-pcs '~{n_pcs}' \
            --protein-cutoff '~{protein_cutoff}' \
            --output-dir adjusted
        echo "[adjust] Completed Olink PC adjustment and outlier calls" >&2
    >>>

    output {
        File z_scores = "adjusted/protein_z_scores.~{metric}.PCs~{n_pcs}.tsv.gz"
        File protein_outliers = "adjusted/protein_outliers.~{metric}.PCs~{n_pcs}.tsv.gz"
        File sample_outliers = "adjusted/sample_outliers.~{metric}.PCs~{n_pcs}.tsv"
        File undefined_residual_assays = "adjusted/undefined_residual_assays.~{metric}.PCs~{n_pcs}.tsv"
        File summary = "adjusted/summary.~{metric}.PCs~{n_pcs}.tsv"
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{effective_disk_gb} SSD"
        bootDiskSizeGb: 20
        maxRetries: max_retries
        preemptible: preemptible
        zones: zones
    }
}

task CollectOlinkSummaries {
    input {
        Array[File] summary_files
        String docker_image
        Int max_retries
        String zones
    }

    command <<<
        set -euo pipefail
        echo "[collect] Starting Olink result summary collection" >&2
        # Write localized paths only when the command is rendered.
        Rscript /opt/olink/olink_outliers.R collect \
            --summary-files '~{write_lines(summary_files)}' \
            --output summary.tsv
        echo "[collect] Completed Olink result summary collection" >&2
    >>>

    output {
        File summary = "summary.tsv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 SSD"
        bootDiskSizeGb: 20
        maxRetries: max_retries
        zones: zones
    }
}
