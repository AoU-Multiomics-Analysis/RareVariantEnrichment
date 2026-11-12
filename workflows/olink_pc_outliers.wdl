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
        String docker_image = "docker.io/epadhi/rarevariantenrichment-olink:root-no-entrypoint"
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
        Int collect_disk_gb = 20
        Int max_retries = 1
        Int preemptible = 0
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    Array[String] metrics = ["NPX", "PCNormalizedNPX"]
    Array[Pair[Int, Int]] adjustment_grid = cross(range(length(metrics)), pc_counts)

    call ValidateOlinkPcGrid {
        input:
            pc_counts = pc_counts,
            docker_image = docker_image,
            max_retries = max_retries,
            zones = zones
    }

    scatter (metric in metrics) {
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

    }

    # The first pair member selects the metric and its typed prepared File.
    # The second member is the requested PC count.
    scatter (adjustment in adjustment_grid) {
        call AdjustOlink {
            input:
                prepared_file = PrepareOlink.prepared_file[adjustment.left],
                metric = metrics[adjustment.left],
                n_pcs = adjustment.right,
                protein_cutoff = protein_cutoff,
                docker_image = docker_image,
                cpu = adjust_cpu,
                memory_gb = adjust_memory_gb,
                disk_gb = adjust_disk_gb,
                max_retries = max_retries,
                preemptible = preemptible,
                zones = zones
        }

    }

    call CollectOlinkSummaries {
        input:
            summary_files = AdjustOlink.summary,
            files_to_collect = flatten([flatten(AdjustOlink.collection_files), flatten(PrepareOlink.collection_files)]),
            adjustment_paths = AdjustOlink.collection_paths,
            preparation_paths = PrepareOlink.collection_paths,
            pc_plan_file = ValidateOlinkPcGrid.pc_plan,
            docker_image = docker_image,
            max_retries = max_retries,
            disk_gb = collect_disk_gb,
            zones = zones
    }

    # Assemble records from the collected File outputs; this starts no tasks.
    scatter (result_index in range(length(adjustment_grid))) {
        Pair[Int, Int] result_setting = adjustment_grid[result_index]
        OlinkPcResult result = object {
            npx_type: metrics[result_setting.left],
            pcs_removed: result_setting.right,
            z_scores: CollectOlinkSummaries.z_score_matrices[result_index],
            protein_outliers: CollectOlinkSummaries.protein_outlier_files[result_index],
            sample_outliers: CollectOlinkSummaries.sample_outlier_files[result_index],
            undefined_residual_assays: CollectOlinkSummaries.undefined_residual_assay_files[result_index],
            summary: CollectOlinkSummaries.run_summaries[result_index]
        }
    }

    output {
        Array[OlinkPcResult] results = result
        Array[File] z_score_matrices = CollectOlinkSummaries.z_score_matrices
        Array[File] protein_outlier_files = CollectOlinkSummaries.protein_outlier_files
        Array[File] sample_outlier_files = CollectOlinkSummaries.sample_outlier_files
        Array[File] undefined_residual_assay_files = CollectOlinkSummaries.undefined_residual_assay_files
        Array[File] run_summaries = CollectOlinkSummaries.run_summaries
        Array[File] sample_qc_files = CollectOlinkSummaries.sample_qc_files
        Array[File] excluded_assay_files = CollectOlinkSummaries.excluded_assay_files
        Array[File] preparation_summaries = CollectOlinkSummaries.preparation_summaries
        Array[File] session_info_files = CollectOlinkSummaries.session_info_files
        Array[File] collected_files = flatten([CollectOlinkSummaries.files, [CollectOlinkSummaries.pc_plan, CollectOlinkSummaries.summary, CollectOlinkSummaries.manifest]])
        File pc_plan = CollectOlinkSummaries.pc_plan
        File summary = CollectOlinkSummaries.summary
        File manifest = CollectOlinkSummaries.manifest
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

    Array[String] file_names = [
        "sample_qc.~{metric}.tsv",
        "excluded_assays.~{metric}.tsv",
        "prepare_summary.~{metric}.tsv",
        "sessionInfo.~{metric}.txt"
    ]

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
        File sample_qc = "prepared/" + file_names[0]
        File excluded_assays_file = "prepared/" + file_names[1]
        File summary = "prepared/" + file_names[2]
        File session_info = "prepared/" + file_names[3]
        Array[File] collection_files = prefix("prepared/", file_names)
        Array[String] collection_paths = prefix(metric + "/", file_names)
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

    Array[String] file_names = [
        "protein_z_scores.~{metric}.PCs~{n_pcs}.tsv.gz",
        "protein_outliers.~{metric}.PCs~{n_pcs}.tsv.gz",
        "sample_outliers.~{metric}.PCs~{n_pcs}.tsv",
        "undefined_residual_assays.~{metric}.PCs~{n_pcs}.tsv",
        "summary.~{metric}.PCs~{n_pcs}.tsv"
    ]

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
        File z_scores = "adjusted/" + file_names[0]
        File protein_outliers = "adjusted/" + file_names[1]
        File sample_outliers = "adjusted/" + file_names[2]
        File undefined_residual_assays = "adjusted/" + file_names[3]
        File summary = "adjusted/" + file_names[4]
        Array[File] collection_files = prefix("adjusted/", file_names)
        Array[String] collection_paths = prefix(metric + "/", file_names)
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
        Array[File] files_to_collect
        Array[Array[String]] adjustment_paths
        Array[Array[String]] preparation_paths
        File pc_plan_file
        String docker_image
        Int max_retries
        Int disk_gb
        String zones
    }

    Array[String] relative_paths = flatten([flatten(adjustment_paths), flatten(preparation_paths)])
    Int calculated_disk_gb = ceil(size(files_to_collect, "GiB") * 2.0 + 20.0)
    Int effective_disk_gb = if calculated_disk_gb > disk_gb then calculated_disk_gb else disk_gb

    command <<<
        set -euo pipefail
        echo "[collect] Starting final Olink file collection" >&2
        # File lists are written only after command-time localization.
        Rscript --vanilla - \
            --files-list '~{write_lines(files_to_collect)}' \
            --relative-paths-list '~{write_lines(relative_paths)}' \
            --summary-files '~{write_lines(summary_files)}' \
            --pc-plan-file '~{sub(pc_plan_file, "'", "'\"'\"'")}' \
            --output-dir collected <<'RSCRIPT'
        source("/opt/olink/olink_outliers.R")
        tryCatch({
          options <- parse_options(
            commandArgs(trailingOnly = TRUE),
            c("files-list", "relative-paths-list", "summary-files", "pc-plan-file", "output-dir"),
            c("files-list", "relative-paths-list", "summary-files", "pc-plan-file", "output-dir")
          )
          walk(options[c("files-list", "relative-paths-list", "summary-files", "pc-plan-file")],
               assert_readable)
          sources <- read_lines(options[["files-list"]])
          relative <- read_lines(options[["relative-paths-list"]])
          if (!length(sources) || length(sources) != length(relative)) {
            stop("Collection file and relative-path counts must match and be nonempty.", call. = FALSE)
          }
          walk(sources, assert_readable)
          if (anyNA(relative) || any(!str_detect(
            relative, "^(NPX|PCNormalizedNPX)/[A-Za-z0-9_][A-Za-z0-9_.-]*$"
          ))) {
            stop("Collection paths must be safe relative paths inside a measurement directory.",
                 call. = FALSE)
          }
          if (anyDuplicated(relative)) {
            stop("Collection paths must not contain duplicates.", call. = FALSE)
          }
          output_dir <- options[["output-dir"]]
          make_output_dir(output_dir)
          copy_one <- function(source, relative_path) {
            destination <- file.path(output_dir, relative_path)
            dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
            if (file.exists(destination)) {
              stop("Collection destination already exists: ", relative_path, call. = FALSE)
            }
            linked <- suppressWarnings(file.link(source, destination))
            if (!linked && !file.copy(source, destination, overwrite = FALSE)) {
              stop("Could not collect: ", relative_path, call. = FALSE)
            }
          }
          walk2(sources, relative, copy_one)
          copy_one(options[["pc-plan-file"]], "pc_plan.tsv")
          summaries <- collect_summaries(
            options[["summary-files"]], file.path(output_dir, "summary.tsv")
          )
          file_columns <- grep("_file$", names(summaries), value = TRUE)
          summaries <- summaries |>
            mutate(across(all_of(file_columns), ~ paste(NPX_type, .x, sep = "/")))
          write_tsv(summaries, file.path(output_dir, "summary.tsv"))
          manifest <- tibble(Relative_path = relative) |>
            mutate(
              Measurement = dirname(Relative_path),
              PCs_removed = as.integer(str_match(basename(Relative_path), "[.]PCs([0-9]+)[.]")[, 2]),
              Output_type = str_remove(basename(Relative_path), "[.](NPX|PCNormalizedNPX)[.].*$")
            ) |>
            bind_rows(tibble(
              Relative_path = c("pc_plan.tsv", "summary.tsv"),
              Measurement = "shared",
              PCs_removed = NA_integer_,
              Output_type = c("pc_plan", "summary")
            )) |>
            dplyr::select(Measurement, PCs_removed, Output_type, Relative_path) |>
            arrange(match(Measurement, c("NPX", "PCNormalizedNPX", "shared")),
                    PCs_removed, Output_type, Relative_path)
          write_tsv(manifest, file.path(output_dir, "manifest.tsv"), na = "")
          log_stage("collect", "Collected ", nrow(manifest), " files and wrote manifest.tsv.")
        }, error = function(error) {
          message("ERROR: ", conditionMessage(error))
          quit(save = "no", status = 1L)
        })
        RSCRIPT
        echo "[collect] Completed final Olink file collection" >&2
    >>>

    output {
        Array[File] files = prefix("collected/", relative_paths)
        Array[File] z_score_matrices = prefix("collected/", transpose(adjustment_paths)[0])
        Array[File] protein_outlier_files = prefix("collected/", transpose(adjustment_paths)[1])
        Array[File] sample_outlier_files = prefix("collected/", transpose(adjustment_paths)[2])
        Array[File] undefined_residual_assay_files = prefix("collected/", transpose(adjustment_paths)[3])
        Array[File] run_summaries = prefix("collected/", transpose(adjustment_paths)[4])
        Array[File] sample_qc_files = prefix("collected/", transpose(preparation_paths)[0])
        Array[File] excluded_assay_files = prefix("collected/", transpose(preparation_paths)[1])
        Array[File] preparation_summaries = prefix("collected/", transpose(preparation_paths)[2])
        Array[File] session_info_files = prefix("collected/", transpose(preparation_paths)[3])
        File pc_plan = "collected/pc_plan.tsv"
        File summary = "collected/summary.tsv"
        File manifest = "collected/manifest.tsv"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ~{effective_disk_gb} SSD"
        bootDiskSizeGb: 20
        maxRetries: max_retries
        zones: zones
    }
}
