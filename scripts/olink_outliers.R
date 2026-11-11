#!/usr/bin/env Rscript
# Analyze each NPX type once, then reuse its QC results and exact PCA fit.
suppressPackageStartupMessages(library(tidyverse))

log_stage <- function(metric, ...) message("[", metric, "] ", ...)

assert_local_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Supply one nonempty local path.", call. = FALSE)
  }
  if (grepl("^[[:alpha:]][[:alnum:]+.-]*://", path)) {
    stop("File localization error: unresolved URI: ", path,
         ". Supply this input as a WDL File.", call. = FALSE)
  }
  invisible(path)
}

assert_readable <- function(path) {
  assert_local_path(path)
  if (!file.exists(path) || dir.exists(path) || file.access(path, 4) != 0L) {
    stop("Input is not readable: ", path, call. = FALSE)
  }
  invisible(path)
}

make_output_dir <- function(path) {
  assert_local_path(path)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path) || file.access(path, 2) != 0L) {
    stop("Output directory is not writable: ", path, call. = FALSE)
  }
  invisible(path)
}

parse_integer <- function(value, name, minimum = 0L) {
  text <- as.character(value)
  if (length(text) != 1L || is.na(text) || !grepl("^[0-9]+$", text) ||
      !is.finite(as.numeric(text)) || as.numeric(text) > .Machine$integer.max ||
      as.numeric(text) < minimum) {
    stop(name, " must be one integer of at least ", minimum, ".", call. = FALSE)
  }
  as.integer(text)
}

parse_number <- function(value, name, positive = FALSE) {
  result <- suppressWarnings(as.numeric(value))
  if (length(result) != 1L || !is.finite(result) || (positive && result <= 0)) {
    stop(name, " must be one ", if (positive) "positive " else "",
         "finite number.", call. = FALSE)
  }
  result
}

parse_pc_counts <- function(value) {
  if (length(value) != 1L || is.na(value) || !nzchar(trimws(value)) ||
      grepl("(^|,)\\s*(,|$)", value)) {
    stop("PC counts must be a nonempty comma-separated integer list.", call. = FALSE)
  }
  fields <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  counts <- map_int(fields, parse_integer, name = "PC count")
  if (anyDuplicated(counts)) stop("PC counts must not contain duplicates.", call. = FALSE)
  counts
}

check_metric <- function(metric) {
  if (length(metric) != 1L || is.na(metric) || !metric %in% c("NPX", "PCNormalizedNPX")) {
    stop("metric must be NPX or PCNormalizedNPX.", call. = FALSE)
  }
  metric
}

pc_filename <- function(stem, metric, n_pcs, extension = ".tsv") {
  paste0(stem, ".", metric, ".PCs", n_pcs, extension)
}

read_exclusions <- function(path = NULL) {
  if (is.null(path)) return(c("GBP1", "MAP2K1"))
  assert_readable(path)
  assays <- readLines(path, warn = FALSE)
  if (any(!nzchar(trimws(assays)))) stop("Excluded assay names must not be empty.", call. = FALSE)
  unique(assays)
}

read_npx <- function(path, metric) {
  assert_readable(path)
  required <- c("Assay", "SampleID", metric)
  header <- names(data.table::fread(path, nrows = 0L, data.table = FALSE))
  missing <- setdiff(required, header)
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(header)) stop("Input column names must be unique.", call. = FALSE)
  data.table::fread(path, select = required,
                   colClasses = list(character = c("Assay", "SampleID")),
                   data.table = FALSE) |>
    as_tibble()
}

drop_constant_assays <- function(values, stage) {
  assay_sd <- apply(values, 2L, stats::sd)
  constant <- !is.finite(assay_sd) | assay_sd == 0
  list(values = values[, !constant, drop = FALSE],
       excluded = tibble(Assay = colnames(values)[constant], Stage = stage))
}

check_pc_capacity <- function(values, max_pcs, metric) {
  capacity <- min(nrow(values) - 2L, ncol(values) - 1L)
  if (max_pcs > capacity) {
    stop(metric, ": requested ", max_pcs, " PCs, but at most ", capacity,
         " can be removed while leaving residual variation. Samples: ", nrow(values),
         "; assays: ", ncol(values), ".", call. = FALSE)
  }
}

sample_connectivity <- function(scaled_values, block_size) {
  profiles <- t(scaled_values)
  blocks <- split(seq_len(ncol(profiles)), ceiling(seq_len(ncol(profiles)) / block_size))
  map_dfr(blocks, function(index) {
    adjacency <- 0.5 + 0.5 * WGCNA::bicor(
      profiles, profiles[, index, drop = FALSE], use = "pairwise.complete.obs"
    )
    if (anyNA(adjacency)) {
      warning("Undefined bicor values were set to zero, as in the original method.", call. = FALSE)
      adjacency[is.na(adjacency)] <- 0
    }
    tibble(SampleID = colnames(profiles)[index],
           Connectivity = colSums(adjacency) - adjacency[cbind(index, seq_along(index))])
  })
}

prepare_analysis <- function(data, metric, max_pcs, excluded_assays = c("GBP1", "MAP2K1"),
                             connectivity_cutoff = -3, block_size = 500L,
                             threads = 8L, strip_sample_suffix = TRUE) {
  metric <- check_metric(metric)
  max_pcs <- parse_integer(max_pcs, "maximum PC count")
  block_size <- parse_integer(block_size, "block size", 1L)
  threads <- parse_integer(threads, "threads", 1L)
  connectivity_cutoff <- parse_number(connectivity_cutoff, "connectivity cutoff")
  if (length(strip_sample_suffix) != 1L || is.na(strip_sample_suffix) ||
      !is.logical(strip_sample_suffix)) stop("strip_sample_suffix must be TRUE or FALSE.", call. = FALSE)
  required <- c("Assay", "SampleID", metric)
  if (!all(required %in% names(data))) {
    stop("Missing columns: ", paste(setdiff(required, names(data)), collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(data[[metric]]) &&
      !(is.logical(data[[metric]]) && all(is.na(data[[metric]])))) {
    stop(metric, " must contain numeric values.", call. = FALSE)
  }
  if (anyNA(excluded_assays) || any(!nzchar(trimws(excluded_assays)))) {
    stop("Excluded assay names must not be empty.", call. = FALSE)
  }
  log_stage(metric, "Prepare assay values; maximum PCs = ", max_pcs, ".")
  selected <- data |>
    dplyr::select(all_of(required)) |>
    mutate(Assay = as.character(Assay), RawSampleID = as.character(SampleID),
           SampleID = if (strip_sample_suffix) str_remove(RawSampleID, "_.*$") else RawSampleID,
           Value = as.numeric(.data[[metric]])) |>
    dplyr::select(Assay, SampleID, RawSampleID, Value)
  configured <- selected |>
    filter(Assay %in% excluded_assays) |>
    distinct(Assay) |>
    mutate(Stage = "configured_exclusion")
  selected <- selected |> filter(!Assay %in% excluded_assays)
  if (!nrow(selected)) stop("No assays remain after the configured exclusions.", call. = FALSE)
  if (anyNA(selected$Assay) || anyNA(selected$SampleID) ||
      any(!nzchar(trimws(selected$Assay))) || any(!nzchar(trimws(selected$SampleID)))) {
    stop("Input must contain nonempty assay names and sample IDs.", call. = FALSE)
  }
  collisions <- selected |>
    distinct(SampleID, RawSampleID) |> count(SampleID) |> filter(n > 1L)
  if (nrow(collisions)) {
    stop("Sample ID suffix removal would merge distinct input IDs: ", collisions$SampleID[1],
         ". Resolve these IDs before analysis.", call. = FALSE)
  }
  duplicates <- selected |> count(Assay, SampleID) |> filter(n > 1L)
  if (nrow(duplicates)) {
    stop("Duplicate assay/sample keys after sample ID normalization: ",
         duplicates$Assay[1], " / ", duplicates$SampleID[1], ".", call. = FALSE)
  }
  values <- selected |>
    dplyr::select(Assay, SampleID, Value) |>
    pivot_wider(names_from = Assay, values_from = Value) |>
    column_to_rownames("SampleID") |>
    as.matrix()
  sample_qc <- tibble(SampleID = rownames(values),
                      Missing_or_nonfinite = rowSums(!is.finite(values)))
  values <- values[sample_qc$Missing_or_nonfinite == 0L, , drop = FALSE]
  if (nrow(values) < 3L) stop(metric, ": fewer than three samples have complete finite values.", call. = FALSE)
  before <- drop_constant_assays(values, "before_connectivity")
  values <- before$values
  if (ncol(values) < 2L) stop("Fewer than two variable assays remain.", call. = FALSE)
  check_pc_capacity(values, max_pcs, metric)
  log_stage(metric, "Calculate sample connectivity for ", nrow(values), " samples and ",
            ncol(values), " assays; block size = ", block_size, "; threads = ", threads, ".")
  if (threads == 1L) WGCNA::disableWGCNAThreads() else WGCNA::enableWGCNAThreads(nThreads = threads)
  connectivity <- sample_connectivity(scale(values), block_size)
  connectivity_sd <- stats::sd(connectivity$Connectivity)
  if (!is.finite(connectivity_sd) || connectivity_sd == 0) {
    stop(metric, ": connectivity has zero or undefined variance.", call. = FALSE)
  }
  connectivity <- connectivity |>
    mutate(Connectivity_Z = (Connectivity - mean(Connectivity)) / connectivity_sd)
  sample_qc <- sample_qc |>
    left_join(connectivity, by = "SampleID") |>
    mutate(NPX_type = metric, Max_PCs_prepared = max_pcs,
           Status = case_when(Missing_or_nonfinite > 0L ~ "missing_values",
                              Connectivity_Z < connectivity_cutoff ~ "connectivity_outlier",
                              TRUE ~ "kept"))
  kept <- sample_qc |> filter(Status == "kept") |> pull(SampleID)
  values <- values[kept, , drop = FALSE]
  if (nrow(values) < 3L) stop("Fewer than three samples remain after QC.", call. = FALSE)
  after <- drop_constant_assays(values, "after_connectivity")
  values <- after$values
  if (ncol(values) < 2L) stop("Fewer than two variable assays remain after QC.", call. = FALSE)
  check_pc_capacity(values, max_pcs, metric)
  scaled <- scale(values)
  scores <- matrix(numeric(), nrow = nrow(scaled), ncol = 0L,
                   dimnames = list(rownames(scaled), NULL))
  loadings <- matrix(numeric(), nrow = ncol(scaled), ncol = 0L,
                     dimnames = list(colnames(scaled), NULL))
  effective_rank <- NA_integer_
  if (max_pcs > 0L) {
    log_stage(metric, "Fit exact PCA once through PC ", max_pcs, ".")
    pca <- stats::prcomp(scaled, center = FALSE, scale. = FALSE, rank. = max_pcs)
    effective_rank <- sum(pca$sdev > max(pca$sdev) * sqrt(.Machine$double.eps))
    if (max_pcs >= effective_rank) {
      stop(metric, ": requested ", max_pcs, " PCs; effective data rank is ", effective_rank,
           ". Leave at least one residual component.", call. = FALSE)
    }
    scores <- pca$x
    loadings <- pca$rotation
  }
  excluded <- bind_rows(configured, before$excluded, after$excluded) |>
    mutate(NPX_type = metric, Max_PCs_prepared = max_pcs)
  summary <- tibble(
    NPX_type = metric, Max_PCs_prepared = max_pcs,
    Connectivity_cutoff = connectivity_cutoff, Strip_sample_suffix = strip_sample_suffix,
    Samples_input = nrow(sample_qc), Samples_missing = sum(sample_qc$Status == "missing_values"),
    Sample_outliers = sum(sample_qc$Status == "connectivity_outlier"),
    Samples_kept = nrow(scaled), Assays_kept = ncol(scaled), Effective_rank = effective_rank,
    Sample_qc_file = paste0("sample_qc.", metric, ".tsv"),
    Excluded_assays_file = paste0("excluded_assays.", metric, ".tsv"),
    Session_info_file = paste0("sessionInfo.", metric, ".txt")
  )
  log_stage(metric, "QC retained ", nrow(scaled), " samples and ", ncol(scaled), " assays.")
  list(format_version = 1L, scaled_values = scaled, scores = scores, loadings = loadings,
       sample_qc = sample_qc, excluded_assays = excluded, prepare_summary = summary,
       metadata = list(metric = metric, max_pcs = max_pcs, effective_rank = effective_rank,
                       connectivity_cutoff = connectivity_cutoff, strip_sample_suffix = strip_sample_suffix))
}

adjust_analysis <- function(state, n_pcs, protein_cutoff = 3, metric = NULL) {
  n_pcs <- parse_integer(n_pcs, "PC count")
  protein_cutoff <- parse_number(protein_cutoff, "protein cutoff", positive = TRUE)
  if (!is.list(state) || !identical(state$format_version, 1L)) {
    stop("Unsupported prepared state.", call. = FALSE)
  }
  if (!is.null(metric) && !identical(check_metric(metric), state$metadata$metric)) {
    stop("Requested metric does not match the prepared state.", call. = FALSE)
  }
  metric <- check_metric(state$metadata$metric)
  if (n_pcs > state$metadata$max_pcs) stop("PC count exceeds the prepared maximum.", call. = FALSE)
  log_stage(metric, "Remove ", n_pcs, " PCs from the prepared fit.")
  residual <- state$scaled_values
  if (n_pcs > 0L) {
    index <- seq_len(n_pcs)
    residual <- residual - tcrossprod(state$scores[, index, drop = FALSE],
                                     state$loadings[, index, drop = FALSE])
  }
  residual <- sweep(residual, 2L, colMeans(residual), "-")
  residual_sd <- sqrt(colMeans(residual^2))
  undefined <- !is.finite(residual_sd) | residual_sd <= sqrt(.Machine$double.eps)
  residual_sd[undefined] <- NA_real_
  z_scores <- sweep(residual, 2L, residual_sd, "/")
  z_scores[!is.finite(z_scores)] <- NA_real_
  index <- which(abs(z_scores) > protein_cutoff, arr.ind = TRUE)
  outliers <- tibble(NPX_type = rep(metric, nrow(index)),
                     Assay = colnames(z_scores)[index[, 2L]],
                     SampleID = rownames(z_scores)[index[, 1L]], Z_score = z_scores[index]) |>
    mutate(Direction = if_else(Z_score < 0, "low", "high"),
           underlier = Z_score < -protein_cutoff, overlier = Z_score > protein_cutoff,
           PCs_removed = n_pcs) |>
    arrange(Assay, desc(abs(Z_score)))
  sample_outliers <- state$sample_qc |>
    filter(Status == "connectivity_outlier") |>
    mutate(PCs_removed = n_pcs)
  undefined_assays <- tibble(NPX_type = rep(metric, sum(undefined)),
                             Assay = colnames(z_scores)[undefined], PCs_removed = n_pcs)
  summary <- state$prepare_summary |>
    mutate(PCs_removed = n_pcs, Protein_cutoff = protein_cutoff,
           Assays_undefined_residual_z = sum(undefined),
           Protein_low_outliers = sum(outliers$underlier), Protein_high_outliers = sum(outliers$overlier),
           Protein_z_scores_file = pc_filename("protein_z_scores", metric, n_pcs, ".tsv.gz"),
           Protein_outliers_file = pc_filename("protein_outliers", metric, n_pcs, ".tsv.gz"),
           Sample_outliers_file = pc_filename("sample_outliers", metric, n_pcs),
           Undefined_residual_assays_file = pc_filename("undefined_residual_assays", metric, n_pcs),
           Summary_file = pc_filename("summary", metric, n_pcs)) |>
    relocate(NPX_type, PCs_removed)
  log_stage(metric, "PCs ", n_pcs, ": found ", nrow(outliers), " protein outlier measurements; ",
            sum(undefined), " assays have undefined residual Z scores.")
  list(z_scores = z_scores, protein_outliers = outliers, sample_outliers = sample_outliers,
       undefined_residual_assays = undefined_assays, summary = summary)
}

write_prepared <- function(state, output_dir) {
  make_output_dir(output_dir)
  metric <- state$metadata$metric
  saveRDS(state, file.path(output_dir, "prepared.rds"))
  write_tsv(state$sample_qc, file.path(output_dir, paste0("sample_qc.", metric, ".tsv")))
  write_tsv(state$excluded_assays, file.path(output_dir, paste0("excluded_assays.", metric, ".tsv")))
  write_tsv(state$prepare_summary, file.path(output_dir, paste0("prepare_summary.", metric, ".tsv")))
  capture.output(sessionInfo(), file = file.path(output_dir, paste0("sessionInfo.", metric, ".txt")))
}

write_adjusted <- function(result, output_dir) {
  make_output_dir(output_dir)
  metric <- result$summary$NPX_type
  n_pcs <- result$summary$PCs_removed
  path <- function(stem, extension = ".tsv") file.path(output_dir, pc_filename(stem, metric, n_pcs, extension))
  t(result$z_scores) |>
    as_tibble(rownames = "Assay", .name_repair = "minimal") |>
    write_tsv(path("protein_z_scores", ".tsv.gz"))
  write_tsv(result$protein_outliers, path("protein_outliers", ".tsv.gz"))
  write_tsv(result$sample_outliers, path("sample_outliers"))
  write_tsv(result$undefined_residual_assays, path("undefined_residual_assays"))
  write_tsv(result$summary, path("summary"))
}

collect_summaries <- function(summary_files, output) {
  assert_readable(summary_files)
  paths <- readLines(summary_files, warn = FALSE)
  if (!length(paths) || any(!nzchar(paths))) stop("Summary file list must not be empty.", call. = FALSE)
  walk(paths, assert_readable)
  summaries <- map_dfr(paths, function(path) {
    row <- read_tsv(path, show_col_types = FALSE)
    if (nrow(row) != 1L || !all(c("NPX_type", "PCs_removed") %in% names(row))) {
      stop("Each summary file must contain one run with NPX_type and PCs_removed.", call. = FALSE)
    }
    check_metric(row$NPX_type)
    parse_integer(row$PCs_removed, "summary PC count")
    row
  })
  if (anyDuplicated(summaries[c("NPX_type", "PCs_removed")])) {
    stop("Summary files contain a duplicate metric/PC run.", call. = FALSE)
  }
  metadata_files <- grep("_file$", names(summaries), value = TRUE)
  if (any(grepl("[/\\\\]", unlist(summaries[metadata_files])))) {
    stop("Summary file metadata must contain filenames, without directory paths.", call. = FALSE)
  }
  assert_local_path(output)
  make_output_dir(dirname(output))
  write_tsv(summaries, output)
  log_stage("collect", "Collected ", nrow(summaries), " metric/PC summaries.")
  invisible(summaries)
}

parse_options <- function(args, allowed, required) {
  if (length(args) %% 2L != 0L) stop("Each named CLI argument needs one value.", call. = FALSE)
  keys <- args[seq.int(1L, length(args), by = 2L)]
  values <- args[seq.int(2L, length(args), by = 2L)]
  if (any(!grepl("^--", keys)) || anyDuplicated(keys) || any(!keys %in% paste0("--", allowed))) {
    stop("CLI arguments must use unique supported --names.", call. = FALSE)
  }
  options <- as.list(setNames(values, sub("^--", "", keys)))
  missing <- setdiff(required, names(options))
  if (length(missing)) stop("Missing CLI arguments: ", paste(paste0("--", missing), collapse = ", "), call. = FALSE)
  options
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || identical(args, "--help")) {
    cat("Usage: Rscript olink_outliers.R <validate|prepare|adjust|collect> --name value ...\n")
    return(invisible(NULL))
  }
  command <- args[1L]
  args <- args[-1L]
  if (command == "validate") {
    options <- parse_options(args, c("pc-counts", "output-dir"), c("pc-counts", "output-dir"))
    counts <- parse_pc_counts(options[["pc-counts"]])
    make_output_dir(options[["output-dir"]])
    writeLines(as.character(max(counts)), file.path(options[["output-dir"]], "max_pcs.txt"))
    write_tsv(tibble(PCs_removed = counts), file.path(options[["output-dir"]], "pc_plan.tsv"))
    log_stage("validate", "Requested PC counts: ", paste(counts, collapse = ", "), ".")
  } else if (command == "prepare") {
    allowed <- c("input", "metric", "max-pcs", "excluded-assays-file", "connectivity-cutoff",
                 "block-size", "threads", "strip-sample-suffix", "output-dir")
    options <- parse_options(args, allowed, c("input", "metric", "max-pcs", "output-dir"))
    metric <- check_metric(options$metric)
    excluded <- read_exclusions(options[["excluded-assays-file"]])
    strip <- options[["strip-sample-suffix"]] %||% "true"
    if (!strip %in% c("true", "false")) stop("--strip-sample-suffix must be true or false.", call. = FALSE)
    state <- prepare_analysis(read_npx(options$input, metric), metric, options[["max-pcs"]],
                              excluded, options[["connectivity-cutoff"]] %||% -3,
                              options[["block-size"]] %||% 500L, options$threads %||% 8L,
                              strip == "true")
    write_prepared(state, options[["output-dir"]])
  } else if (command == "adjust") {
    options <- parse_options(args, c("prepared-file", "n-pcs", "metric", "protein-cutoff", "output-dir"),
                              c("prepared-file", "n-pcs", "output-dir"))
    assert_readable(options[["prepared-file"]])
    state <- readRDS(options[["prepared-file"]])
    result <- adjust_analysis(state, options[["n-pcs"]], options[["protein-cutoff"]] %||% 3, options$metric)
    write_adjusted(result, options[["output-dir"]])
  } else if (command == "collect") {
    options <- parse_options(args, c("summary-files", "output"), c("summary-files", "output"))
    collect_summaries(options[["summary-files"]], options$output)
  } else {
    stop("Unknown command: ", command, call. = FALSE)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  tryCatch(main(), error = function(error) {
    message("ERROR: ", conditionMessage(error))
    quit(save = "no", status = 1L)
  })
}
