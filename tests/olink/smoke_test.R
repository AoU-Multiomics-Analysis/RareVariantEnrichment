#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(tidyverse))

script <- Sys.getenv("OLINK_SCRIPT", "scripts/olink_outliers.R")
this_script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
fixture_script <- file.path(dirname(this_script), "make_fixture.R")
source(fixture_script)
work <- tempfile("olink space's ")
dir.create(work)
run_cli <- function(...) {
  args <- c(script, ...)
  result <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                    args = vapply(args, shQuote, character(1)),
                                    stdout = TRUE, stderr = TRUE))
  list(status = attr(result, "status") %||% 0L, text = paste(result, collapse = "\n"))
}

test_that("validate preserves the requested PC order and writes a typed plan", {
  result <- run_cli("validate", "--pc-counts", "2,0,4", "--output-dir", file.path(work, "validation"))
  expect_equal(result$status, 0L, info = result$text)
  expect_equal(readLines(file.path(work, "validation", "max_pcs.txt")), "4")
  expect_equal(readr::read_tsv(file.path(work, "validation", "pc_plan.tsv"), show_col_types = FALSE)$PCs_removed,
               c(2, 0, 4))
})

source(script)
fixture <- make_olink_fixture()
input <- file.path(work, "input's data.tsv")
data.table::fwrite(fixture, input, sep = "\t", na = "NA")
exclusions <- file.path(work, "excluded assays.txt")
writeLines(c("GBP1", "MAP2K1"), exclusions)

# Independent, dense reference: full sample adjacency and QR regression on PCs.
# It does not call any production analysis helper.
reference_result <- function(data, metric, pcs, cutoff = -3) {
  data <- data[!data$Assay %in% c("GBP1", "MAP2K1"), ]
  sample <- sub("_.*$", "", data$SampleID)
  sample_names <- unique(sample)
  assay_names <- unique(data$Assay)
  values <- matrix(NA_real_, length(sample_names), length(assay_names),
                   dimnames = list(sample_names, assay_names))
  values[cbind(match(sample, sample_names), match(data$Assay, assay_names))] <- data[[metric]]
  values <- values[rowSums(!is.finite(values)) == 0L, , drop = FALSE]
  values <- values[, apply(values, 2, sd) > 0, drop = FALSE]
  adjacency <- 0.5 + 0.5 * WGCNA::bicor(t(scale(values)), use = "pairwise.complete.obs")
  adjacency[is.na(adjacency)] <- 0
  connectivity <- colSums(adjacency) - diag(adjacency)
  z <- (connectivity - mean(connectivity)) / sd(connectivity)
  values <- values[z >= cutoff, , drop = FALSE]
  values <- values[, apply(values, 2, sd) > 0, drop = FALSE]
  scaled <- scale(values)
  if (pcs > 0) {
    pca <- prcomp(scaled, center = FALSE, scale. = FALSE)
    design <- cbind(1, pca$x[, seq_len(pcs), drop = FALSE])
    residual <- qr.resid(qr(design), scaled)
  } else {
    residual <- sweep(scaled, 2, colMeans(scaled), "-")
  }
  residual <- sweep(residual, 2, colMeans(residual), "-")
  population_sd <- sqrt(colMeans(residual^2))
  population_sd[population_sd <= sqrt(.Machine$double.eps)] <- NA_real_
  list(z = sweep(residual, 2, population_sd, "/"), connectivity = connectivity)
}

test_that("PC grids reject empty, repeated, fractional, negative and malformed entries", {
  for (bad in c("", "1,", ",1", "1,,2", "1,1", "-1,2", "1.5,2", "x,2", "1e2", "2147483648")) {
    expect_error(parse_pc_counts(bad), "PC|pc")
  }
  expect_equal(parse_pc_counts(" 2,0,4 "), c(2L, 0L, 4L))
})

prepared <- list()
for (metric in c("NPX", "PCNormalizedNPX")) {
  test_that(paste(metric, "prepare and all PC prefixes match dense connectivity and QR residual Z scores"), {
    directory <- file.path(work, paste0("prepared ", metric))
    result <- run_cli("prepare", "--input", input, "--metric", metric,
                      "--max-pcs", "4", "--excluded-assays-file", exclusions,
                      "--connectivity-cutoff", "-3", "--block-size", "7",
                      "--threads", "1", "--strip-sample-suffix", "true",
                      "--output-dir", directory)
    expect_equal(result$status, 0L, info = result$text)
    state <- readRDS(file.path(directory, "prepared.rds"))
    prepared[[metric]] <<- state
    expect_equal(ncol(state$scores), 4L)
    expect_equal(ncol(state$loadings), 4L)
    expect_setequal(state$excluded_assays$Assay, c("GBP1", "MAP2K1", "Constant"))
    expect_equal(state$sample_qc$SampleID[state$sample_qc$Status == "missing_values"],
                 if (metric == "NPX") "001" else "002")
    reference <- reference_result(fixture, metric, 0)
    observed <- state$sample_qc |>
      filter(Status != "missing_values") |>
      dplyr::select(SampleID, Connectivity)
    expect_equal(unname(observed$Connectivity), unname(reference$connectivity), tolerance = 1e-10)
    expect_equal(sd(state$sample_qc$Connectivity_Z, na.rm = TRUE), 1)
    for (pcs in c(0L, 2L, 4L)) {
      adjusted <- adjust_analysis(state, pcs, 3, metric)
      reference <- reference_result(fixture, metric, pcs)
      expect_equal(unname(adjusted$z_scores), unname(reference$z), tolerance = 1e-9)
      expect_equal(dimnames(adjusted$z_scores), dimnames(reference$z))
      expect_true(all(adjusted$protein_outliers$PCs_removed == pcs))
      expect_true(all(adjusted$protein_outliers$NPX_type == metric))
      expect_equal(adjusted$summary$PCs_removed, pcs)
      expect_equal(adjusted$summary$NPX_type, metric)
      calls <- which(abs(reference$z) > 3, arr.ind = TRUE)
      expect_setequal(paste(adjusted$protein_outliers$Assay, adjusted$protein_outliers$SampleID),
                      paste(colnames(reference$z)[calls[, 2]], rownames(reference$z)[calls[, 1]]))
      expect_equal(adjusted$summary$Protein_low_outliers, sum(reference$z < -3, na.rm = TRUE))
      expect_equal(adjusted$summary$Protein_high_outliers, sum(reference$z > 3, na.rm = TRUE))
      expect_true(all(adjusted$protein_outliers$underlier == (adjusted$protein_outliers$Z_score < -3)))
      expect_true(all(adjusted$protein_outliers$overlier == (adjusted$protein_outliers$Z_score > 3)))
    }
    adjusted_dir <- file.path(work, paste0("adjusted ", metric))
    result <- run_cli("adjust", "--prepared-file", file.path(directory, "prepared.rds"),
                      "--metric", metric, "--n-pcs", "2", "--protein-cutoff", "3",
                      "--output-dir", adjusted_dir)
    expect_equal(result$status, 0L, info = result$text)
    expected <- paste0(c("protein_z_scores", "protein_outliers", "sample_outliers",
                         "undefined_residual_assays", "summary"), ".", metric, ".PCs2",
                       c(".tsv.gz", ".tsv.gz", ".tsv", ".tsv", ".tsv"))
    expect_setequal(list.files(adjusted_dir), expected)
    written <- readr::read_tsv(file.path(adjusted_dir, expected[1]), show_col_types = FALSE,
                               name_repair = "minimal")
    expect_identical(names(written)[1], "Assay")
    expect_identical(names(written)[-1], rownames(reference_result(fixture, metric, 2)$z))
    expected_z <- t(reference_result(fixture, metric, 2)$z)
    expect_equal(dim(as.matrix(written[, -1])), dim(expected_z))
    expect_equal(as.vector(as.matrix(written[, -1])), as.vector(expected_z), tolerance = 1e-9)
    expect_false(any(grepl("clean_values|residual_values|pca", list.files(adjusted_dir))))
  })
}

test_that("prepare CLI reads gzipped NPX and preserves sample IDs", {
  gzip_input <- file.path(work, "input's compressed data.tsv.gz")
  gzip_connection <- gzfile(gzip_input, "wt")
  writeLines(readLines(input), gzip_connection)
  close(gzip_connection)
  directory <- file.path(work, "prepared gzip NPX")
  result <- run_cli("prepare", "--input", gzip_input, "--metric", "NPX",
                    "--max-pcs", "0", "--excluded-assays-file", exclusions,
                    "--threads", "1", "--output-dir", directory)
  expect_equal(result$status, 0L, info = result$text)
  state <- readRDS(file.path(directory, "prepared.rds"))
  expect_identical(state$sample_qc$SampleID, prepared$NPX$sample_qc$SampleID)
  expect_true(all(c("002", "003's") %in% rownames(state$scaled_values)))
  expect_identical(state$scaled_values, prepared$NPX$scaled_values)
})

test_that("configured exclusions precede duplicate and sample ID checks", {
  baseline <- prepare_analysis(fixture, "NPX", 4L, c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE)
  repeated <- fixture |>
    filter(Assay %in% c("GBP1", "MAP2K1")) |>
    mutate(NPX = NPX + 100, PCNormalizedNPX = PCNormalizedNPX - 100)
  data <- bind_rows(fixture, repeated)
  state <- NULL
  expect_no_error(state <- prepare_analysis(data, "NPX", 4L, c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE))
  if (!is.null(state)) {
    expect_identical(state$scaled_values, baseline$scaled_values)
    expect_identical(state$sample_qc, baseline$sample_qc)
  }
  expect_error(prepare_analysis(data, "NPX", 0L, character(), -3, 7L, 1L, TRUE),
               "Duplicate assay/sample")
  excluded_ids <- tibble(Assay = c("GBP1", "MAP2K1"),
                         SampleID = c("excluded-only_visit1", "004_visit2"),
                         NPX = c(100, 99), PCNormalizedNPX = c(99, 100))
  state <- NULL
  expect_no_error(state <- prepare_analysis(bind_rows(fixture, excluded_ids), "NPX", 4L,
                                            c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE))
  if (!is.null(state)) {
    expect_identical(state$sample_qc, baseline$sample_qc)
    expect_identical(state$prepare_summary$Samples_input, baseline$prepare_summary$Samples_input)
  }
})

test_that("sample suffix rules reject collisions and duplicate assay/sample keys", {
  duplicated <- bind_rows(fixture, fixture[1, ])
  expect_error(prepare_analysis(duplicated, "NPX", 0L, c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE),
               "Duplicate assay/sample")
  collided <- fixture
  collided$SampleID[collided$Assay == "Protein01" & collided$SampleID == "004_visit1"] <- "004_visit2"
  expect_error(prepare_analysis(collided, "NPX", 0L, c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE),
               "merge distinct input IDs")
  result <- prepare_analysis(fixture, "NPX", 0L, c("GBP1", "MAP2K1"), -3, 7L, 1L, FALSE)
  expect_true("003's_visit1" %in% rownames(result$scaled_values))
})

test_that("sample QC removes new constant assays and rescales the retained samples", {
  spike <- tibble(Assay = "AfterQCConstant", SampleID = unique(fixture$SampleID),
                  NPX = as.numeric(sub("_.*$", "", unique(fixture$SampleID)) == "024"),
                  PCNormalizedNPX = 0)
  data <- bind_rows(fixture, spike)
  state <- prepare_analysis(data, "NPX", 2L, c("GBP1", "MAP2K1"), -1, 7L, 1L, TRUE)
  expect_true("024" %in% state$sample_qc$SampleID[state$sample_qc$Status == "connectivity_outlier"])
  expect_equal(state$excluded_assays$Stage[state$excluded_assays$Assay == "AfterQCConstant"],
               "after_connectivity")
  expect_equal(unname(colMeans(state$scaled_values)), rep(0, ncol(state$scaled_values)), tolerance = 1e-12)
  expect_equal(unname(apply(state$scaled_values, 2, sd)), rep(1, ncol(state$scaled_values)))
  reference <- reference_result(data, "NPX", 2, cutoff = -1)
  expect_equal(as.vector(adjust_analysis(state, 2L)$z_scores), as.vector(reference$z), tolerance = 1e-9)
})

test_that("prepare rejects impossible PC counts and insufficient matrix rank", {
  expect_error(prepare_analysis(fixture, "NPX", 1000L, c("GBP1", "MAP2K1"), -3, 7L, 1L, TRUE),
               "at most.*residual")
  low_rank <- fixture
  i <- match(low_rank$SampleID, unique(low_rank$SampleID))
  j <- match(low_rank$Assay, unique(low_rank$Assay))
  low_rank$NPX <- sin(i / 3) * cos(j / 2) + cos(i / 5) * sin(j / 3)
  expect_error(prepare_analysis(low_rank, "NPX", 2L, c("GBP1", "MAP2K1"), -100, 7L, 1L, TRUE),
               "effective data rank is 2")
  expect_error(adjust_analysis(prepared$NPX, 5L, 3, "NPX"), "prepared|maximum")
  expect_error(adjust_analysis(prepared$NPX, 2L, 3, "PCNormalizedNPX"), "metric.*match|does not match")
})

test_that("residual Z uses population SD and strict protein cutoffs", {
  state <- prepared$NPX
  adjusted <- adjust_analysis(state, 0L, 3, "NPX")
  expect_equal(unname(colMeans(adjusted$z_scores^2)), rep(1, ncol(adjusted$z_scores)))
  edge <- max(abs(adjusted$z_scores))
  expect_equal(nrow(adjust_analysis(state, 0L, edge, "NPX")$protein_outliers), 0L)
  x <- c(-1, 1, -1, 1)
  y <- c(-1, -1, 1, 1)
  state$scaled_values <- scale(cbind(A = x, B = x, C = y))
  rownames(state$scaled_values) <- letters[1:4]
  pca <- prcomp(state$scaled_values, center = FALSE, scale. = FALSE, rank. = 1)
  state$scores <- pca$x
  state$loadings <- pca$rotation
  undefined <- adjust_analysis(state, 1L, 3, "NPX")
  expect_true(all(is.na(undefined$z_scores[, c("A", "B")])))
  expect_setequal(undefined$undefined_residual_assays$Assay, c("A", "B"))
})

test_that("input checks report cloud localization errors and missing files before opening", {
  for (path in c("gs://bucket/input.tsv", "s3://bucket/input.tsv", "https://example.org/input.tsv")) {
    expect_error(assert_readable(path), "localization error")
  }
  expect_error(assert_readable(file.path(work, "missing's input.tsv")), "not readable")
  result <- run_cli("prepare", "--input", "gs://bucket/input.tsv", "--metric", "NPX",
                    "--max-pcs", "0", "--output-dir", file.path(work, "bad cloud"))
  expect_true(result$status != 0L)
  expect_match(result$text, "localization error")
})

test_that("empty exclusions files replace the default assay exclusion list", {
  empty <- file.path(work, "empty exclusions.txt")
  file.create(empty)
  expect_identical(read_exclusions(empty), character())
  expect_equal(read_exclusions(NULL), c("GBP1", "MAP2K1"))
  state <- prepare_analysis(fixture, "NPX", 0L, character(), -3, 7L, 1L, TRUE)
  expect_true(all(c("GBP1", "MAP2K1") %in% colnames(state$scaled_values)))
})

test_that("collect reads localized summary files and keeps filenames as metadata", {
  summaries <- file.path(work, paste0("adjusted ", c("NPX", "PCNormalizedNPX")),
                         paste0("summary.", c("NPX", "PCNormalizedNPX"), ".PCs2.tsv"))
  file_list <- file.path(work, "summary paths's list.txt")
  writeLines(summaries, file_list)
  output <- file.path(work, "combined summary.tsv")
  result <- run_cli("collect", "--summary-files", file_list, "--output", output)
  expect_equal(result$status, 0L, info = result$text)
  combined <- readr::read_tsv(output, show_col_types = FALSE)
  expect_equal(combined$NPX_type, c("NPX", "PCNormalizedNPX"))
  expect_equal(combined$PCs_removed, c(2, 2))
  file_columns <- grep("_file$", names(combined), value = TRUE)
  expect_true(length(file_columns) > 0)
  expect_false(any(grepl("/", unlist(combined[file_columns]))))
  writeLines("gs://bucket/unlocalized.tsv", file_list)
  expect_error(collect_summaries(file_list, output), "localization error")
})

unlink(work, recursive = TRUE)
cat("Olink smoke tests passed.\n")
