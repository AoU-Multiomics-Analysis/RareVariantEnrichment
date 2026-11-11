#!/usr/bin/env Rscript
# Synthetic data only. This fixture has different missing samples for each metric.
make_olink_fixture <- function() {
  set.seed(2701)
  samples <- sprintf("%03d_visit1", seq_len(48L))
  samples[3] <- "003's_visit1"
  assays <- c(sprintf("Protein%02d", seq_len(12L)), "Constant", "GBP1", "MAP2K1")
  factors <- matrix(rnorm(length(samples) * 3L), nrow = length(samples))
  values <- factors %*% matrix(rnorm(3L * length(assays)), nrow = 3L) +
    matrix(rnorm(length(samples) * length(assays)), nrow = length(samples))
  normalized <- values * 0.7 + matrix(rnorm(length(values)), nrow = nrow(values))
  values[, match("Constant", assays)] <- 5
  normalized[, match("Constant", assays)] <- 8
  values[1, 1] <- NA_real_
  normalized[2, 2] <- Inf
  data.frame(
    Assay = rep(assays, each = length(samples)),
    SampleID = rep(samples, times = length(assays)),
    NPX = as.vector(values), PCNormalizedNPX = as.vector(normalized),
    check.names = FALSE
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply one output TSV path.", call. = FALSE)
  data.table::fwrite(make_olink_fixture(), args[1], sep = "\t", na = "NA")
}
