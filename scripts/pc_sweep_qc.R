suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

argument <- function(name) {
  index <- match(name, args)
  if (is.na(index) || index == length(args)) {
    stop("Missing command-line argument: ", name)
  }
  args[[index + 1L]]
}

results_input <- argument("--results-input")
summary_output <- argument("--summary-output")
plot_output <- argument("--plot-output")
thresholds <- as.numeric(strsplit(argument("--selection-z-thresholds"), ",", fixed = TRUE)[[1L]])

if (length(thresholds) == 0L || any(!is.finite(thresholds)) ||
    any(thresholds >= 0) || anyDuplicated(thresholds)) {
  stop("--selection-z-thresholds must contain unique finite negative values")
}

definitions <- c("HC", "any_lof")
required_columns <- c("pc_count", "z_threshold", "carrier_definition", "odds_ratio")
results <- fread(results_input)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0L) {
  stop(
    "Results TSV is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

results[, c("pc_count", "z_threshold", "odds_ratio") := list(
  as.integer(pc_count),
  as.numeric(z_threshold),
  as.numeric(odds_ratio)
)]
filtered <- results[
  z_threshold %in% thresholds &
    carrier_definition %in% definitions &
    is.finite(odds_ratio) & odds_ratio > 0
]
if (nrow(filtered) == 0L) {
  stop("No finite odds-ratio rows matched the requested thresholds and definitions")
}

summary_table <- filtered[
  , .(
    n_thresholds = uniqueN(z_threshold),
    thresholds_present = paste(sort(unique(z_threshold), decreasing = TRUE), collapse = ","),
    median_log_odds_ratio = median(log(odds_ratio)),
    median_odds_ratio = exp(median(log(odds_ratio)))
  ),
  by = .(carrier_definition, pc_count)
]
incomplete <- summary_table[n_thresholds != length(thresholds)]
if (nrow(incomplete) > 0L) {
  stop(
    "Some PC/definition combinations do not contain all requested thresholds: ",
    paste0(incomplete$carrier_definition, ":PC", incomplete$pc_count, collapse = ", ")
  )
}

maximums <- summary_table[
  order(carrier_definition, -median_odds_ratio, pc_count),
  .(
    maximum_pc = pc_count[1L],
    maximum_median_odds_ratio = median_odds_ratio[1L]
  ),
  by = carrier_definition
]
summary_table <- merge(summary_table, maximums, by = "carrier_definition", sort = FALSE)
summary_table[, percent_of_maximum := 100 * median_odds_ratio / maximum_median_odds_ratio]
summary_table[, carrier_definition := factor(carrier_definition, levels = definitions)]
setorder(summary_table, carrier_definition, pc_count)
fwrite(summary_table, summary_output, sep = "\t", quote = FALSE, na = "NA")

plot_table <- summary_table[pc_count <= maximum_pc]
label_pcs <- unique(c(0L, 100L, 500L, 1000L, 1500L, 2000L, maximums$maximum_pc))
plot_labels <- plot_table[pc_count %in% label_pcs & percent_of_maximum < 100]
plot_labels[, label := sprintf("OR = %.1f", median_odds_ratio)]

plot_maximums <- unique(maximums[, .(
  carrier_definition,
  maximum_pc,
  maximum_median_odds_ratio
)])
plot_maximums[, carrier_definition := factor(carrier_definition, levels = definitions)]
plot_maximums[, label := paste0(
  "maximum PC = ", maximum_pc,
  "\nOR = ", formatC(maximum_median_odds_ratio, format = "f", digits = 1),
  " (100%)"
)]
plot_maximums[, label_y := 106]

plot <- ggplot(plot_table, aes(x = pc_count, y = percent_of_maximum)) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey45") +
  geom_line(color = "#2563EB", linewidth = 0.8) +
  geom_point(color = "#1D4ED8", size = 1.7) +
  geom_text_repel(
    data = plot_labels,
    aes(label = label),
    direction = "both",
    force = 1.5,
    force_pull = 0.5,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 123,
    size = 2.7
  ) +
  geom_vline(
    data = plot_maximums,
    aes(xintercept = maximum_pc),
    linetype = "dotted",
    color = "#111827"
  ) +
  geom_label_repel(
    data = plot_maximums,
    aes(x = maximum_pc, y = label_y, label = label),
    inherit.aes = FALSE,
    direction = "both",
    force = 1.5,
    force_pull = 0.5,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 123,
    size = 3,
    label.size = 0.2
  ) +
  facet_wrap(~carrier_definition, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    name = "Median enrichment as percent of maximum",
    limits = c(0, 108),
    breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  scale_x_continuous(name = "Number of PCs included (k)") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "grey65"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

ggsave(plot_output, plot, width = 11, height = 6, dpi = 300)
