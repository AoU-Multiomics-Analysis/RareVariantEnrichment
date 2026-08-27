suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

argument <- function(name, default = NULL) {
  index <- match(name, args)
  if (is.na(index)) {
    if (!is.null(default)) {
      return(default)
    }
    stop("Missing command-line argument: ", name)
  }
  if (index == length(args)) {
    stop("Missing value for command-line argument: ", name)
  }
  args[[index + 1L]]
}

parse_csv <- function(value, name) {
  values <- str_split(value, ",", simplify = TRUE) |>
    as.character() |>
    str_trim()
  if (length(values) == 0L || any(values == "") || anyDuplicated(values)) {
    stop(name, " must contain unique non-empty values")
  }
  values
}

results_input <- argument("--results-input")
summary_output <- argument("--summary-output")
plot_output <- argument("--plot-output")
thresholds <- parse_csv(
  argument("--selection-z-thresholds"),
  "--selection-z-thresholds"
) |>
  as.numeric()
definitions <- parse_csv(
  argument("--carrier-definitions", "HC,any_lof"),
  "--carrier-definitions"
)

if (length(thresholds) == 0L || any(!is.finite(thresholds)) ||
    any(thresholds >= 0) || anyDuplicated(thresholds)) {
  stop("--selection-z-thresholds must contain unique finite negative values")
}

required_columns <- c("pc_count", "z_threshold", "carrier_definition", "odds_ratio")
results <- readr::read_tsv(
  results_input,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("NA", "NaN", "")
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0L) {
  stop(
    "Results TSV is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

selected_results <- results |>
  transmute(
    pc_count = as.integer(pc_count),
    z_threshold = as.numeric(z_threshold),
    carrier_definition = as.character(carrier_definition),
    odds_ratio = as.numeric(odds_ratio)
  ) |>
  filter(
    z_threshold %in% thresholds,
    carrier_definition %in% definitions
  ) |>
  mutate(
    finite_odds_ratio = is.finite(odds_ratio) & odds_ratio > 0,
    carrier_definition = factor(carrier_definition, levels = definitions)
  )

if (nrow(selected_results) == 0L) {
  stop("No enrichment rows matched the requested thresholds and definitions")
}

summary_table <- selected_results |>
  group_by(carrier_definition, pc_count) |>
  summarise(
    n_thresholds = n_distinct(z_threshold[finite_odds_ratio]),
    thresholds_present = paste(
      sort(unique(z_threshold[finite_odds_ratio]), decreasing = TRUE),
      collapse = ","
    ),
    median_log_odds_ratio = if_else(
      any(finite_odds_ratio),
      median(log(odds_ratio[finite_odds_ratio])),
      NA_real_
    ),
    median_odds_ratio = if_else(
      any(finite_odds_ratio),
      exp(median(log(odds_ratio[finite_odds_ratio]))),
      NA_real_
    ),
    .groups = "drop"
  )

maximums <- summary_table |>
  filter(is.finite(median_odds_ratio)) |>
  arrange(carrier_definition, desc(median_odds_ratio), pc_count) |>
  group_by(carrier_definition) |>
  slice_head(n = 1L) |>
  transmute(
    carrier_definition,
    maximum_pc = pc_count,
    maximum_median_odds_ratio = median_odds_ratio
  ) |>
  ungroup()

summary_table <- summary_table |>
  left_join(maximums, by = "carrier_definition") |>
  mutate(
    percent_of_maximum = 100 * median_odds_ratio / maximum_median_odds_ratio,
    carrier_definition = factor(carrier_definition, levels = definitions)
  ) |>
  arrange(carrier_definition, pc_count)

readr::write_tsv(summary_table, summary_output, na = "NA")

plot_table <- summary_table |>
  filter(is.na(maximum_pc) | pc_count <= maximum_pc)
label_pcs <- unique(c(
  0L, 100L, 500L, 1000L, 1500L, 2000L,
  maximums$maximum_pc
))
plot_labels <- plot_table |>
  filter(
    pc_count %in% label_pcs,
    is.finite(percent_of_maximum),
    percent_of_maximum < 100
  ) |>
  mutate(label = sprintf("OR = %.1f", median_odds_ratio))

plot_maximums <- maximums |>
  mutate(
    carrier_definition = factor(carrier_definition, levels = definitions),
    label = paste0(
      "maximum PC = ", maximum_pc,
      "\nOR = ",
      formatC(maximum_median_odds_ratio, format = "f", digits = 1),
      " (100%)"
    ),
    label_y = 106
  )

plot <- ggplot(
  plot_table,
  aes(x = pc_count, y = percent_of_maximum, group = carrier_definition)
) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey45") +
  geom_line(color = "#2563EB", linewidth = 0.8, na.rm = TRUE) +
  geom_point(color = "#1D4ED8", size = 1.7, na.rm = TRUE) +
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
  facet_wrap(
    ~carrier_definition,
    nrow = 1,
    scales = "free_x",
    drop = FALSE
  ) +
  scale_y_continuous(
    name = "Median enrichment as percent of maximum",
    limits = c(0, 108),
    breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  scale_x_continuous(name = "Number of PCs included (k)") +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey75"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )

ggsave(plot_output, plot, width = 11, height = 6, dpi = 300)
