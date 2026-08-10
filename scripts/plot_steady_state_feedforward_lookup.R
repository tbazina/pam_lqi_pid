#!/usr/bin/env Rscript

# Plot the active branchwise steady-state feedforward lookup in a compact,
# manuscript-ready two-panel layout.

script_args <- commandArgs(trailingOnly = FALSE)
script_match <- script_args[grepl("^--file=", script_args)]
script_dir <- if (length(script_match) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_match[1]), mustWork = TRUE))
} else {
  file.path(getwd(), "scripts")
}

source(file.path(script_dir, "lib", "setup.R"))
load_analysis_packages()
context <- initialize_analysis_context()

required_columns <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")

read_lookup <- function(path, expected_branch) {
  if (!file.exists(path)) {
    stop(sprintf("Missing lookup table: %s", path))
  }

  lookup <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(lookup) <- trimws(names(lookup))
  if (!identical(names(lookup), required_columns)) {
    stop(sprintf(
      "Unexpected columns in %s. Expected: %s",
      path,
      paste(required_columns, collapse = ", ")
    ))
  }

  if (nrow(lookup) < 2L || any(lookup$branch != expected_branch)) {
    stop(sprintf("Invalid or incomplete %s lookup: %s", expected_branch, path))
  }

  numeric_columns <- lookup[c("s_ref_mm", "u_ff_v", "p_ff_bar")]
  if (!all(vapply(numeric_columns, function(x) all(is.finite(x)), logical(1)))) {
    stop(sprintf("Non-finite values found in %s", path))
  }
  if (anyDuplicated(lookup$s_ref_mm)) {
    stop(sprintf("Duplicate reference values found in %s", path))
  }

  lookup %>%
    dplyr::arrange(s_ref_mm) %>%
    dplyr::mutate(branch = expected_branch)
}

up_lookup <- read_lookup(
  file.path(context$output_dir, "lqr_steady_state_lookup_up.csv"),
  "up"
)
down_lookup <- read_lookup(
  file.path(context$output_dir, "lqr_steady_state_lookup_down.csv"),
  "down"
)

lookup_long <- dplyr::bind_rows(up_lookup, down_lookup) %>%
  tidyr::pivot_longer(
    cols = c(u_ff_v, p_ff_bar),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    branch_label = factor(
      dplyr::recode(branch, up = "Upward branch", down = "Downward branch"),
      levels = c("Upward branch", "Downward branch")
    ),
    metric_label = factor(
      dplyr::recode(
        metric,
        u_ff_v = "Feedforward voltage [V]",
        p_ff_bar = "Feedforward pressure [bar]"
      ),
      levels = c("Feedforward voltage [V]", "Feedforward pressure [bar]")
    )
  )

plot_theme <- create_analysis_plot_theme() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
    panel.spacing.x = grid::unit(3, "mm"),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_text(
      face = "plain",
      size = 10,
      angle = 0,
      hjust = 0.5,
      margin = ggplot2::margin(0, 0, 1, 0, "mm")
    ),
    plot.margin = ggplot2::margin(1.5, 2.5, 1.5, 2.5, "mm")
  )

lookup_plot <- ggplot2::ggplot(
  lookup_long,
  ggplot2::aes(x = value, y = s_ref_mm, color = branch_label, group = branch_label)
) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::facet_grid(
    . ~ metric_label, scales = "free_x"
  ) +
  ggplot2::scale_color_manual(
    name = "Branch",
    values = c("Upward branch" = "#D55E00", "Downward branch" = "#0072B2")
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Reference displacement [mm]"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      override.aes = list(linewidth = 0.45)
    )
  ) +
  plot_theme

output_path <- file.path(context$output_dir, "steady_state_feedforward_lookup.png")
ggplot2::ggsave(
  output_path,
  lookup_plot,
  width = 15,
  height = 5.1724,
  units = "cm",
  dpi = 600,
  bg = "white"
)

cat(sprintf("Saved steady-state feedforward lookup plot to: %s\n", output_path))
cat(sprintf("Upward lookup points: %d; downward lookup points: %d\n", nrow(up_lookup), nrow(down_lookup)))
