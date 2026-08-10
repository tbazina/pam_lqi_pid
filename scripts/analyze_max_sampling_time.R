#!/usr/bin/env Rscript

# Resolve the entry script directory for both `Rscript` and sourced execution.
find_entry_script_dir <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) {
    script_path <- sub("^--file=", "", script_match[1])
    return(dirname(normalizePath(script_path)))
  }

  current_candidate <- file.path(getwd(), "scripts")
  if (dir.exists(current_candidate)) {
    return(normalizePath(current_candidate))
  }

  normalizePath(getwd())
}

# Source one helper file from the local analysis library.
source_analysis_helper <- function(script_dir, helper_file) {
  source(file.path(script_dir, "lib", helper_file), local = FALSE)
}

# Read the loop-cycle timing experiment and enforce the expected schema.
load_sampling_time_data <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing sampling-time experiment file: %s", path))
  }

  raw_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  required_cols <- "t[ms]"
  missing_cols <- setdiff(required_cols, names(raw_df))

  if (length(missing_cols) > 0) {
    stop("Sampling-time experiment file must contain column: t[ms]")
  }

  raw_df %>%
    transmute(
      sequence_id = row_number(),
      t_ms = as.numeric(`t[ms]`)
    )
}

# Flag rows that should contribute to the summary and identify the worst-case sample.
annotate_sampling_time_points <- function(df) {
  annotated_df <- df %>%
    mutate(
      finite_row = is.finite(t_ms),
      positive_row = finite_row & t_ms > 0,
      keep_for_summary = positive_row
    )

  if (!any(annotated_df$keep_for_summary)) {
    stop("No finite positive sampling-time rows remain for analysis.")
  }

  max_sequence_id <- annotated_df %>%
    filter(keep_for_summary) %>%
    arrange(desc(t_ms), sequence_id) %>%
    slice(1) %>%
    pull(sequence_id)

  annotated_df %>%
    mutate(
      is_max_row = keep_for_summary & sequence_id == max_sequence_id,
      exclude_reason = case_when(
        !finite_row ~ "non_finite",
        finite_row & t_ms <= 0 ~ "non_positive",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(sequence_id, t_ms, keep_for_summary, is_max_row, exclude_reason)
}

# Build the global timing summary used later by excitation design.
create_sampling_time_summary <- function(points_df) {
  used_df <- points_df %>%
    filter(keep_for_summary)

  max_row <- used_df %>%
    arrange(desc(t_ms), sequence_id) %>%
    slice(1)

  tibble(
    level = "global",
    source_file = "experiment/max_sampling_time.csv",
    n_raw = nrow(points_df),
    n_used = nrow(used_df),
    n_non_finite = sum(points_df$exclude_reason == "non_finite", na.rm = TRUE),
    n_non_positive = sum(points_df$exclude_reason == "non_positive", na.rm = TRUE),
    min_sampling_time_ms = min(used_df$t_ms, na.rm = TRUE),
    mean_sampling_time_ms = mean(used_df$t_ms, na.rm = TRUE),
    median_sampling_time_ms = median(used_df$t_ms, na.rm = TRUE),
    p90_sampling_time_ms = as.numeric(stats::quantile(
      used_df$t_ms,
      probs = 0.90,
      na.rm = TRUE,
      names = FALSE
    )),
    p95_sampling_time_ms = as.numeric(stats::quantile(
      used_df$t_ms,
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    )),
    p99_sampling_time_ms = as.numeric(stats::quantile(
      used_df$t_ms,
      probs = 0.99,
      na.rm = TRUE,
      names = FALSE
    )),
    observed_max_sampling_time_ms = max_row$t_ms[[1]],
    max_sampling_time_sequence_id = max_row$sequence_id[[1]],
    recommended_max_sampling_time_ms = max_row$t_ms[[1]],
    recommended_sampling_period_s = max_row$t_ms[[1]] / 1000,
    max_supported_sampling_rate_hz = 1000 / max_row$t_ms[[1]],
    mean_supported_sampling_rate_hz = 1000 / mean(used_df$t_ms, na.rm = TRUE),
    recommendation_basis = "observed_max_loop_cycle_time"
  )
}

# Save the point-level timing annotations for later inspection.
write_sampling_time_points <- function(path, points_df) {
  points_df %>%
    write_csv(path)
}

# Save the global timing summary.
write_sampling_time_summary <- function(path, summary_df) {
  summary_df %>%
    write_csv(path)
}

# Save a diagnostic plot of loop-cycle timing over the experiment.
save_sampling_time_plot <- function(path, points_df, summary_df, analysis_plot_theme) {
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)

  plot_obj <- ggplot(points_df, aes(x = sequence_id, y = t_ms)) +
    geom_hline(
      yintercept = global_row$mean_sampling_time_ms[[1]],
      color = "#1F1F1F",
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_hline(
      yintercept = global_row$p95_sampling_time_ms[[1]],
      color = "#7A7A7A",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_line(color = "#4C72B0", linewidth = 0.35, alpha = 0.8) +
    geom_point(
      data = points_df %>% filter(keep_for_summary),
      color = "#4C72B0",
      size = 0.8,
      alpha = 0.8
    ) +
    geom_point(
      data = points_df %>% filter(is_max_row),
      color = "#C44E52",
      size = 1.5
    ) +
    labs(
      x = "Loop-cycle index",
      y = "Sampling time [ms]"
    ) +
    analysis_plot_theme

  ggsave(
    filename = path,
    plot = plot_obj,
    width = 180,
    height = 110,
    units = "mm",
    dpi = 300
  )
}

# Print a concise summary for terminal use.
print_sampling_time_summary <- function(summary_df, output_dir) {
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)

  cat("\nMaximum sampling time analysis\n")
  cat(sprintf(
    "Observed max sampling time: %.6f ms at sequence %d\n",
    global_row$observed_max_sampling_time_ms[[1]],
    global_row$max_sampling_time_sequence_id[[1]]
  ))
  cat(sprintf(
    "Recommended sampling period: %.9f s\n",
    global_row$recommended_sampling_period_s[[1]]
  ))
  cat(sprintf(
    "Fastest guaranteed loop rate from observed max: %.3f Hz\n",
    global_row$max_supported_sampling_rate_hz[[1]]
  ))
  cat(sprintf(
    "Distribution summary: mean %.6f ms, p95 %.6f ms, p99 %.6f ms\n",
    global_row$mean_sampling_time_ms[[1]],
    global_row$p95_sampling_time_ms[[1]],
    global_row$p99_sampling_time_ms[[1]]
  ))
  cat(sprintf("Outputs written to: %s\n", output_dir))
}

# Load helper modules in dependency order, then initialize shared context.
script_dir <- find_entry_script_dir()

source_analysis_helper(script_dir, "setup.R")
load_analysis_packages()

analysis_context <- initialize_analysis_context()
analysis_plot_theme <- create_analysis_plot_theme()

repo_dir <- analysis_context$repo_dir
output_dir <- analysis_context$output_dir

input_path <- file.path(repo_dir, "experiment", "max_sampling_time.csv")
summary_path <- file.path(output_dir, "max_sampling_time_summary.csv")
points_path <- file.path(output_dir, "max_sampling_time_points.csv")
plot_path <- file.path(output_dir, "max_sampling_time_diagnostics.png")

points_df <- load_sampling_time_data(input_path) %>%
  annotate_sampling_time_points()
summary_df <- create_sampling_time_summary(points_df)

write_sampling_time_points(points_path, points_df)
write_sampling_time_summary(summary_path, summary_df)
save_sampling_time_plot(plot_path, points_df, summary_df, analysis_plot_theme)
print_sampling_time_summary(summary_df, output_dir)
