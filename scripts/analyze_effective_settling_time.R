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

# Read the dedicated settling-time experiment and enforce the expected schema.
load_settling_time_data <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing settling-time experiment file: %s", path))
  }

  raw_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  required_cols <- c("u[V]", "t[ms]")
  missing_cols <- setdiff(required_cols, names(raw_df))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Settling-time experiment file must contain columns: %s",
      paste(required_cols, collapse = ", ")
    ))
  }

  raw_df %>%
    transmute(
      sequence_id = row_number(),
      u = as.numeric(`u[V]`),
      t_ms = as.numeric(`t[ms]`)
    )
}

# Read the global control window if it exists so the settling test can be annotated.
read_control_window_annotation <- function(output_dir, margin_fraction = 0.05) {
  path <- file.path(output_dir, "recommended_control_window.csv")

  if (!file.exists(path)) {
    return(list(
      available = FALSE,
      margin_fraction = margin_fraction,
      window_low_v = NA_real_,
      window_high_v = NA_real_,
      excitation_low_v = NA_real_,
      excitation_high_v = NA_real_
    ))
  }

  window_df <- read_csv(path, show_col_types = FALSE)
  global_row <- window_df %>%
    filter(level == "global") %>%
    slice(1)

  if (nrow(global_row) == 0) {
    return(list(
      available = FALSE,
      margin_fraction = margin_fraction,
      window_low_v = NA_real_,
      window_high_v = NA_real_,
      excitation_low_v = NA_real_,
      excitation_high_v = NA_real_
    ))
  }

  window_low_v <- as.numeric(global_row$window_low_v[[1]])
  window_high_v <- as.numeric(global_row$window_high_v[[1]])
  margin_u <- margin_fraction * (window_high_v - window_low_v)

  list(
    available = is.finite(window_low_v) && is.finite(window_high_v),
    margin_fraction = margin_fraction,
    window_low_v = window_low_v,
    window_high_v = window_high_v,
    excitation_low_v = window_low_v + margin_u,
    excitation_high_v = window_high_v - margin_u
  )
}

# Describe whether a target lies below, inside, or above a given window.
classify_against_window <- function(u_value, low_v, high_v) {
  if (!is.finite(u_value) || !is.finite(low_v) || !is.finite(high_v)) {
    return(NA_character_)
  }

  if (u_value < low_v) {
    return("below")
  }

  if (u_value > high_v) {
    return("above")
  }

  "inside"
}

# Flag startup rows and robust target-level outliers.
flag_settling_time_points <- function(df) {
  base_df <- df %>%
    mutate(
      finite_row = is.finite(u) & is.finite(t_ms),
      analysis_candidate = finite_row & sequence_id > 1 & u > 0,
      outlier_flag = FALSE
    )

  design_df <- base_df %>%
    filter(analysis_candidate) %>%
    group_by(u) %>%
    mutate(
      median_t_ms = median(t_ms, na.rm = TRUE),
      mad_t_ms = stats::mad(t_ms, constant = 1, na.rm = TRUE),
      outlier_flag = if_else(
        is.finite(mad_t_ms) & mad_t_ms > 0,
        abs(t_ms - median_t_ms) > 5 * mad_t_ms,
        FALSE
      )
    ) %>%
    ungroup() %>%
    dplyr::select(sequence_id, outlier_flag)

  base_df %>%
    left_join(design_df, by = "sequence_id", suffix = c("", "_design")) %>%
    mutate(
      outlier_flag = coalesce(outlier_flag_design, outlier_flag),
      keep_for_recommendation = analysis_candidate & !outlier_flag,
      exclude_reason = case_when(
        !finite_row ~ "non_finite",
        sequence_id == 1 ~ "initial_warmup_row",
        u <= 0 ~ "startup_or_non_design",
        outlier_flag ~ "robust_outlier",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(sequence_id, u, t_ms, keep_for_recommendation, outlier_flag, exclude_reason)
}

# Build per-target and global settling-time summaries.
create_effective_settling_time_summary <- function(points_df, control_window_info) {
  target_points <- points_df %>%
    filter(is.finite(u), is.finite(t_ms), u > 0)

  if (nrow(target_points) == 0) {
    stop("No positive, finite settling-time rows remain after startup filtering.")
  }

  target_summary <- target_points %>%
    group_by(u) %>%
    summarise(
      n_raw = n(),
      n_used = sum(keep_for_recommendation),
      n_outliers = sum(outlier_flag),
      raw_median_t_ms = median(t_ms, na.rm = TRUE),
      raw_p90_t_ms = as.numeric(stats::quantile(t_ms, probs = 0.90, na.rm = TRUE, names = FALSE)),
      raw_p95_t_ms = as.numeric(stats::quantile(t_ms, probs = 0.95, na.rm = TRUE, names = FALSE)),
      raw_max_t_ms = max(t_ms, na.rm = TRUE),
      trimmed_median_t_ms = median(t_ms[keep_for_recommendation], na.rm = TRUE),
      trimmed_p90_t_ms = as.numeric(stats::quantile(
        t_ms[keep_for_recommendation],
        probs = 0.90,
        na.rm = TRUE,
        names = FALSE
      )),
      trimmed_p95_t_ms = as.numeric(stats::quantile(
        t_ms[keep_for_recommendation],
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )),
      trimmed_max_t_ms = max(t_ms[keep_for_recommendation], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      level = "target",
      target_u_v = u,
      recommended_t_settle_ms = trimmed_p95_t_ms,
      recommended_t_hold_s = recommended_t_settle_ms / 1000,
      recommendation_basis = "per_target_p95_trimmed",
      control_window_available = control_window_info$available,
      control_window_low_v = control_window_info$window_low_v,
      control_window_high_v = control_window_info$window_high_v,
      excitation_margin_fraction = control_window_info$margin_fraction,
      excitation_low_v = control_window_info$excitation_low_v,
      excitation_high_v = control_window_info$excitation_high_v,
      relation_to_control_window = purrr::map_chr(
        target_u_v,
        ~classify_against_window(
          .x,
          control_window_info$window_low_v,
          control_window_info$window_high_v
        )
      ),
      relation_to_excitation_window = purrr::map_chr(
        target_u_v,
        ~classify_against_window(
          .x,
          control_window_info$excitation_low_v,
          control_window_info$excitation_high_v
        )
      ),
      inside_control_window = relation_to_control_window == "inside",
      inside_excitation_window = relation_to_excitation_window == "inside",
      driver_target_u_v = NA_real_
    ) %>%
    dplyr::select(
      level,
      target_u_v,
      n_raw,
      n_used,
      n_outliers,
      raw_median_t_ms,
      raw_p90_t_ms,
      raw_p95_t_ms,
      raw_max_t_ms,
      trimmed_median_t_ms,
      trimmed_p90_t_ms,
      trimmed_p95_t_ms,
      trimmed_max_t_ms,
      recommended_t_settle_ms,
      recommended_t_hold_s,
      recommendation_basis,
      control_window_available,
      control_window_low_v,
      control_window_high_v,
      excitation_margin_fraction,
      excitation_low_v,
      excitation_high_v,
      relation_to_control_window,
      relation_to_excitation_window,
      inside_control_window,
      inside_excitation_window,
      driver_target_u_v
    )

  if (any(!is.finite(target_summary$recommended_t_settle_ms))) {
    stop("At least one target level has no retained rows left for settling-time recommendation.")
  }

  driver_row <- target_summary %>%
    arrange(desc(recommended_t_settle_ms), target_u_v) %>%
    slice(1)

  bind_rows(
    target_summary,
    tibble(
      level = "global",
      target_u_v = NA_real_,
      n_raw = sum(target_summary$n_raw),
      n_used = sum(target_summary$n_used),
      n_outliers = sum(target_summary$n_outliers),
      raw_median_t_ms = NA_real_,
      raw_p90_t_ms = NA_real_,
      raw_p95_t_ms = NA_real_,
      raw_max_t_ms = max(target_summary$raw_max_t_ms, na.rm = TRUE),
      trimmed_median_t_ms = NA_real_,
      trimmed_p90_t_ms = NA_real_,
      trimmed_p95_t_ms = NA_real_,
      trimmed_max_t_ms = max(target_summary$trimmed_max_t_ms, na.rm = TRUE),
      recommended_t_settle_ms = driver_row$recommended_t_settle_ms[[1]],
      recommended_t_hold_s = driver_row$recommended_t_hold_s[[1]],
      recommendation_basis = "max_target_p95_trimmed",
      control_window_available = control_window_info$available,
      control_window_low_v = control_window_info$window_low_v,
      control_window_high_v = control_window_info$window_high_v,
      excitation_margin_fraction = control_window_info$margin_fraction,
      excitation_low_v = control_window_info$excitation_low_v,
      excitation_high_v = control_window_info$excitation_high_v,
      relation_to_control_window = driver_row$relation_to_control_window[[1]],
      relation_to_excitation_window = driver_row$relation_to_excitation_window[[1]],
      inside_control_window = driver_row$inside_control_window[[1]],
      inside_excitation_window = driver_row$inside_excitation_window[[1]],
      driver_target_u_v = driver_row$target_u_v[[1]]
    )
  )
}

# Save the raw/flagged settling-time points.
write_effective_settling_time_points <- function(path, points_df) {
  points_df %>%
    write_csv(path)
}

# Save the target/global settling-time summaries.
write_effective_settling_time_summary <- function(path, summary_df) {
  summary_df %>%
    write_csv(path)
}

# Save a two-panel diagnostics figure without adding a new plotting dependency.
save_effective_settling_time_diagnostics <- function(
  path,
  points_df,
  summary_df,
  analysis_plot_theme
) {
  target_summary <- summary_df %>%
    filter(level == "target")

  points_plot_df <- points_df %>%
    filter(is.finite(u), is.finite(t_ms), sequence_id > 1, u > 0) %>%
    mutate(
      target_label = factor(sprintf("u = %.3f V", u)),
      point_status = if_else(keep_for_recommendation, "Retained", "Excluded")
    )

  target_line_half_width <- 0.18

  sequence_plot <- ggplot(
    points_plot_df,
    aes(x = sequence_id, y = t_ms, color = target_label, shape = point_status)
  ) +
    geom_point(alpha = 0.9, size = 1.5) +
    scale_color_npg(name = "Target level") +
    scale_shape_manual(values = c(Retained = 16, Excluded = 4), name = "Point status") +
    labs(
      x = "Observation index",
      y = "Measured settling time [ms]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 2)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", size = 2))
    ) +
    analysis_plot_theme

  target_plot <- ggplot(
    points_plot_df,
    aes(x = u, y = t_ms, color = target_label, shape = point_status)
  ) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0, seed = 123),
      alpha = 0.9,
      size = 1.5
    ) +
    geom_segment(
      data = target_summary,
      aes(
        x = target_u_v - target_line_half_width,
        xend = target_u_v + target_line_half_width,
        y = trimmed_median_t_ms,
        yend = trimmed_median_t_ms
      ),
      inherit.aes = FALSE,
      color = "#1F1F1F",
      linewidth = 0.45
    ) +
    geom_segment(
      data = target_summary,
      aes(
        x = target_u_v - target_line_half_width,
        xend = target_u_v + target_line_half_width,
        y = trimmed_p95_t_ms,
        yend = trimmed_p95_t_ms
      ),
      inherit.aes = FALSE,
      color = "#1F1F1F",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    scale_color_npg(name = "Target level") +
    scale_shape_manual(values = c(Retained = 16, Excluded = 4), name = "Point status") +
    labs(
      x = "Target voltage [V]",
      y = "Measured settling time [ms]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 2)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", size = 2))
    ) +
    analysis_plot_theme

  png(
    filename = path,
    width = 180,
    height = 180,
    units = "mm",
    res = 300
  )
  on.exit(dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
  print(
    sequence_plot,
    vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  print(
    target_plot,
    vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1)
  )
}

# Print a concise console summary with the recommended hold time.
print_effective_settling_time_summary <- function(summary_df, points_df, control_window_info) {
  target_summary <- summary_df %>%
    filter(level == "target") %>%
    arrange(target_u_v)
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)

  cat("\nEffective settling-time analysis\n")
  if (isTRUE(control_window_info$available)) {
    cat(sprintf(
      "Control window: [%.6f, %.6f] V; excitation window: [%.6f, %.6f] V\n",
      control_window_info$window_low_v,
      control_window_info$window_high_v,
      control_window_info$excitation_low_v,
      control_window_info$excitation_high_v
    ))
  } else {
    cat("Control-window annotation unavailable: recommended_control_window.csv was not found or had no global row.\n")
  }

  cat(sprintf(
    "Excluded rows: %d total (%d initial warmup, %d startup/non-design, %d robust outliers)\n",
    sum(!points_df$keep_for_recommendation),
    sum(points_df$exclude_reason == "initial_warmup_row", na.rm = TRUE),
    sum(points_df$exclude_reason == "startup_or_non_design", na.rm = TRUE),
    sum(points_df$exclude_reason == "robust_outlier", na.rm = TRUE)
  ))

  cat("Per-target retained 95th percentiles:\n")
  purrr::pwalk(
    list(
      target_summary$target_u_v,
      target_summary$trimmed_p95_t_ms,
      target_summary$relation_to_control_window,
      target_summary$relation_to_excitation_window
    ),
    function(target_u_v, trimmed_p95_t_ms, relation_ctrl, relation_exc) {
      cat(sprintf(
        "  u = %.3f V: p95 = %.2f ms (control window: %s, excitation window: %s)\n",
        target_u_v,
        trimmed_p95_t_ms,
        relation_ctrl %||% "NA",
        relation_exc %||% "NA"
      ))
    }
  )

  cat(sprintf(
    "Recommended effective settling time: %.2f ms\n",
    global_row$recommended_t_settle_ms[[1]]
  ))
  cat(sprintf(
    "Recommended excitation hold time: %.4f s\n",
    global_row$recommended_t_hold_s[[1]]
  ))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

# Load helpers and initialize the shared analysis context.
script_dir <- find_entry_script_dir()
source_analysis_helper(script_dir, "setup.R")
load_analysis_packages()

analysis_context <- initialize_analysis_context()
analysis_plot_theme <- create_analysis_plot_theme()

input_path <- file.path(
  analysis_context$input_dir,
  "settling_time_10_cycl_ma_p_0.005V.csv"
)
summary_path <- file.path(
  analysis_context$output_dir,
  "effective_settling_time_summary.csv"
)
points_path <- file.path(
  analysis_context$output_dir,
  "effective_settling_time_points.csv"
)
plot_path <- file.path(
  analysis_context$output_dir,
  "effective_settling_time_diagnostics.png"
)

settling_time_df <- load_settling_time_data(input_path)
control_window_info <- read_control_window_annotation(analysis_context$output_dir)
points_df <- flag_settling_time_points(settling_time_df)
summary_df <- create_effective_settling_time_summary(points_df, control_window_info)

write_effective_settling_time_summary(summary_path, summary_df)
write_effective_settling_time_points(points_path, points_df)
save_effective_settling_time_diagnostics(
  plot_path,
  points_df,
  summary_df,
  analysis_plot_theme
)
print_effective_settling_time_summary(summary_df, points_df, control_window_info)
