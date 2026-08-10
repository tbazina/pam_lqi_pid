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

# Use the manufacturer-provided Schmalz VP8 inverse calibration curve with
# a nonnegative pressure clamp for this experiment family.
manufacturer_calibration <- function() {
  list(
    source = "Schmalz VP8 datasheet",
    slope_bar_per_v = 0.9,
    intercept_bar = -1.0,
    pressure_floor_bar = 0,
    pressure_ceiling_bar = 8,
    sensor_voltage_min_v = 0,
    sensor_voltage_max_v = 10,
    zero_pressure_voltage_v = 1 / 0.9
  )
}

# Return `NA` instead of failing when correlation is undefined.
safe_correlation <- function(x, y) {
  if (length(x) < 2 || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }

  stats::cor(x, y)
}

# Parse run parameters from one raw pressure-voltage experiment file.
parse_raw_voltage_run_parameters <- function(path) {
  stripped_path <- sub("_p_voltage\\.csv$", ".csv", path)
  params <- parse_run_parameters(stripped_path)
  params$file <- basename(path)
  params
}

# Discover only raw pressure-voltage experiment files.
discover_raw_voltage_runs <- function(input_dir) {
  raw_files <- list.files(
    input_dir,
    pattern = "_p_voltage\\.csv$",
    full.names = TRUE
  )

  if (length(raw_files) == 0) {
    stop("No *_p_voltage.csv files were found for pressure calibration.")
  }

  tibble(raw_path = raw_files) %>%
    mutate(
      pair_key = sub("_p_voltage\\.csv$", "", basename(raw_path)),
      params = purrr::map(raw_path, parse_raw_voltage_run_parameters)
    ) %>%
    tidyr::unnest(params) %>%
    mutate(
      run_id = sprintf("step_%0.2fV_settle_%dms", step_v, settle_ms),
      run_id = factor(run_id, levels = run_id)
    ) %>%
    arrange(step_v, settle_ms)
}

# Read the global recommended control-window upper bound used for annotations.
read_window_high_v <- function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }

  window_df <- read_csv(path, show_col_types = FALSE)
  global_row <- window_df %>%
    filter(level == "global") %>%
    slice(1)

  if (nrow(global_row) == 0 || !is.finite(global_row$window_high_v[[1]])) {
    return(NA_real_)
  }

  as.numeric(global_row$window_high_v[[1]])
}

# Convert raw sensor voltage to nonnegative pressure in bar.
calibrate_pressure_bar <- function(p_raw_v, calibration_model) {
  pmax(
    calibration_model$pressure_floor_bar,
    calibration_model$intercept_bar + calibration_model$slope_bar_per_v * p_raw_v
  )
}

# Load one raw pressure-voltage file and apply the shared cycle preprocessing.
load_raw_voltage_run <- function(path, params) {
  raw_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  raw_df %>%
    rename(
      u = `u[V]`,
      p = `p[V]`,
      s = `s[mm]`
    ) %>%
    mutate(row_id = row_number()) %>%
    bind_cols(params[rep(1, nrow(raw_df)), c("file", "step_v", "settle_ms", "run_id")]) %>%
    infer_direction_and_cycles() %>%
    prune_incomplete_segments()
}

# Collapse one run to direction-by-voltage settled means.
aggregate_voltage_pressure_points <- function(df, calibration_model) {
  df %>%
    group_by(run_id, file, step_v, settle_ms, direction, u) %>%
    summarise(
      p_raw_v_mean = mean(p, na.rm = TRUE),
      p_raw_v_min = min(p, na.rm = TRUE),
      p_raw_v_max = max(p, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    mutate(
      p_bar_mean = calibrate_pressure_bar(p_raw_v_mean, calibration_model),
      p_bar_min = calibrate_pressure_bar(p_raw_v_min, calibration_model),
      p_bar_max = calibrate_pressure_bar(p_raw_v_max, calibration_model),
      clamped = p_bar_mean <= calibration_model$pressure_floor_bar + 1e-12,
      in_sensor_voltage_range = dplyr::between(
        p_raw_v_mean,
        calibration_model$sensor_voltage_min_v,
        calibration_model$sensor_voltage_max_v
      ),
      in_sensor_pressure_range = dplyr::between(
        p_bar_mean,
        calibration_model$pressure_floor_bar,
        calibration_model$pressure_ceiling_bar
      ),
      above_zero_pressure_threshold = p_raw_v_mean >= calibration_model$zero_pressure_voltage_v
    ) %>%
    arrange(run_id, direction, u)
}

# Save the settled voltage-to-pressure points used for review and plotting.
write_calibration_points <- function(path, points_df) {
  points_df %>%
    transmute(
      run_id,
      file,
      step_v,
      settle_ms,
      direction,
      u,
      n_samples,
      p_raw_v_mean,
      p_raw_v_min,
      p_raw_v_max,
      p_bar_mean,
      p_bar_min,
      p_bar_max,
      clamped,
      in_sensor_voltage_range,
      in_sensor_pressure_range,
      above_zero_pressure_threshold
    ) %>%
    write_csv(path)
}

# Summarize per-run and global raw-signal/calibration quality.
create_calibration_summary <- function(
  points_df,
  calibration_model,
  window_high_v
) {
  summarize_group <- function(df, level_label, direction_value = NA_character_) {
    tibble(
      level = level_label,
      run_id = as.character(df$run_id[[1]]),
      file = as.character(df$file[[1]]),
      step_v = as.numeric(df$step_v[[1]]),
      settle_ms = as.integer(df$settle_ms[[1]]),
      direction = direction_value,
      calibration_source = calibration_model$source,
      slope_bar_per_v = calibration_model$slope_bar_per_v,
      intercept_bar = calibration_model$intercept_bar,
      pressure_floor_bar = calibration_model$pressure_floor_bar,
      n_points = nrow(df),
      u_min_v = min(df$u, na.rm = TRUE),
      u_max_v = max(df$u, na.rm = TRUE),
      window_high_v_used = window_high_v,
      raw_voltage_mean_v = mean(df$p_raw_v_mean, na.rm = TRUE),
      raw_voltage_min_v = min(df$p_raw_v_mean, na.rm = TRUE),
      raw_voltage_max_v = max(df$p_raw_v_mean, na.rm = TRUE),
      calibrated_pressure_mean_bar = mean(df$p_bar_mean, na.rm = TRUE),
      calibrated_pressure_min_bar = min(df$p_bar_mean, na.rm = TRUE),
      calibrated_pressure_max_bar = max(df$p_bar_mean, na.rm = TRUE),
      calibrated_pressure_span_bar = max(df$p_bar_mean, na.rm = TRUE) -
        min(df$p_bar_mean, na.rm = TRUE),
      clamped_point_frac = mean(df$clamped, na.rm = TRUE),
      above_zero_pressure_threshold_frac = mean(df$above_zero_pressure_threshold, na.rm = TRUE),
      in_sensor_voltage_range_frac = mean(df$in_sensor_voltage_range, na.rm = TRUE),
      in_sensor_pressure_range_frac = mean(df$in_sensor_pressure_range, na.rm = TRUE),
      monotonic_correlation_u_vs_p_bar = safe_correlation(df$u, df$p_bar_mean)
    )
  }

  by_direction <- points_df %>%
    group_by(run_id, file, step_v, settle_ms, direction) %>%
    group_split(.keep = TRUE) %>%
    purrr::map_dfr(~summarize_group(
      .x,
      level_label = "run_direction",
      direction_value = as.character(.x$direction[[1]])
    ))

  by_run <- points_df %>%
    group_by(run_id, file, step_v, settle_ms) %>%
    group_split(.keep = TRUE) %>%
    purrr::map_dfr(~summarize_group(
      .x,
      level_label = "run",
      direction_value = "combined"
    ))

  global_row <- summarize_group(
    points_df,
    level_label = "global",
    direction_value = "combined"
  ) %>%
    mutate(
      run_id = "global_manufacturer_curve",
      file = NA_character_,
      step_v = NA_real_,
      settle_ms = NA_integer_
    )

  bind_rows(by_direction, by_run, global_row) %>%
    mutate(
      voltage_range_ok = in_sensor_voltage_range_frac >= 0.999,
      pressure_range_ok = in_sensor_pressure_range_frac >= 0.999,
      monotonic_ok = is.na(monotonic_correlation_u_vs_p_bar) |
        monotonic_correlation_u_vs_p_bar >= 0.98,
      validation_ok = voltage_range_ok & pressure_range_ok & monotonic_ok,
      review_needed = !validation_ok
    )
}

# Save the summary CSV describing run-level and global calibration quality.
write_calibration_summary <- function(path, summary_df) {
  summary_df %>%
    write_csv(path)
}

# Save validation statistics for later inspection.
write_voltage_pressure_validation_summary <- function(path, summary_df) {
  summary_df %>%
    write_csv(path)
}

# Create a validation summary focused on operating-window settled points.
create_voltage_pressure_validation_summary <- function(
  points_df,
  calibration_model,
  window_high_v
) {
  validation_points <- points_df
  if (is.finite(window_high_v)) {
    validation_points <- validation_points %>%
      filter(u <= window_high_v)
  }

  if (nrow(validation_points) == 0) {
    stop("No settled raw-voltage points remain for validation.")
  }

  create_calibration_summary(
    points_df = validation_points,
    calibration_model = calibration_model,
    window_high_v = window_high_v
  )
}

# Plot the manufacturer map and the measured settled raw-voltage means.
save_calibration_plot <- function(
  path,
  points_df,
  calibration_model,
  analysis_plot_theme
) {
  x_range <- range(points_df$p_raw_v_mean, na.rm = TRUE)
  map_line <- tibble(
    p_raw_v_mean = seq(x_range[1], x_range[2], length.out = 300)
  ) %>%
    mutate(
      p_bar_fit = calibrate_pressure_bar(p_raw_v_mean, calibration_model)
    )

  plot_obj <- ggplot(points_df, aes(x = p_raw_v_mean, y = p_bar_mean, color = run_id)) +
    geom_vline(
      xintercept = calibration_model$zero_pressure_voltage_v,
      color = "#1F1F1F",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_point(alpha = 0.8, size = 1.4) +
    geom_line(
      data = map_line,
      aes(x = p_raw_v_mean, y = p_bar_fit),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.6
    ) +
    scale_color_npg(name = "Run", labels = humanize_run_id) +
    labs(
      x = "Raw pressure sensor voltage [V]",
      y = "Calibrated pressure [bar]"
    ) +
    analysis_plot_theme

  ggsave(
    filename = path,
    plot = plot_obj,
    width = 160,
    height = 120,
    units = "mm",
    dpi = 300
  )
}

# Save a faceted raw-voltage/pressure-vs-command validation panel.
save_voltage_pressure_validation_plot <- function(
  path,
  points_df,
  calibration_model,
  window_high_v,
  analysis_plot_theme
) {
  plot_df <- points_df %>%
    dplyr::select(run_id, direction, u, p_raw_v_mean, p_bar_mean) %>%
    pivot_longer(
      cols = c(p_raw_v_mean, p_bar_mean),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = recode(
        metric,
        p_raw_v_mean = "Raw pressure sensor voltage [V]",
        p_bar_mean = "Calibrated pressure [bar]"
      ),
      metric = factor(
        metric,
        levels = c("Raw pressure sensor voltage [V]", "Calibrated pressure [bar]")
      )
    )

  hline_df <- tibble(
    metric = factor(
      c("Raw pressure sensor voltage [V]", "Calibrated pressure [bar]"),
      levels = levels(plot_df$metric)
    ),
    yintercept = c(
      calibration_model$zero_pressure_voltage_v,
      calibration_model$pressure_floor_bar
    )
  )

  plot_obj <- ggplot(plot_df, aes(x = u, y = value, color = direction, group = direction)) +
    {
      if (is.finite(window_high_v)) {
        geom_vline(
          xintercept = window_high_v,
          color = "#1F1F1F",
          linetype = "dashed",
          linewidth = 0.35
        )
      }
    } +
    geom_hline(
      data = hline_df,
      aes(yintercept = yintercept),
      inherit.aes = FALSE,
      color = "#1F1F1F",
      linetype = "dotted",
      linewidth = 0.3
    ) +
    geom_point(alpha = 0.85, size = 1.1) +
    geom_line(linewidth = 0.45) +
    facet_grid(metric ~ run_id, scales = "free_y", labeller = as_labeller(humanize_run_id)) +
    scale_color_npg(name = "Direction", labels = c(up = "Upward", down = "Downward")) +
    labs(
      x = "Command voltage [V]",
      y = NULL
    ) +
    analysis_plot_theme

  ggsave(
    filename = path,
    plot = plot_obj,
    width = 180,
    height = 120,
    units = "mm",
    dpi = 300
  )
}

# Append calibrated pressure in bar to one raw sensor-voltage CSV file.
write_calibrated_experiment_file <- function(raw_path, output_path, calibration_model) {
  raw_df <- read_csv(raw_path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  raw_df <- raw_df %>%
    mutate(`p[bar]` = calibrate_pressure_bar(`p[V]`, calibration_model))

  write_csv(raw_df, output_path)
}

# Print a concise final calibration summary and recommendation guidance.
print_calibration_summary <- function(
  summary_df,
  validation_summary_df,
  output_dir,
  calibration_model
) {
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)
  validation_global_row <- validation_summary_df %>%
    filter(level == "global") %>%
    slice(1)

  cat("\nPressure sensor calibration\n")
  cat(sprintf(
    "Manufacturer map: p[bar] = max(0, %.6f * p[V] %+0.6f)\n",
    calibration_model$slope_bar_per_v,
    calibration_model$intercept_bar
  ))
  if (is.finite(global_row$window_high_v_used[[1]])) {
    cat(sprintf(
      "Operating-window annotation: u[V] <= %.6f\n",
      global_row$window_high_v_used[[1]]
    ))
  }
  cat(sprintf(
    "Global settled points: raw voltage %.3f to %.3f V, calibrated pressure %.3f to %.3f bar.\n",
    global_row$raw_voltage_min_v[[1]],
    global_row$raw_voltage_max_v[[1]],
    global_row$calibrated_pressure_min_bar[[1]],
    global_row$calibrated_pressure_max_bar[[1]]
  ))
  cat(sprintf(
    "Clamp usage: %.2f%% of settled points are at zero pressure after clamping.\n",
    100 * global_row$clamped_point_frac[[1]]
  ))
  cat(sprintf(
    "Monotonicity check: corr(u, p[bar]) = %.6f",
    validation_global_row$monotonic_correlation_u_vs_p_bar[[1]]
  ))
  if (isTRUE(validation_global_row$validation_ok[[1]])) {
    cat("; validation result: raw-voltage data are consistent with the manufacturer map.\n")
  } else {
    cat("; validation result: review the raw-voltage range or monotonicity.\n")
  }
  cat(sprintf("Outputs written to: %s\n", output_dir))
}

# Load helper modules in dependency order, then initialize shared context.
script_dir <- find_entry_script_dir()

source_analysis_helper(script_dir, "setup.R")
load_analysis_packages()
source_analysis_helper(script_dir, "io_preprocess.R")

analysis_context <- initialize_analysis_context()
analysis_plot_theme <- create_analysis_plot_theme()

input_dir <- analysis_context$input_dir
output_dir <- analysis_context$output_dir

calibration_model <- manufacturer_calibration()
window_high_v <- read_window_high_v(
  file.path(output_dir, "recommended_control_window.csv")
)
runs_df <- discover_raw_voltage_runs(input_dir)

points_df <- runs_df %>%
  mutate(
    raw_data = purrr::map(
      seq_len(n()),
      ~load_raw_voltage_run(
        path = raw_path[.x],
        params = runs_df[.x, ]
      )
    ),
    settled_points = purrr::map(
      raw_data,
      aggregate_voltage_pressure_points,
      calibration_model = calibration_model
    )
  ) %>%
  pull(settled_points) %>%
  bind_rows() %>%
  arrange(run_id, direction, u)

if (nrow(points_df) == 0) {
  stop("No settled raw-voltage points were produced from the *_p_voltage.csv files.")
}

summary_df <- create_calibration_summary(
  points_df = points_df,
  calibration_model = calibration_model,
  window_high_v = window_high_v
)
validation_summary_df <- create_voltage_pressure_validation_summary(
  points_df = points_df,
  calibration_model = calibration_model,
  window_high_v = window_high_v
)

summary_path <- file.path(output_dir, "pressure_sensor_calibration_summary.csv")
points_path <- file.path(output_dir, "pressure_sensor_calibration_points.csv")
plot_path <- file.path(output_dir, "pressure_sensor_calibration_fit.png")
validation_summary_path <- file.path(
  output_dir,
  "pressure_sensor_calibration_validation_summary.csv"
)
validation_plot_path <- file.path(
  output_dir,
  "pressure_sensor_calibration_validation_panel.png"
)

write_calibration_summary(summary_path, summary_df)
write_calibration_points(points_path, points_df)
write_voltage_pressure_validation_summary(
  validation_summary_path,
  validation_summary_df
)
save_calibration_plot(
  plot_path,
  points_df,
  calibration_model,
  analysis_plot_theme
)
save_voltage_pressure_validation_plot(
  validation_plot_path,
  points_df,
  calibration_model,
  window_high_v,
  analysis_plot_theme
)

purrr::walk2(
  .x = runs_df$raw_path,
  .y = file.path(
    input_dir,
    sprintf("%s_p_voltage_bar.csv", runs_df$pair_key)
  ),
  .f = ~write_calibrated_experiment_file(
    raw_path = .x,
    output_path = .y,
    calibration_model = calibration_model
  )
)

print_calibration_summary(
  summary_df = summary_df,
  validation_summary_df = validation_summary_df,
  output_dir = output_dir,
  calibration_model = calibration_model
)
