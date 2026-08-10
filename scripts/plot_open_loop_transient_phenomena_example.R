#!/usr/bin/env Rscript

# Plot representative upward and downward holds from the 4 s excitation
# experiment. The figure illustrates pressure-reference transients and the
# later displacement relaxation that motivates creep compensation.

script_args_full <- commandArgs(trailingOnly = FALSE)
script_match <- script_args_full[grepl("^--file=", script_args_full)]
script_dir <- if (length(script_match) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_match[[1]])))
} else {
  file.path(getwd(), "scripts")
}
source(file.path(script_dir, "lib", "setup.R"))

script_args <- commandArgs(trailingOnly = TRUE)
input_csv <- "experiment/excitation_experiment_100hz_15min_4s_hold.csv"
output_png <- "analysis_outputs/open_loop_transient_phenomena_example.png"
for (arg in script_args) {
  pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  if (length(pieces) < 2) {
    stop(sprintf("Arguments must use --name=value: %s", arg))
  }
  if (pieces[[1]] == "input_csv") {
    input_csv <- paste(pieces[-1], collapse = "=")
  } else if (pieces[[1]] == "output_png") {
    output_png <- paste(pieces[-1], collapse = "=")
  } else {
    stop(sprintf("Unknown argument: %s", arg))
  }
}

load_analysis_packages()
repo_dir <- find_repo_dir()
input_path <- file.path(repo_dir, input_csv)
output_path <- file.path(repo_dir, output_png)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop(sprintf("Missing excitation file: %s", input_path))
}
raw_df <- readr::read_csv(input_path, show_col_types = FALSE, trim_ws = TRUE)
names(raw_df) <- trimws(names(raw_df))
required <- c("u[V]", "s[mm]", "p[bar]", "t[ms]")
if (!identical(names(raw_df), required)) {
  stop(sprintf("Expected columns: %s", paste(required, collapse = ", ")))
}

df <- raw_df %>%
  dplyr::transmute(
    sample_id = dplyr::row_number(),
    u_cmd_v = as.numeric(.data[["u[V]"]]),
    s_meas_mm = as.numeric(.data[["s[mm]"]]),
    p_meas_bar = as.numeric(.data[["p[bar]"]]),
    loop_dt_ms = as.numeric(.data[["t[ms]"]])
  )
if (any(!is.finite(unlist(df)))) {
  stop("The excitation file contains non-finite values.")
}

valid_dt_ms <- df$loop_dt_ms[
  df$sample_id > 2L & df$loop_dt_ms > 0 & df$loop_dt_ms < 1000
]
if (length(valid_dt_ms) == 0) {
  stop("Could not infer the sampling period.")
}
ts_s <- stats::median(valid_dt_ms) / 1000

df <- df %>%
  dplyr::mutate(
    hold_id = cumsum(
      dplyr::row_number() == 1L |
        abs(u_cmd_v - dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v))) >
          1e-12
    )
  )

hold_summary <- df %>%
  dplyr::group_by(hold_id) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    u_cmd_v = dplyr::first(u_cmd_v),
    s_start_mm = dplyr::first(s_meas_mm),
    p_start_bar = dplyr::first(p_meas_bar),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    u_prev_v = dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v)),
    delta_u_v = u_cmd_v - u_prev_v,
    branch = dplyr::case_when(
      delta_u_v > 1e-12 ~ "up",
      delta_u_v < -1e-12 ~ "down",
      TRUE ~ "start"
    )
  )

modal_n <- hold_summary %>%
  dplyr::filter(branch %in% c("up", "down"), n_samples > 1L) %>%
  dplyr::count(n_samples, sort = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::pull(n_samples)
if (length(modal_n) != 1 || modal_n < 300) {
  stop("Could not identify full 4 s holds.")
}

df <- df %>%
  dplyr::left_join(hold_summary %>% dplyr::select(-u_cmd_v), by = "hold_id") %>%
  dplyr::group_by(hold_id) %>%
  dplyr::mutate(
    sample_in_hold = dplyr::row_number(),
    time_since_step_s = (sample_in_hold - 1L) * ts_s
  ) %>%
  dplyr::ungroup()

tail_n <- max(5L, ceiling(0.05 * modal_n))
hold_metrics <- df %>%
  dplyr::group_by(
    hold_id,
    branch,
    n_samples,
    u_cmd_v,
    u_prev_v,
    delta_u_v,
    s_start_mm,
    p_start_bar
  ) %>%
  dplyr::summarise(
    s_terminal_mm = mean(utils::tail(s_meas_mm, tail_n)),
    p_terminal_bar = mean(utils::tail(p_meas_bar, tail_n)),
    terminal_s_sd_mm = stats::sd(utils::tail(s_meas_mm, tail_n)),
    terminal_p_sd_bar = stats::sd(utils::tail(p_meas_bar, tail_n)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    displacement_step_mm = s_terminal_mm - s_start_mm,
    pressure_step_bar = p_terminal_bar - p_start_bar,
    full_hold = n_samples == modal_n,
    accepted = branch %in%
      c("up", "down") &
      full_hold &
      abs(displacement_step_mm) >= 0.50 &
      abs(pressure_step_bar) >= 0.02 &
      is.finite(terminal_s_sd_mm) &
      is.finite(terminal_p_sd_bar) &
      terminal_s_sd_mm <= 0.60 &
      terminal_p_sd_bar <= 0.12
  )

candidate_ids <- hold_metrics %>%
  dplyr::filter(accepted) %>%
  dplyr::group_by(branch) %>%
  dplyr::mutate(
    step_magnitude_distance = abs(
      abs(displacement_step_mm) - stats::median(abs(displacement_step_mm))
    )
  ) %>%
  dplyr::arrange(step_magnitude_distance, hold_id, .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
if (nrow(candidate_ids) != 2) {
  stop("Could not find one accepted representative hold per branch.")
}

selected <- df %>%
  dplyr::inner_join(
    candidate_ids %>%
      dplyr::select(
        hold_id,
        branch,
        s_terminal_mm,
        p_terminal_bar,
        displacement_step_mm,
        pressure_step_bar
      ),
    by = c("hold_id", "branch")
  ) %>%
  dplyr::filter(sample_in_hold <= modal_n) %>%
  dplyr::mutate(
    displacement_progress_normalized = (s_meas_mm - s_start_mm) /
      displacement_step_mm,
    pressure_progress_normalized = (p_meas_bar - p_start_bar) /
      pressure_step_bar,
    displacement_progress_display = dplyr::if_else(
      branch == "down",
      -displacement_progress_normalized,
      displacement_progress_normalized
    ),
    pressure_progress_display = dplyr::if_else(
      branch == "down",
      -pressure_progress_normalized,
      pressure_progress_normalized
    ),
    branch_label = factor(
      dplyr::recode(branch, up = "Upward step", down = "Downward step"),
      levels = c("Upward step", "Downward step")
    )
  )

dynamic_summary_path <- file.path(
  repo_dir,
  "analysis_outputs/open_loop_dynamic_pressure_reference_summary.csv"
)
creep_summary_path <- file.path(
  repo_dir,
  "analysis_outputs/open_loop_creep_summary.csv"
)
if (!file.exists(dynamic_summary_path) || !file.exists(creep_summary_path)) {
  stop(
    "Run the open-loop creep and dynamic-pressure analyses before this plot script."
  )
}
dynamic_summary <- readr::read_csv(dynamic_summary_path, show_col_types = FALSE)
creep_summary <- readr::read_csv(creep_summary_path, show_col_types = FALSE)

selected_models <- dynamic_summary %>%
  dplyr::filter(selected_model, branch %in% c("up", "down")) %>%
  dplyr::transmute(branch, theta_s, tau_fast_s, tau_slow_s, weight_fast)
creep_onsets <- creep_summary %>%
  dplyr::filter(selected_model, branch %in% c("up", "down")) %>%
  dplyr::transmute(branch, creep_onset_s = pressure_plateau_start_s)

model_progress <- function(t, theta_s, tau_fast_s, tau_slow_s, weight_fast) {
  elapsed <- pmax(t - theta_s, 0)
  response <- 1 -
    weight_fast * exp(-elapsed / tau_fast_s) -
    (1 - weight_fast) * exp(-elapsed / tau_slow_s)
  ifelse(t < theta_s, 0, response)
}

region_limits <- selected_models %>%
  dplyr::left_join(creep_onsets, by = "branch") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    pressure_95_s = {
      grid <- seq(0, 3.99, by = 0.001)
      progress <- model_progress(
        grid,
        theta_s,
        tau_fast_s,
        tau_slow_s,
        weight_fast
      )
      if (any(progress >= 0.95)) grid[which(progress >= 0.95)[1]] else 3.99
    },
    calibrated_horizon_s = 3.99
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    branch_label = factor(
      dplyr::recode(branch, up = "Upward step", down = "Downward step"),
      levels = c("Upward step", "Downward step")
    )
  )

plot_df <- selected %>%
  dplyr::left_join(
    region_limits %>%
      dplyr::select(branch, theta_s, pressure_95_s, creep_onset_s),
    by = "branch"
  )
region_df <- region_limits %>%
  dplyr::transmute(
    branch,
    branch_label,
    pressure_start = 0,
    pressure_end = pressure_95_s,
    creep_start = creep_onset_s,
    creep_end = 3.99
  )

vline_df <- region_limits %>%
  tidyr::pivot_longer(
    c(theta_s, pressure_95_s, creep_onset_s),
    names_to = "marker",
    values_to = "time_s"
  ) %>%
  dplyr::mutate(
    label = dplyr::recode(
      marker,
      theta_s = "Pressure lag onset",
      pressure_95_s = "95% pressure",
      creep_onset_s = "Creep onset"
    ),
    linetype = label,
    branch_label = factor(
      dplyr::recode(branch, up = "Upward step", down = "Downward step"),
      levels = c("Upward step", "Downward step")
    )
  )

# geom_vline creates vertical legend keys. Keep the event lines vertical in the
# panels, but use clipped horizontal dummy segments for readable legend keys.
horizontal_legend_df <- tibble::tibble(
  linetype = c("Pressure lag onset", "95% pressure", "Creep onset"),
  x = 4.01,
  xend = 4.02,
  y = 0,
  yend = 0
)

region_labels <- dplyr::bind_rows(
  region_df %>%
    dplyr::mutate(
      x = pressure_end / 2,
      y = dplyr::if_else(branch == "up", 0.20, -0.52),
      label = "Pressure\ntransient",
      angle = dplyr::if_else(branch == "up", 90, 0)
    ) %>%
    dplyr::select(branch_label, x, y, label, angle),
  region_df %>%
    dplyr::mutate(
      x = (creep_start + creep_end) / 2,
      y = dplyr::if_else(branch == "up", 0.12, -0.82),
      label = "Post-plateau\ndisplacement creep",
      angle = 0
    ) %>%
    dplyr::select(branch_label, x, y, label, angle)
)

plot_theme <- create_wide_analysis_plot_theme() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 8.5),
    legend.title = ggplot2::element_text(size = 8.5),
    legend.box.spacing = grid::unit(0.5, "mm"),
    strip.text = ggplot2::element_text(face = "plain", size = 9.5),
    plot.margin = ggplot2::margin(1, 1, 7, 1, "mm")
  )

plot <- ggplot2::ggplot(plot_df, ggplot2::aes(time_since_step_s)) +
  ggplot2::geom_rect(
    data = region_df,
    ggplot2::aes(
      xmin = pressure_start,
      xmax = pressure_end,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "#E69F00",
    alpha = 0.10
  ) +
  ggplot2::geom_rect(
    data = region_df,
    ggplot2::aes(xmin = creep_start, xmax = creep_end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#56B4E9",
    alpha = 0.12
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = displacement_progress_display,
      color = "Measured displacement"
    ),
    linewidth = 0.45
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = pressure_progress_display,
      color = "Measured pressure"
    ),
    linewidth = 0.45
  ) +
  ggplot2::geom_vline(
    data = vline_df,
    ggplot2::aes(xintercept = time_s, linetype = linetype),
    color = "grey25",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  ggplot2::geom_segment(
    data = horizontal_legend_df,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, linetype = linetype),
    inherit.aes = FALSE,
    color = "grey25",
    linewidth = 0.45,
    alpha = 0,
    show.legend = TRUE
  ) +
  ggplot2::geom_label(
    data = region_labels,
    ggplot2::aes(x = x, y = y, label = label, angle = angle),
    inherit.aes = FALSE,
    size = 3.0,
    lineheight = 0.9,
    fontface = "plain",
    color = "grey20",
    fill = scales::alpha("white", 0.75),
    linewidth = 0
  ) +
  ggplot2::facet_wrap(~branch_label, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(
    name = NULL,
    labels = c(
      "Measured displacement" = "Displacement",
      "Measured pressure" = "Pressure"
    ),
    values = c(
      "Measured displacement" = "#0072B2",
      "Measured pressure" = "#D55E00"
    )
  ) +
  ggplot2::scale_linetype_manual(
    name = NULL,
    labels = c(
      "Pressure lag onset" = "Pressure lag",
      "95% pressure" = "95% pressure",
      "Creep onset" = "Creep onset"
    ),
    values = c(
      "Pressure lag onset" = "dashed",
      "95% pressure" = "dotted",
      "Creep onset" = "dotdash"
    )
  ) +
  ggplot2::guides(
    linetype = ggplot2::guide_legend(
      keywidth = grid::unit(12, "mm"),
      override.aes = list(alpha = 1, color = "grey25", linewidth = 0.45)
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 4.05),
    breaks = seq(0, 4, by = 0.5),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    name = "Normalized response [-]",
    breaks = seq(-1, 1, by = 0.2),
    labels = function(x) abs(x),
    expand = c(0, 0)
  ) +
  ggplot2::labs(x = "Time since voltage step [s]") +
  plot_theme

ggplot2::ggsave(
  output_path,
  plot,
  width = 17,
  height = 9.5625,
  units = "cm",
  dpi = 600
)
message(sprintf("Saved transient-phenomena example plot to: %s", output_path))
message(sprintf(
  "Selected holds: upward %d (%.3f V -> %.3f V), downward %d (%.3f V -> %.3f V)",
  candidate_ids$hold_id[candidate_ids$branch == "up"],
  candidate_ids$u_prev_v[candidate_ids$branch == "up"],
  candidate_ids$u_cmd_v[candidate_ids$branch == "up"],
  candidate_ids$hold_id[candidate_ids$branch == "down"],
  candidate_ids$u_prev_v[candidate_ids$branch == "down"],
  candidate_ids$u_cmd_v[candidate_ids$branch == "down"]
))
