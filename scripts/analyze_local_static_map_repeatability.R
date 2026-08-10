#!/usr/bin/env Rscript

# Assess local terminal repeatability in the 4 s open-loop excitation run.
# This is a binned static-map diagnostic, not an exact-cycle repeatability test.

VOLTAGE_BIN_WIDTH_V <- 0.10
TERMINAL_TAIL_FRACTION <- 0.05
MIN_REPEATS_PER_BIN <- 3L
# Match the steady-state extraction criteria documented for the active lookup.
MAX_TERMINAL_S_SD_MM <- 0.60
MAX_TERMINAL_S_SLOPE_MM_S <- 1.20
MAX_TERMINAL_P_SD_BAR <- 0.12
MAX_TERMINAL_P_SLOPE_BAR_S <- 0.25

find_entry_script_dir <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_match[[1]]))))
  }
  candidate <- file.path(getwd(), "scripts")
  if (dir.exists(candidate)) {
    normalizePath(candidate)
  } else {
    normalizePath(getwd())
  }
}

script_dir <- find_entry_script_dir()
source(file.path(script_dir, "lib", "setup.R"))
load_analysis_packages()

parse_settings <- function(defaults) {
  settings <- defaults
  for (arg in commandArgs(trailingOnly = TRUE)) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) < 2) {
      stop(sprintf("Arguments must use --name=value: %s", arg))
    }
    key <- pieces[[1]]
    value <- paste(pieces[-1], collapse = "=")
    if (!key %in% names(settings)) {
      stop(sprintf("Unknown setting override: %s", key))
    }
    if (is.integer(settings[[key]])) {
      settings[[key]] <- as.integer(value)
    } else if (is.numeric(settings[[key]])) {
      settings[[key]] <- as.numeric(value)
    } else {
      settings[[key]] <- value
    }
  }
  settings
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else stats::sd(x)
}

safe_mad <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    NA_real_
  } else {
    stats::mad(x, center = stats::median(x), constant = 1.4826)
  }
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    NA_real_
  } else {
    as.numeric(stats::quantile(x, probability, names = FALSE, type = 7))
  }
}

safe_slope <- function(y, time_s) {
  keep <- is.finite(y) & is.finite(time_s)
  if (sum(keep) < 2 || length(unique(time_s[keep])) < 2) {
    return(NA_real_)
  }
  as.numeric(stats::coef(stats::lm(y[keep] ~ time_s[keep]))[[2]])
}

load_excitation <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing excitation experiment: %s", path))
  }
  df <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(df) <- trimws(names(df))
  required <- c("u[V]", "s[mm]", "p[bar]", "t[ms]")
  if (!identical(names(df), required)) {
    stop(sprintf(
      "Expected columns %s in %s",
      paste(required, collapse = ", "),
      path
    ))
  }
  df <- df %>%
    dplyr::transmute(
      sample_id = dplyr::row_number(),
      u_cmd_v = as.numeric(`u[V]`),
      s_meas_mm = as.numeric(`s[mm]`),
      p_meas_bar = as.numeric(`p[bar]`),
      loop_dt_ms = as.numeric(`t[ms]`)
    )
  required_numeric <- c("u_cmd_v", "s_meas_mm", "p_meas_bar", "loop_dt_ms")
  if (any(!is.finite(unlist(df[required_numeric])))) {
    stop("The excitation file contains non-finite values.")
  }
  df
}

infer_sampling_time_s <- function(df) {
  valid <- df$loop_dt_ms[
    df$sample_id > 2L & df$loop_dt_ms > 0 & df$loop_dt_ms < 1000
  ]
  if (length(valid) == 0) {
    stop("Could not infer a sampling period from loop times.")
  }
  stats::median(valid) / 1000
}

segment_holds <- function(df, ts_s) {
  df %>%
    dplyr::mutate(
      hold_id = cumsum(
        dplyr::row_number() == 1L |
          abs(u_cmd_v - dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v))) >
            1e-12
      )
    ) %>%
    dplyr::group_by(hold_id) %>%
    dplyr::mutate(
      sample_in_hold = dplyr::row_number(),
      time_since_step_s = (sample_in_hold - 1L) * ts_s
    ) %>%
    dplyr::ungroup()
}

make_hold_table <- function(segmented_df, ts_s, settings) {
  hold_table <- segmented_df %>%
    dplyr::group_by(hold_id) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      u_cmd_v = dplyr::first(u_cmd_v),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      u_previous_v = dplyr::lag(u_cmd_v),
      delta_u_v = u_cmd_v - u_previous_v,
      branch = dplyr::case_when(
        is.na(delta_u_v) ~ "start",
        delta_u_v > 1e-12 ~ "up",
        delta_u_v < -1e-12 ~ "down",
        TRUE ~ "start"
      )
    )

  modal_counts <- hold_table %>%
    dplyr::filter(branch %in% c("up", "down"), n_samples > 1L) %>%
    dplyr::count(n_samples, sort = TRUE)
  if (nrow(modal_counts) == 0) {
    stop("Could not determine the modal full-hold length.")
  }
  modal_n <- modal_counts$n_samples[[1]]

  hold_rows <- lapply(seq_len(nrow(hold_table)), function(i) {
    hold_info <- hold_table[i, ]
    hold_df <- segmented_df %>% dplyr::filter(hold_id == hold_info$hold_id[[1]])
    tail_n <- max(5L, ceiling(settings$terminal_tail_fraction * nrow(hold_df)))
    tail_df <- utils::tail(hold_df, tail_n)
    is_full <- hold_info$branch[[1]] %in%
      c("up", "down") &&
      abs(hold_info$n_samples[[1]] - modal_n) <=
        settings$full_hold_count_tolerance
    s_terminal <- stats::median(tail_df$s_meas_mm)
    p_terminal <- stats::median(tail_df$p_meas_bar)
    s_sd <- safe_sd(tail_df$s_meas_mm)
    p_sd <- safe_sd(tail_df$p_meas_bar)
    s_slope <- safe_slope(tail_df$s_meas_mm, tail_df$time_since_step_s)
    p_slope <- safe_slope(tail_df$p_meas_bar, tail_df$time_since_step_s)
    terminal_stable <- is.finite(s_sd) &&
      is.finite(p_sd) &&
      is.finite(s_slope) &&
      is.finite(p_slope) &&
      s_sd <= settings$max_terminal_s_sd_mm &&
      abs(s_slope) <= settings$max_terminal_s_slope_mm_s &&
      p_sd <= settings$max_terminal_p_sd_bar &&
      abs(p_slope) <= settings$max_terminal_p_slope_bar_s
    exclusion_reason <- dplyr::case_when(
      hold_info$branch[[1]] == "start" ~ "startup_or_zero_command_hold",
      !is_full ~ "incomplete_hold",
      !terminal_stable ~ "unstable_terminal_tail",
      TRUE ~ ""
    )
    tibble::tibble(
      hold_id = hold_info$hold_id[[1]],
      branch = hold_info$branch[[1]],
      n_samples = hold_info$n_samples[[1]],
      is_full_hold = is_full,
      hold_duration_s = hold_info$n_samples[[1]] * ts_s,
      terminal_tail_n_samples = tail_n,
      terminal_tail_start_s = max(
        0,
        (hold_info$n_samples[[1]] - tail_n) * ts_s
      ),
      u_previous_v = hold_info$u_previous_v[[1]],
      u_cmd_v = hold_info$u_cmd_v[[1]],
      delta_u_v = hold_info$delta_u_v[[1]],
      abs_delta_u_v = abs(hold_info$delta_u_v[[1]]),
      s_start_mm = hold_df$s_meas_mm[[1]],
      s_terminal_mm = s_terminal,
      p_start_bar = hold_df$p_meas_bar[[1]],
      p_terminal_bar = p_terminal,
      s_terminal_sd_mm = s_sd,
      p_terminal_sd_bar = p_sd,
      s_terminal_slope_mm_s = s_slope,
      p_terminal_slope_bar_s = p_slope,
      terminal_stable = terminal_stable,
      eligible_for_repeatability = is_full && terminal_stable,
      exclusion_reason = dplyr::na_if(exclusion_reason, "")
    )
  })

  dplyr::bind_rows(hold_rows) %>%
    dplyr::mutate(
      full_modal_n = as.integer(modal_n),
      voltage_bin_width_v = settings$voltage_bin_width_v,
      voltage_bin_lower_v = floor(
        (u_cmd_v + 1e-10) / settings$voltage_bin_width_v
      ) *
        settings$voltage_bin_width_v,
      voltage_bin_upper_v = voltage_bin_lower_v + settings$voltage_bin_width_v,
      voltage_bin_center_v = (voltage_bin_lower_v + voltage_bin_upper_v) / 2,
      voltage_bin = sprintf(
        "[%.2f, %.2f)",
        voltage_bin_lower_v,
        voltage_bin_upper_v
      )
    )
}

make_local_bin_summary <- function(hold_table, settings) {
  eligible <- hold_table %>%
    dplyr::filter(eligible_for_repeatability, branch %in% c("up", "down"))
  if (nrow(eligible) == 0) {
    stop(
      "No complete stable directional holds remain for repeatability analysis."
    )
  }

  branch_spans <- eligible %>%
    dplyr::group_by(branch) %>%
    dplyr::summarise(
      branch_s_span_mm = max(s_terminal_mm) - min(s_terminal_mm),
      branch_p_span_bar = max(p_terminal_bar) - min(p_terminal_bar),
      .groups = "drop"
    )

  by_bin <- eligible %>%
    dplyr::group_by(
      branch,
      voltage_bin,
      voltage_bin_lower_v,
      voltage_bin_upper_v,
      voltage_bin_center_v
    ) %>%
    dplyr::summarise(
      n_holds = dplyr::n(),
      n_unique_voltage_values = dplyr::n_distinct(u_cmd_v),
      u_min_v = min(u_cmd_v),
      u_median_v = stats::median(u_cmd_v),
      u_max_v = max(u_cmd_v),
      s_terminal_mean_mm = mean(s_terminal_mm),
      s_terminal_median_mm = stats::median(s_terminal_mm),
      s_terminal_sd_mm = safe_sd(s_terminal_mm),
      s_terminal_mad_mm = safe_mad(s_terminal_mm),
      s_terminal_p95_abs_dev_mm = safe_quantile(
        abs(s_terminal_mm - stats::median(s_terminal_mm)),
        0.95
      ),
      p_terminal_mean_bar = mean(p_terminal_bar),
      p_terminal_median_bar = stats::median(p_terminal_bar),
      p_terminal_sd_bar = safe_sd(p_terminal_bar),
      p_terminal_mad_bar = safe_mad(p_terminal_bar),
      p_terminal_p95_abs_dev_bar = safe_quantile(
        abs(p_terminal_bar - stats::median(p_terminal_bar)),
        0.95
      ),
      mean_terminal_s_sd_mm = mean(s_terminal_sd_mm),
      mean_terminal_p_sd_bar = mean(p_terminal_sd_bar),
      mean_terminal_slope_mm_s = mean(s_terminal_slope_mm_s),
      mean_terminal_p_slope_bar_s = mean(p_terminal_slope_bar_s),
      .groups = "drop"
    ) %>%
    dplyr::left_join(branch_spans, by = "branch") %>%
    dplyr::mutate(
      s_repeatability_sd_pct_branch_span = 100 *
        s_terminal_sd_mm /
        branch_s_span_mm,
      p_repeatability_sd_pct_branch_span = 100 *
        p_terminal_sd_bar /
        branch_p_span_bar,
      s_repeatability_p95_pct_branch_span = 100 *
        s_terminal_p95_abs_dev_mm /
        branch_s_span_mm,
      p_repeatability_p95_pct_branch_span = 100 *
        p_terminal_p95_abs_dev_bar /
        branch_p_span_bar,
      sufficient_repeats = n_holds >= settings$min_repeats_per_bin,
      repeatability_assessment = dplyr::case_when(
        !sufficient_repeats ~ "insufficient_repeats",
        s_repeatability_sd_pct_branch_span <= 3 &
          p_repeatability_sd_pct_branch_span <= 1 ~ "likely_sufficient",
        s_repeatability_sd_pct_branch_span <= 5 &
          p_repeatability_sd_pct_branch_span <= 2 ~ "possibly_sufficient",
        TRUE ~ "poor_repeatability_or_local_dynamics"
      )
    ) %>%
    dplyr::arrange(branch, voltage_bin_lower_v)

  by_bin
}

make_overview <- function(hold_table, by_bin, settings, input_path, ts_s) {
  branch_rows <- lapply(c("up", "down"), function(branch_name) {
    all_branch <- hold_table %>% dplyr::filter(branch == branch_name)
    stable_branch <- all_branch %>% dplyr::filter(eligible_for_repeatability)
    eligible_bins <- by_bin %>%
      dplyr::filter(branch == branch_name, sufficient_repeats)
    branch_span_s <- if (nrow(stable_branch) > 0) {
      max(stable_branch$s_terminal_mm) - min(stable_branch$s_terminal_mm)
    } else {
      NA_real_
    }
    branch_span_p <- if (nrow(stable_branch) > 0) {
      max(stable_branch$p_terminal_bar) - min(stable_branch$p_terminal_bar)
    } else {
      NA_real_
    }
    tibble::tibble(
      branch = branch_name,
      input_dataset = input_path,
      voltage_bin_width_v = settings$voltage_bin_width_v,
      nominal_ts_s = ts_s,
      max_terminal_s_sd_mm = settings$max_terminal_s_sd_mm,
      max_terminal_s_slope_mm_s = settings$max_terminal_s_slope_mm_s,
      max_terminal_p_sd_bar = settings$max_terminal_p_sd_bar,
      max_terminal_p_slope_bar_s = settings$max_terminal_p_slope_bar_s,
      modal_hold_samples = unique(hold_table$full_modal_n),
      modal_hold_duration_s = unique(hold_table$full_modal_n) * ts_s,
      total_directional_holds = nrow(all_branch),
      complete_directional_holds = sum(all_branch$is_full_hold),
      stable_directional_holds = nrow(stable_branch),
      excluded_directional_holds = sum(!all_branch$eligible_for_repeatability),
      voltage_bins_present = nrow(
        by_bin %>% dplyr::filter(branch == branch_name)
      ),
      voltage_bins_with_min_repeats = nrow(eligible_bins),
      bins_likely_sufficient = sum(
        eligible_bins$repeatability_assessment == "likely_sufficient"
      ),
      bins_possibly_sufficient = sum(
        eligible_bins$repeatability_assessment == "possibly_sufficient"
      ),
      bins_poor_or_dynamic = sum(
        eligible_bins$repeatability_assessment ==
          "poor_repeatability_or_local_dynamics"
      ),
      branch_s_span_mm = branch_span_s,
      branch_p_span_bar = branch_span_p,
      s_bin_sd_p95_mm = safe_quantile(eligible_bins$s_terminal_sd_mm, 0.95),
      p_bin_sd_p95_bar = safe_quantile(eligible_bins$p_terminal_sd_bar, 0.95),
      s_repeatability_sd_p95_pct_branch_span = safe_quantile(
        eligible_bins$s_repeatability_sd_pct_branch_span,
        0.95
      ),
      p_repeatability_sd_p95_pct_branch_span = safe_quantile(
        eligible_bins$p_repeatability_sd_pct_branch_span,
        0.95
      ),
      assessment = dplyr::case_when(
        nrow(eligible_bins) == 0 ~ "insufficient_repeated_voltage_bins",
        safe_quantile(eligible_bins$s_repeatability_sd_pct_branch_span, 0.95) <=
          3 &&
          safe_quantile(
            eligible_bins$p_repeatability_sd_pct_branch_span,
            0.95
          ) <=
            1 ~ "likely_sufficient",
        safe_quantile(eligible_bins$s_repeatability_sd_pct_branch_span, 0.95) <=
          5 &&
          safe_quantile(
            eligible_bins$p_repeatability_sd_pct_branch_span,
            0.95
          ) <=
            2 ~ "possibly_sufficient",
        TRUE ~ "poor_repeatability_or_local_dynamics"
      )
    )
  })
  dplyr::bind_rows(branch_rows)
}

create_repeatability_plot <- function(by_bin, summary_df, output_dir) {
  plot_data <- by_bin %>%
    dplyr::filter(sufficient_repeats) %>%
    dplyr::transmute(
      branch,
      voltage_bin_center_v,
      settled_displacement_sd = s_repeatability_sd_pct_branch_span,
      settled_pressure_sd = p_repeatability_sd_pct_branch_span
    ) %>%
    tidyr::pivot_longer(
      cols = c(settled_displacement_sd, settled_pressure_sd),
      names_to = "metric",
      values_to = "repeatability_pct"
    )

  final_metric_data <- summary_df %>%
    dplyr::transmute(
      branch,
      settled_displacement_sd = s_repeatability_sd_p95_pct_branch_span,
      settled_pressure_sd = p_repeatability_sd_p95_pct_branch_span
    ) %>%
    tidyr::pivot_longer(
      cols = c(settled_displacement_sd, settled_pressure_sd),
      names_to = "metric",
      values_to = "final_value_pct"
    ) %>%
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c("settled_displacement_sd", "settled_pressure_sd")
      ),
      label = sprintf(
        "P95 settled SD = %.2f%%",
        final_value_pct
      )
    )

  x_min <- min(plot_data$voltage_bin_center_v, na.rm = TRUE)
  x_max <- max(plot_data$voltage_bin_center_v, na.rm = TRUE)
  x_span <- x_max - x_min
  final_metric_data <- final_metric_data %>%
    dplyr::mutate(
      label_x = ifelse(
        metric == "settled_displacement_sd",
        x_min + 0.02 * x_span,
        x_max - 0.02 * x_span
      ),
      label_hjust = ifelse(metric == "settled_displacement_sd", 0, 1),
      label_vjust = dplyr::case_when(
        metric == "settled_pressure_sd" & branch == "up" ~ 1.5,
        TRUE ~ -0.45
      )
    )

  facet_labels <- c(
    settled_displacement_sd = "Displacement",
    settled_pressure_sd = "Pressure"
  )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = voltage_bin_center_v,
      y = repeatability_pct,
      color = branch
    )
  ) +
    ggplot2::geom_hline(
      data = final_metric_data,
      ggplot2::aes(yintercept = final_value_pct, color = branch),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = final_metric_data,
      ggplot2::aes(
        x = label_x,
        y = final_value_pct,
        label = label,
        color = branch,
        hjust = label_hjust,
        vjust = label_vjust
      ),
      inherit.aes = FALSE,
      size = 3.1,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 1.8,
      alpha = 0.9
    ) +
    ggplot2::facet_grid(
      metric ~ ., 
      scales = "free_y",
      labeller = ggplot2::as_labeller(facet_labels)
    ) +
    ggsci::scale_color_nejm(
      name = "Branch",
      labels = c(down = "Downward", up = "Upward")
    ) +
    ggplot2::labs(
      x = "Voltage-bin center [V]",
      y = "Repeatability: settled SD / branch span [%]"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, -2, 0),
      legend.spacing.x = grid::unit(2, "pt"),
      plot.margin = ggplot2::margin(2, 3, 2, 3),
      strip.text = ggplot2::element_text(face = "plain", size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    filename = file.path(output_dir, "local_static_map_repeatability.png"),
    plot = plot,
    width = 15,
    height = 10,
    units = "cm",
    dpi = 600
  )
}

repo_dir <- normalizePath(file.path(script_dir, ".."))
settings <- parse_settings(list(
  input_csv = file.path(
    "experiment",
    "excitation_experiment_100hz_15min_4s_hold.csv"
  ),
  output_dir = "analysis_outputs",
  voltage_bin_width_v = VOLTAGE_BIN_WIDTH_V,
  terminal_tail_fraction = TERMINAL_TAIL_FRACTION,
  min_repeats_per_bin = MIN_REPEATS_PER_BIN,
  max_terminal_s_sd_mm = MAX_TERMINAL_S_SD_MM,
  max_terminal_s_slope_mm_s = MAX_TERMINAL_S_SLOPE_MM_S,
  max_terminal_p_sd_bar = MAX_TERMINAL_P_SD_BAR,
  max_terminal_p_slope_bar_s = MAX_TERMINAL_P_SLOPE_BAR_S,
  full_hold_count_tolerance = 2L
))

if (
  !is.finite(settings$voltage_bin_width_v) || settings$voltage_bin_width_v <= 0
) {
  stop("voltage_bin_width_v must be positive.")
}

input_path <- if (grepl("^/", settings$input_csv)) {
  settings$input_csv
} else {
  file.path(repo_dir, settings$input_csv)
}
output_dir <- if (grepl("^/", settings$output_dir)) {
  settings$output_dir
} else {
  file.path(repo_dir, settings$output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw_df <- load_excitation(input_path)
ts_s <- infer_sampling_time_s(raw_df)
segmented_df <- segment_holds(raw_df, ts_s)
hold_table <- make_hold_table(segmented_df, ts_s, settings)
by_bin <- make_local_bin_summary(hold_table, settings)
summary_df <- make_overview(
  hold_table,
  by_bin,
  settings,
  settings$input_csv,
  ts_s
)

hold_audit <- hold_table %>%
  dplyr::mutate(
    analysis_name = "local_static_map_repeatability",
    terminal_tail_fraction = settings$terminal_tail_fraction,
    min_repeats_per_bin = settings$min_repeats_per_bin,
    max_terminal_s_sd_mm = settings$max_terminal_s_sd_mm,
    max_terminal_s_slope_mm_s = settings$max_terminal_s_slope_mm_s,
    max_terminal_p_sd_bar = settings$max_terminal_p_sd_bar,
    max_terminal_p_slope_bar_s = settings$max_terminal_p_slope_bar_s
  ) %>%
  dplyr::select(
    analysis_name,
    hold_id,
    branch,
    n_samples,
    full_modal_n,
    is_full_hold,
    hold_duration_s,
    terminal_tail_n_samples,
    terminal_tail_start_s,
    u_previous_v,
    u_cmd_v,
    delta_u_v,
    abs_delta_u_v,
    voltage_bin_width_v,
    voltage_bin_lower_v,
    voltage_bin_upper_v,
    voltage_bin_center_v,
    voltage_bin,
    s_start_mm,
    s_terminal_mm,
    p_start_bar,
    p_terminal_bar,
    s_terminal_sd_mm,
    p_terminal_sd_bar,
    s_terminal_slope_mm_s,
    p_terminal_slope_bar_s,
    terminal_stable,
    eligible_for_repeatability,
    exclusion_reason,
    terminal_tail_fraction,
    min_repeats_per_bin,
    max_terminal_s_sd_mm,
    max_terminal_s_slope_mm_s,
    max_terminal_p_sd_bar,
    max_terminal_p_slope_bar_s
  )

readr::write_csv(
  hold_audit,
  file.path(output_dir, "local_static_map_repeatability_hold_points.csv")
)
readr::write_csv(
  by_bin,
  file.path(output_dir, "local_static_map_repeatability_by_bin.csv")
)
readr::write_csv(
  summary_df,
  file.path(output_dir, "local_static_map_repeatability_summary.csv")
)
create_repeatability_plot(by_bin, summary_df, output_dir)

cat(sprintf(
  "Local static-map repeatability complete: %d holds, %d stable directional holds, %.3f V bins, outputs in %s\n",
  nrow(hold_table),
  sum(hold_table$eligible_for_repeatability),
  settings$voltage_bin_width_v,
  output_dir
))
print(summary_df)
