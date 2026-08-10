#!/usr/bin/env Rscript

# Resolve the script directory for both `Rscript` and sourced execution.
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

# Parse simple `--name=value` overrides using the default value type.
parse_settings <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    return(defaults)
  }

  settings <- defaults

  for (arg in args) {
    arg_clean <- sub("^--", "", arg)
    pieces <- strsplit(arg_clean, "=", fixed = TRUE)[[1]]

    if (length(pieces) < 2) {
      stop(sprintf("Arguments must use --name=value format: %s", arg))
    }

    key <- pieces[1]
    value <- paste(pieces[-1], collapse = "=")

    if (!key %in% names(settings)) {
      stop(sprintf("Unknown setting override: %s", key))
    }

    default_value <- settings[[key]]

    if (is.integer(default_value)) {
      settings[[key]] <- as.integer(value)
    } else if (is.numeric(default_value)) {
      settings[[key]] <- as.numeric(value)
    } else if (is.logical(default_value)) {
      settings[[key]] <- tolower(value) %in% c("true", "1", "yes")
    } else {
      settings[[key]] <- value
    }
  }

  settings
}

# Load the global conservative window from the existing analysis output.
read_global_control_window <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing control-window file: %s", path))
  }

  window_df <- read_csv(path, show_col_types = FALSE)
  global_row <- window_df %>%
    filter(level == "global") %>%
    slice(1)

  if (nrow(global_row) == 0) {
    stop("No global control window found in recommended_control_window.csv.")
  }

  global_row
}

# Read the calibrated static-analysis summaries used to plan excitation.
read_calibrated_analysis_outputs <- function(output_dir) {
  control_window_path <- file.path(output_dir, "recommended_control_window.csv")

  control_window_df <- read_csv(control_window_path, show_col_types = FALSE)
  control_window <- control_window_df %>%
    filter(level == "global") %>%
    slice(1)
  run_rows <- control_window_df %>%
    filter(level == "run")

  if (nrow(run_rows) == 0) {
    stop("No run-level rows found in recommended_control_window.csv.")
  }

  if (!all(grepl("_p_voltage_bar\\.csv$", run_rows$file))) {
    stop(
      paste(
        "recommended_control_window.csv does not appear to come from calibrated",
        "*_p_voltage_bar.csv experiments. Re-run analyze_static_threshold_hysteresis.R first."
      )
    )
  }

  list(
    control_window = control_window
  )
}

# Read the global effective settling-time recommendation if available.
read_effective_settling_time_summary <- function(output_dir) {
  path <- file.path(output_dir, "effective_settling_time_summary.csv")

  if (!file.exists(path)) {
    return(NULL)
  }

  summary_df <- read_csv(path, show_col_types = FALSE)
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)

  if (
    nrow(global_row) != 1 ||
      !is.finite(global_row$recommended_t_settle_ms[[1]]) ||
      !is.finite(global_row$recommended_t_hold_s[[1]]) ||
      !nzchar(global_row$recommendation_basis[[1]])
  ) {
    stop(
      paste(
        "effective_settling_time_summary.csv exists but does not contain one valid global row.",
        "Re-run analyze_effective_settling_time.R."
      )
    )
  }

  global_row
}

# Read the global maximum sampling-time analysis if available.
read_max_sampling_time_summary <- function(output_dir) {
  path <- file.path(output_dir, "max_sampling_time_summary.csv")

  if (!file.exists(path)) {
    return(NULL)
  }

  summary_df <- read_csv(path, show_col_types = FALSE)
  global_row <- summary_df %>%
    filter(level == "global") %>%
    slice(1)

  if (
    nrow(global_row) != 1 ||
      !is.finite(global_row$recommended_max_sampling_time_ms[[1]]) ||
      !is.finite(global_row$max_supported_sampling_rate_hz[[1]])
  ) {
    stop(
      paste(
        "max_sampling_time_summary.csv exists but does not contain one valid global row.",
        "Re-run analyze_max_sampling_time.R."
      )
    )
  }

  global_row
}

# Resolve the commanded sample rate with a conservative safety margin.
resolve_sampling_rate_inputs <- function(
  output_dir,
  fs_override_hz,
  safety_factor = 8,
  rounding_grid_hz = 50,
  fallback_hz = 100
) {
  if (is.finite(fs_override_hz)) {
    fs_hz <- as.integer(round(fs_override_hz))
    if (fs_hz <= 0) {
      stop("fs_hz override must be positive.")
    }

    return(list(
      fs_hz = fs_hz,
      fs_hz_source = "user_override",
      sampling_rate_margin_factor = NA_real_,
      observed_max_sampling_time_ms = NA_real_,
      max_supported_sampling_rate_hz = NA_real_,
      recommended_max_sampling_time_ms = NA_real_,
      recommendation_basis = "manual_fs_hz_override"
    ))
  }

  max_sampling_row <- read_max_sampling_time_summary(output_dir)
  if (!is.null(max_sampling_row)) {
    max_rate_hz <- as.numeric(max_sampling_row$max_supported_sampling_rate_hz[[1]])
    target_rate_hz <- max_rate_hz / safety_factor
    fs_hz <- max(
      rounding_grid_hz,
      as.integer(round(target_rate_hz / rounding_grid_hz) * rounding_grid_hz)
    )

    return(list(
      fs_hz = fs_hz,
      fs_hz_source = "max_sampling_time_analysis",
      sampling_rate_margin_factor = safety_factor,
      observed_max_sampling_time_ms = as.numeric(max_sampling_row$observed_max_sampling_time_ms[[1]]),
      max_supported_sampling_rate_hz = max_rate_hz,
      recommended_max_sampling_time_ms = as.numeric(max_sampling_row$recommended_max_sampling_time_ms[[1]]),
      recommendation_basis = sprintf(
        "rounded_to_%dHz_grid_from_measured_max_rate_divided_by_%g",
        rounding_grid_hz,
        safety_factor
      )
    ))
  }

  list(
    fs_hz = as.integer(fallback_hz),
    fs_hz_source = "temporary_fallback_default",
    sampling_rate_margin_factor = NA_real_,
    observed_max_sampling_time_ms = NA_real_,
    max_supported_sampling_rate_hz = NA_real_,
    recommended_max_sampling_time_ms = NA_real_,
    recommendation_basis = "temporary sampling-rate fallback because max sampling-time summary not available"
  )
}

# Resolve the settling-time and hold-time inputs used for excitation planning.
resolve_hold_time_inputs <- function(output_dir, t_settle_override_s, hold_factor_override) {
  if (is.finite(t_settle_override_s)) {
    t_settle_slow_s <- as.numeric(t_settle_override_s)
    hold_factor_applied <- if (is.finite(hold_factor_override)) {
      as.numeric(hold_factor_override)
    } else {
      NA_real_
    }
    T_hold_s <- if (is.finite(hold_factor_applied)) {
      hold_factor_applied * t_settle_slow_s
    } else {
      t_settle_slow_s
    }

    return(list(
      t_settle_slow_s = t_settle_slow_s,
      recommended_t_settle_ms = 1000 * t_settle_slow_s,
      recommended_t_settle_ms_source = "user_override",
      recommended_t_hold_s = T_hold_s,
      T_hold_s = T_hold_s,
      T_hold_source = if (is.finite(hold_factor_applied)) {
        "user_override_with_hold_factor"
      } else {
        "user_override_direct"
      },
      hold_factor_applied = hold_factor_applied,
      recommendation_basis = if (is.finite(hold_factor_applied)) {
        "manual_override_with_explicit_hold_factor"
      } else {
        "manual_override_direct"
      }
    ))
  }

  settling_row <- read_effective_settling_time_summary(output_dir)
  if (!is.null(settling_row)) {
    return(list(
      t_settle_slow_s = as.numeric(settling_row$recommended_t_settle_ms[[1]]) / 1000,
      recommended_t_settle_ms = as.numeric(settling_row$recommended_t_settle_ms[[1]]),
      recommended_t_settle_ms_source = "effective_settling_time_analysis",
      recommended_t_hold_s = as.numeric(settling_row$recommended_t_hold_s[[1]]),
      T_hold_s = as.numeric(settling_row$recommended_t_hold_s[[1]]),
      T_hold_source = "effective_settling_time_analysis",
      hold_factor_applied = NA_real_,
      recommendation_basis = as.character(settling_row$recommendation_basis[[1]])
    ))
  }

  list(
    t_settle_slow_s = 1.0,
    recommended_t_settle_ms = 1000,
    recommended_t_settle_ms_source = "temporary_fallback_default",
    recommended_t_hold_s = 1.0,
    T_hold_s = 1.0,
    T_hold_source = "temporary_fallback_default",
    hold_factor_applied = NA_real_,
    recommendation_basis = "temporary fallback default because effective settling-time summary not available"
  )
}

# Allocate integer counts that stay close to the requested proportions.
allocate_counts <- function(total_n, target_props) {
  raw_counts <- total_n * target_props
  counts <- floor(raw_counts)
  remainder <- total_n - sum(counts)

  if (remainder > 0) {
    order_idx <- order(raw_counts - counts, decreasing = TRUE)
    counts[order_idx[seq_len(remainder)]] <- counts[order_idx[seq_len(
      remainder
    )]] +
      1L
  }

  stats::setNames(as.integer(counts), names(target_props))
}

# Split each step-size bin approximately evenly across upward and downward steps.
allocate_direction_bin_counts <- function(bin_counts_target) {
  direction_bin_df <- tibble(
    step_bin = names(bin_counts_target),
    total_count = as.integer(bin_counts_target)
  ) %>%
    arrange(desc(total_count), step_bin) %>%
    mutate(
      up_count = total_count %/% 2L,
      down_count = total_count %/% 2L
    )

  odd_idx <- which(direction_bin_df$total_count %% 2L == 1L)

  if (length(odd_idx) > 0) {
    total_up <- sum(direction_bin_df$up_count)
    total_down <- sum(direction_bin_df$down_count)

    for (idx in odd_idx) {
      if (total_up <= total_down) {
        direction_bin_df$up_count[idx] <- direction_bin_df$up_count[idx] + 1L
        total_up <- total_up + 1L
      } else {
        direction_bin_df$down_count[idx] <- direction_bin_df$down_count[idx] + 1L
        total_down <- total_down + 1L
      }
    }
  }

  direction_bin_df %>%
    dplyr::select(step_bin, up_count, down_count) %>%
    pivot_longer(
      cols = c(up_count, down_count),
      names_to = "direction",
      values_to = "count"
    ) %>%
    mutate(direction = recode(direction, up_count = "up", down_count = "down"))
}

# Derive transition count from hold time and the requested duration budget.
plan_transition_count <- function(
  T_hold_s,
  target_duration_min,
  hard_max_duration_min,
  min_transitions,
  max_transitions
) {
  target_duration_s <- target_duration_min * 60
  hard_max_duration_s <- hard_max_duration_min * 60

  target_transitions <- max(1L, floor(target_duration_s / T_hold_s) - 1L)
  feasible_transitions_max <- max(
    1L,
    floor(hard_max_duration_s / T_hold_s) - 1L
  )

  n_transitions <- min(max_transitions, target_transitions)
  used_hard_cap <- FALSE
  planning_warning <- NA_character_

  if (n_transitions < min_transitions) {
    used_hard_cap <- TRUE
    n_transitions <- min(max_transitions, feasible_transitions_max)

    if (n_transitions < min_transitions) {
      planning_warning <- paste(
        "Minimum transition target could not be reached within the hard duration cap.",
        "Using the maximum feasible number of transitions instead."
      )
    }
  }

  n_holds <- n_transitions + 1L
  total_duration_s <- n_holds * T_hold_s

  list(
    n_transitions = as.integer(n_transitions),
    n_holds = as.integer(n_holds),
    total_duration_s = total_duration_s,
    used_hard_cap = used_hard_cap,
    planning_warning = planning_warning
  )
}

# Generate one bounded random excitation sequence that matches target quotas.
generate_transition_sequence <- function(
  n_transitions,
  u_low_exc,
  u_high_exc,
  seed,
  direction_bin_counts_target,
  step_bin_defs
) {
  u_span_exc <- u_high_exc - u_low_exc

  for (attempt in seq_len(200)) {
    set.seed(seed + attempt - 1L)

    remaining_targets <- direction_bin_counts_target %>%
      left_join(step_bin_defs, by = "step_bin") %>%
      mutate(
        total_count = count,
        remaining_strata = purrr::map(count, function(n) {
          if (n <= 0) {
            integer(0)
          } else {
            seq_len(n)
          }
        })
      )
    current_u <- mean(c(u_low_exc, u_high_exc))
    transition_rows <- vector("list", n_transitions)
    success <- TRUE

    for (i in seq_len(n_transitions)) {
      available_span <- c(
        up = u_high_exc - current_u,
        down = current_u - u_low_exc
      )

      candidate_df <- remaining_targets %>%
        mutate(candidate_id = row_number()) %>%
        filter(count > 0) %>%
        mutate(
          available_v = available_span[direction],
          available_frac = available_v / u_span_exc,
          feasible_strata = pmap(
            list(lower_frac, upper_frac, total_count, remaining_strata, available_frac),
            function(lower_frac, upper_frac, total_count, remaining_strata, available_frac) {
              if (length(remaining_strata) == 0 || !is.finite(available_frac)) {
                return(integer(0))
              }

              stratum_edges <- seq(lower_frac, upper_frac, length.out = total_count + 1L)
              remaining_strata[
                stratum_edges[remaining_strata] <= available_frac + 1e-9
              ]
            }
          ),
          feasible_count = lengths(feasible_strata)
        ) %>%
        filter(feasible_count > 0) %>%
        mutate(
          weight = feasible_count
        )

      if (nrow(candidate_df) == 0) {
        success <- FALSE
        break
      }

      pick_idx <- sample(
        seq_len(nrow(candidate_df)),
        size = 1,
        prob = candidate_df$weight
      )
      pick <- candidate_df[pick_idx, ]
      chosen_stratum <- sample(pick$feasible_strata[[1]], size = 1)
      stratum_edges <- seq(
        pick$lower_frac,
        pick$upper_frac,
        length.out = pick$total_count + 1L
      )
      stratum_low <- stratum_edges[chosen_stratum]
      stratum_high <- min(stratum_edges[chosen_stratum + 1L], pick$available_frac)

      if (!is.finite(stratum_high) || stratum_high < stratum_low) {
        success <- FALSE
        break
      }

      step_frac <- runif(1, min = stratum_low, max = stratum_high)
      step_sign <- ifelse(pick$direction == "up", 1, -1)
      next_u <- current_u + step_sign * step_frac * u_span_exc
      next_u <- min(max(next_u, u_low_exc), u_high_exc)
      delta_u_v <- next_u - current_u
      delta_u_pct_span <- 100 * abs(delta_u_v) / u_span_exc

      transition_rows[[i]] <- tibble(
        trans_id = i,
        u_prev_v = current_u,
        u_cmd_v = next_u,
        delta_u_v = delta_u_v,
        delta_u_pct_span = delta_u_pct_span,
        direction = pick$direction,
        step_bin = pick$step_bin
      )

      remaining_targets$count[pick$candidate_id] <- remaining_targets$count[pick$candidate_id] - 1L
      remaining_targets$remaining_strata[[pick$candidate_id]] <-
        setdiff(remaining_targets$remaining_strata[[pick$candidate_id]], chosen_stratum)
      current_u <- next_u
    }

    if (success && all(remaining_targets$count == 0L)) {
      return(bind_rows(transition_rows))
    }
  }

  stop(
    "Could not generate a bounded excitation sequence that satisfies the requested quotas."
  )
}

# Convert transitions into a hold-based schedule.
build_excitation_schedule <- function(
  transition_df,
  T_hold_s,
  u_start_v,
  u_span_exc
) {
  initial_hold <- tibble(
    hold_id = 1L,
    t_start_s = 0,
    t_end_s = T_hold_s,
    u_cmd_v = u_start_v,
    delta_u_v = 0,
    delta_u_pct_span = 0,
    direction = "start",
    step_bin = "start"
  )

  transition_holds <- transition_df %>%
    transmute(
      hold_id = trans_id + 1L,
      t_start_s = trans_id * T_hold_s,
      t_end_s = (trans_id + 1L) * T_hold_s,
      u_cmd_v = u_cmd_v,
      delta_u_v = delta_u_v,
      delta_u_pct_span = delta_u_pct_span,
      direction = direction,
      step_bin = step_bin
    )

  bind_rows(initial_hold, transition_holds)
}

# Expand the hold schedule to one row per sample at the commanded sampling rate.
expand_schedule_to_samples <- function(schedule_df, fs_hz) {
  samples_per_hold <- max(
    1L,
    as.integer(round(
      (schedule_df$t_end_s[[1]] - schedule_df$t_start_s[[1]]) * fs_hz
    ))
  )

  sample_df <- map_dfr(seq_len(nrow(schedule_df)), function(i) {
    hold_row <- schedule_df[i, ]
    tibble(
      time_s = hold_row$t_start_s + (seq_len(samples_per_hold) - 1L) / fs_hz,
      hold_id = hold_row$hold_id,
      u_cmd_v = hold_row$u_cmd_v
    )
  }) %>%
    mutate(sample_id = row_number()) %>%
    relocate(sample_id)

  sample_df
}

# Summarize actual signal properties and constraint checks.
create_excitation_summary <- function(
  control_window,
  hold_time_resolution,
  sampling_rate_resolution,
  schedule_df,
  sample_df,
  settings,
  T_hold_s,
  margin_u,
  u_low_exc,
  u_high_exc,
  step_bin_defs,
  planning_info
) {
  transition_df <- schedule_df %>%
    filter(hold_id > 1)

  u_span_exc <- u_high_exc - u_low_exc
  up_count <- sum(transition_df$direction == "up")
  down_count <- sum(transition_df$direction == "down")
  up_down_imbalance_pct <- 100 *
    abs(up_count - down_count) /
    max(1, up_count + down_count)

  bin_summary <- transition_df %>%
    count(step_bin, name = "count") %>%
    right_join(step_bin_defs, by = "step_bin") %>%
    mutate(
      count = replace_na(count, 0L),
      prop = count / max(1, nrow(transition_df)),
      prop_ok = prop >= target_low_prop & prop <= target_high_prop
    )

  bounds_ok <- all(
    sample_df$u_cmd_v >= u_low_exc - 1e-9 &
      sample_df$u_cmd_v <= u_high_exc + 1e-9
  )
  balance_ok <- up_down_imbalance_pct <= 10
  near_zero_ok <- all(abs(transition_df$delta_u_v) >= 0.05 * u_span_exc - 1e-9)
  all_constraints_ok <- bounds_ok &&
    balance_ok &&
    near_zero_ok &&
    all(bin_summary$prop_ok)

  summary_base <- tibble(
    seed = settings$seed,
    fs_hz = sampling_rate_resolution$fs_hz,
    fs_hz_source = sampling_rate_resolution$fs_hz_source,
    sampling_rate_margin_factor = sampling_rate_resolution$sampling_rate_margin_factor,
    observed_max_sampling_time_ms = sampling_rate_resolution$observed_max_sampling_time_ms,
    recommended_max_sampling_time_ms = sampling_rate_resolution$recommended_max_sampling_time_ms,
    max_supported_sampling_rate_hz = sampling_rate_resolution$max_supported_sampling_rate_hz,
    t_settle_slow_s = hold_time_resolution$t_settle_slow_s,
    recommended_t_settle_ms = hold_time_resolution$recommended_t_settle_ms,
    recommended_t_settle_ms_source = hold_time_resolution$recommended_t_settle_ms_source,
    recommended_t_hold_s = hold_time_resolution$recommended_t_hold_s,
    hold_factor = settings$hold_factor,
    hold_factor_applied = hold_time_resolution$hold_factor_applied,
    T_hold_source = hold_time_resolution$T_hold_source,
    recommendation_basis = hold_time_resolution$recommendation_basis,
    T_hold_s = T_hold_s,
    target_duration_min = settings$target_duration_min,
    hard_max_duration_min = settings$hard_max_duration_min,
    min_transitions = settings$min_transitions,
    max_transitions = settings$max_transitions,
    planned_transitions = nrow(transition_df),
    planned_holds = nrow(schedule_df),
    total_duration_min = max(schedule_df$t_end_s) / 60,
    used_hard_cap = planning_info$used_hard_cap,
    planning_warning = planning_info$planning_warning %||% NA_character_,
    control_window_low_v = control_window$window_low_v[[1]],
    control_window_high_v = control_window$window_high_v[[1]],
    control_window_width_v = control_window$window_width_v[[1]],
    source_runs_used = control_window$runs_used[[1]],
    source_any_direction_unstable = control_window$any_direction_unstable[[1]],
    source_any_stroke_endpoint_flag = control_window$any_stroke_endpoint_flag[[
      1
    ]],
    margin_fraction = settings$margin_fraction,
    margin_u = margin_u,
    excitation_low_v = u_low_exc,
    excitation_high_v = u_high_exc,
    excitation_span_v = u_span_exc,
    up_count = up_count,
    down_count = down_count,
    up_down_imbalance_pct = up_down_imbalance_pct,
    bounds_ok = bounds_ok,
    balance_ok = balance_ok,
    near_zero_ok = near_zero_ok
  )

  bin_wide <- bin_summary %>%
    transmute(
      step_bin,
      count,
      prop,
      prop_ok
    ) %>%
    tidyr::pivot_wider(
      names_from = step_bin,
      values_from = c(count, prop, prop_ok),
      names_sep = "_"
    )

  bind_cols(
    summary_base,
    bin_wide,
    tibble(all_constraints_ok = all_constraints_ok)
  )
}

# Build regular axis breaks over the plotted data range.
regular_breaks <- function(step, anchor = 0) {
  force(step)
  force(anchor)

  function(x) {
    finite_x <- x[is.finite(x)]
    if (length(finite_x) == 0) {
      return(numeric(0))
    }

    lower <- anchor + floor((min(finite_x) - anchor) / step) * step
    upper <- anchor + ceiling((max(finite_x) - anchor) / step) * step
    seq(lower, upper, by = step)
  }
}

# Plot the generated command signal over time.
create_voltage_time_plot <- function(schedule_df, analysis_plot_theme) {
  ggplot(schedule_df, aes(x = t_start_s / 60, y = u_cmd_v)) +
    geom_step(color = "#1F1F1F", linewidth = 0.28) +
    scale_x_continuous(breaks = regular_breaks(1)) +
    scale_y_continuous(breaks = regular_breaks(1)) +
    labs(
      x = "Time [min]",
      y = "Command voltage [V]"
    ) +
    analysis_plot_theme +
    theme(legend.position = "none")
}

# Plot the sampled command-voltage histogram.
create_voltage_histogram_plot <- function(sample_df, analysis_plot_theme) {
  ggplot(sample_df, aes(x = u_cmd_v)) +
    geom_histogram(
      bins = 40,
      fill = "#3C5488",
      color = "white",
      linewidth = 0.1
    ) +
    scale_x_continuous(breaks = regular_breaks(0.5)) +
    scale_y_continuous(breaks = regular_breaks(1000)) +
    labs(
      x = "Command voltage [V]",
      y = "Count"
    ) +
    analysis_plot_theme +
    theme(legend.position = "none")
}

# Plot the signed step-change histogram.
create_signed_step_histogram_plot <- function(
  schedule_df,
  analysis_plot_theme
) {
  schedule_df %>%
    filter(hold_id > 1) %>%
    ggplot(aes(x = delta_u_v, fill = direction)) +
    geom_histogram(bins = 40, alpha = 0.85, color = "white", linewidth = 0.1) +
    scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    scale_x_continuous(breaks = regular_breaks(0.5)) +
    scale_y_continuous(breaks = regular_breaks(10)) +
    labs(
      x = expression(Delta * "u [V]"),
      y = "Count",
      fill = "Direction"
    ) +
    analysis_plot_theme
}

# Plot actual step-bin proportions against their target ranges.
create_step_bin_plot <- function(
  schedule_df,
  step_bin_defs,
  analysis_plot_theme
) {
  actual_df <- schedule_df %>%
    filter(hold_id > 1) %>%
    count(step_bin, name = "count") %>%
    right_join(step_bin_defs, by = "step_bin") %>%
    mutate(
      count = replace_na(count, 0L),
      prop = count / max(1, sum(count)),
      step_bin = factor(
        step_bin,
        levels = c("xsmall", "small", "medium", "large", "xlarge")
      ),
      step_bin_label = factor(
        step_bin_label,
        levels = step_bin_defs %>%
          arrange(match(
            step_bin,
            c("xsmall", "small", "medium", "large", "xlarge")
          )) %>%
          pull(step_bin_label)
      )
    )

  ggplot(actual_df, aes(x = step_bin_label, y = prop)) +
    geom_linerange(
      aes(ymin = target_low_prop, ymax = target_high_prop),
      color = "grey40",
      linewidth = 1.5
    ) +
    geom_point(size = 2.2, color = "#B24745") +
    scale_y_continuous(
      breaks = regular_breaks(0.1),
      labels = function(x) sprintf("%.0f%%", 100 * x)
    ) +
    labs(
      x = expression(
        "Step-size bin and " *
          "|" *
          Delta *
          "u| / " *
          Delta *
          "u"[exc] *
          " bounds"
      ),
      y = "Proportion of transitions"
    ) +
    analysis_plot_theme +
    theme(legend.position = "none")
}

# Plot cumulative up/down transition counts through the schedule.
create_direction_balance_plot <- function(schedule_df, analysis_plot_theme) {
  balance_df <- schedule_df %>%
    filter(hold_id > 1) %>%
    mutate(
      trans_id = row_number(),
      cum_up = cumsum(direction == "up"),
      cum_down = cumsum(direction == "down")
    ) %>%
    dplyr::select(trans_id, cum_up, cum_down) %>%
    pivot_longer(
      cols = c(cum_up, cum_down),
      names_to = "series",
      values_to = "count"
    ) %>%
    mutate(
      series = recode(
        series,
        cum_up = "Upward transitions",
        cum_down = "Downward transitions"
      )
    )

  ggplot(balance_df, aes(x = trans_id, y = count, color = series)) +
    geom_line(linewidth = 0.4) +
    scale_color_nejm(name = "Series") +
    scale_x_continuous(breaks = regular_breaks(50)) +
    scale_y_continuous(breaks = regular_breaks(100)) +
    labs(
      x = "Transition index",
      y = "Cumulative count",
      color = "Series"
    ) +
    analysis_plot_theme
}

# Save all excitation plots to the analysis output directory.
save_excitation_plots <- function(
  output_dir,
  voltage_time_plot,
  voltage_histogram_plot,
  signed_step_histogram_plot,
  step_bin_plot,
  direction_balance_plot
) {
  ggsave(
    filename = file.path(output_dir, "excitation_voltage_time.png"),
    plot = voltage_time_plot,
    width = 15,
    height = 5,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "excitation_voltage_histogram.png"),
    plot = voltage_histogram_plot,
    width = 15,
    height = 3.75,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "excitation_signed_step_histogram.png"),
    plot = signed_step_histogram_plot,
    width = 15,
    height = 3.75,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "excitation_step_bin_proportions.png"),
    plot = step_bin_plot,
    width = 15,
    height = 5,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "excitation_direction_balance.png"),
    plot = direction_balance_plot,
    width = 15,
    height = 3.75,
    units = "cm",
    dpi = 600
  )
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

defaults <- list(
  t_settle_slow_s = NA_real_,
  hold_factor = NA_real_,
  fs_hz = NA_real_,
  target_duration_min = 15.0,
  hard_max_duration_min = 25.0,
  min_transitions = 300L,
  max_transitions = 1000L,
  margin_fraction = 0.05,
  max_excitation_high_v = 7.2,
  seed = 12345L
)

settings <- parse_settings(defaults)

analysis_outputs <- read_calibrated_analysis_outputs(
  analysis_context$output_dir
)
control_window <- analysis_outputs$control_window
hold_time_resolution <- resolve_hold_time_inputs(
  output_dir = analysis_context$output_dir,
  t_settle_override_s = settings$t_settle_slow_s,
  hold_factor_override = settings$hold_factor
)
sampling_rate_resolution <- resolve_sampling_rate_inputs(
  output_dir = analysis_context$output_dir,
  fs_override_hz = settings$fs_hz
)

Delta_u_ctrl <- control_window$window_high_v[[1]] -
  control_window$window_low_v[[1]]
margin_u <- settings$margin_fraction * Delta_u_ctrl
u_low_exc <- control_window$window_low_v[[1]] + margin_u
u_high_exc <- min(
  control_window$window_high_v[[1]] - margin_u,
  settings$max_excitation_high_v
)

if (
  !is.finite(u_low_exc) || !is.finite(u_high_exc) || u_low_exc >= u_high_exc
) {
  stop("Excitation bounds are invalid after applying the inner voltage margin.")
}

T_hold_s <- hold_time_resolution$T_hold_s
planning_info <- plan_transition_count(
  T_hold_s = T_hold_s,
  target_duration_min = settings$target_duration_min,
  hard_max_duration_min = settings$hard_max_duration_min,
  min_transitions = settings$min_transitions,
  max_transitions = settings$max_transitions
)

step_bin_defs <- tibble(
  step_bin = c("xsmall", "small", "medium", "large", "xlarge"),
  lower_frac = c(0.05, 0.10, 0.25, 0.40, 0.60),
  upper_frac = c(0.10, 0.25, 0.40, 0.60, 0.85),
  target_prop = c(0.10, 0.50, 0.25, 0.12, 0.03),
  target_low_prop = c(0.05, 0.45, 0.20, 0.10, 0.00),
  target_high_prop = c(0.15, 0.55, 0.30, 0.15, 0.05)
) %>%
  mutate(
    step_bin_label = sprintf(
      "%s\n%.0f-%.0f%% of span",
      recode(
        step_bin,
        xsmall = "Very small",
        .default = str_to_title(step_bin)
      ),
      100 * lower_frac,
      100 * upper_frac
    )
  )

bin_counts_target <- allocate_counts(
  planning_info$n_transitions,
  setNames(step_bin_defs$target_prop, step_bin_defs$step_bin)
)
direction_bin_counts_target <- allocate_direction_bin_counts(bin_counts_target)

u_start_v <- mean(c(u_low_exc, u_high_exc))
transition_df <- generate_transition_sequence(
  n_transitions = planning_info$n_transitions,
  u_low_exc = u_low_exc,
  u_high_exc = u_high_exc,
  seed = settings$seed,
  direction_bin_counts_target = direction_bin_counts_target,
  step_bin_defs = step_bin_defs
)

schedule_df <- build_excitation_schedule(
  transition_df = transition_df,
  T_hold_s = T_hold_s,
  u_start_v = u_start_v,
  u_span_exc = u_high_exc - u_low_exc
)

sample_df <- expand_schedule_to_samples(
  schedule_df = schedule_df,
  fs_hz = sampling_rate_resolution$fs_hz
)

# Reuse the exact same transitions with a fixed 4 s hold for an extended-hold experiment.
T_hold_4s_s <- 4.0
schedule_4s_hold_df <- build_excitation_schedule(
  transition_df = transition_df,
  T_hold_s = T_hold_4s_s,
  u_start_v = u_start_v,
  u_span_exc = u_high_exc - u_low_exc
)
sample_4s_hold_df <- expand_schedule_to_samples(
  schedule_df = schedule_4s_hold_df,
  fs_hz = sampling_rate_resolution$fs_hz
)

hold_time_resolution_4s <- hold_time_resolution
hold_time_resolution_4s$recommended_t_hold_s <- T_hold_4s_s
hold_time_resolution_4s$T_hold_s <- T_hold_4s_s
hold_time_resolution_4s$T_hold_source <- "fixed_4s_hold_override"
hold_time_resolution_4s$recommendation_basis <- paste(
  "same transition sequence as the standard excitation with a fixed 4 s hold"
)
settings_4s_hold <- settings
settings_4s_hold$target_duration_min <- max(schedule_4s_hold_df$t_end_s) / 60
settings_4s_hold$hard_max_duration_min <- settings_4s_hold$target_duration_min
planning_info_4s_hold <- planning_info
planning_info_4s_hold$used_hard_cap <- FALSE
planning_info_4s_hold$planning_warning <- NA_character_

summary_df <- create_excitation_summary(
  control_window = control_window,
  hold_time_resolution = hold_time_resolution,
  sampling_rate_resolution = sampling_rate_resolution,
  schedule_df = schedule_df,
  sample_df = sample_df,
  settings = settings,
  T_hold_s = T_hold_s,
  margin_u = margin_u,
  u_low_exc = u_low_exc,
  u_high_exc = u_high_exc,
  step_bin_defs = step_bin_defs,
  planning_info = planning_info
)

summary_4s_hold_df <- create_excitation_summary(
  control_window = control_window,
  hold_time_resolution = hold_time_resolution_4s,
  sampling_rate_resolution = sampling_rate_resolution,
  schedule_df = schedule_4s_hold_df,
  sample_df = sample_4s_hold_df,
  settings = settings_4s_hold,
  T_hold_s = T_hold_4s_s,
  margin_u = margin_u,
  u_low_exc = u_low_exc,
  u_high_exc = u_high_exc,
  step_bin_defs = step_bin_defs,
  planning_info = planning_info_4s_hold
)

write_csv(
  schedule_df,
  file.path(analysis_context$output_dir, "excitation_schedule.csv")
)
write_csv(
  sample_df,
  file.path(analysis_context$output_dir, "excitation_samples.csv")
)
write_csv(
  summary_df,
  file.path(analysis_context$output_dir, "excitation_signal_summary.csv")
)
write_csv(
  schedule_4s_hold_df,
  file.path(analysis_context$output_dir, "excitation_schedule_4s_hold.csv")
)
write_csv(
  sample_4s_hold_df,
  file.path(analysis_context$output_dir, "excitation_samples_4s_hold.csv")
)
write_csv(
  summary_4s_hold_df,
  file.path(analysis_context$output_dir, "excitation_signal_summary_4s_hold.csv")
)

voltage_time_plot <- create_voltage_time_plot(schedule_df, analysis_plot_theme)
voltage_histogram_plot <- create_voltage_histogram_plot(
  sample_df,
  analysis_plot_theme
)
signed_step_histogram_plot <- create_signed_step_histogram_plot(
  schedule_df,
  analysis_plot_theme
)
step_bin_plot <- create_step_bin_plot(
  schedule_df,
  step_bin_defs,
  analysis_plot_theme
)
direction_balance_plot <- create_direction_balance_plot(
  schedule_df,
  analysis_plot_theme
)

save_excitation_plots(
  output_dir = analysis_context$output_dir,
  voltage_time_plot = voltage_time_plot,
  voltage_histogram_plot = voltage_histogram_plot,
  signed_step_histogram_plot = signed_step_histogram_plot,
  step_bin_plot = step_bin_plot,
  direction_balance_plot = direction_balance_plot
)

cat("\nExcitation signal summary\n")
print(summary_df)
cat(sprintf(
  "\nSaved excitation CSVs, 4 s hold CSVs, and plots to: %s\n",
  analysis_context$output_dir
))
