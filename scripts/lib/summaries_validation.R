# Read lower and upper usable-window bounds from different summary table shapes.
get_window_bounds <- function(cutoffs) {
  lower_bound <- NA_real_
  upper_bound <- NA_real_

  if ("window_low_v" %in% names(cutoffs)) {
    lower_bound <- cutoffs$window_low_v[[1]]
  } else if ("overlap_low_v" %in% names(cutoffs)) {
    lower_bound <- cutoffs$overlap_low_v[[1]]
  } else if ("cutoff_low_v" %in% names(cutoffs)) {
    lower_bound <- cutoffs$cutoff_low_v[[1]]
  }

  if ("window_high_v" %in% names(cutoffs)) {
    upper_bound <- cutoffs$window_high_v[[1]]
  } else if ("overlap_high_v" %in% names(cutoffs)) {
    upper_bound <- cutoffs$overlap_high_v[[1]]
  } else if ("cutoff_high_v" %in% names(cutoffs)) {
    upper_bound <- cutoffs$cutoff_high_v[[1]]
  }

  list(lower = lower_bound, upper = upper_bound)
}

# Summarize cycle-to-cycle repeatability inside the global control window.
assess_settling <- function(df, cutoffs) {
  window_bounds <- get_window_bounds(cutoffs)
  lower_bound <- window_bounds$lower
  upper_bound <- window_bounds$upper

  df_use <- df %>%
    filter(u >= lower_bound, u <= upper_bound)

  if (
    nrow(df_use) < 20 ||
      anyNA(lower_bound) ||
      anyNA(upper_bound)
  ) {
    return(tibble(
      n_usable = nrow(df_use),
      window_low_v = lower_bound,
      window_high_v = upper_bound,
      p_span_global = NA_real_,
      s_span_global = NA_real_,
      p_sd_bar_p95 = NA_real_,
      s_sd_mm_p95 = NA_real_,
      p_mad_bar_median = NA_real_,
      s_mad_mm_median = NA_real_,
      p_repeatability_pct = NA_real_,
      s_repeatability_pct = NA_real_,
      settle_assessment = "insufficient usable data"
    ))
  }

  cycle_means <- df_use %>%
    group_by(direction, u, cycle_id) %>%
    summarise(
      p_mean = mean(p, na.rm = TRUE),
      s_mean = mean(s, na.rm = TRUE),
      .groups = "drop"
    )

  repeatability_by_level <- cycle_means %>%
    group_by(direction, u) %>%
    summarise(
      n_cycles = n_distinct(cycle_id),
      p_sd_bar = if (n_cycles >= 3) sd(p_mean, na.rm = TRUE) else NA_real_,
      s_sd_mm = if (n_cycles >= 3) sd(s_mean, na.rm = TRUE) else NA_real_,
      p_mad_bar = if (n_cycles >= 3) mad(p_mean, constant = 1, na.rm = TRUE) else NA_real_,
      s_mad_mm = if (n_cycles >= 3) mad(s_mean, constant = 1, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  repeatability_use <- repeatability_by_level %>%
    filter(n_cycles >= 3)

  if (nrow(repeatability_use) < 6) {
    return(tibble(
      n_usable = nrow(df_use),
      window_low_v = lower_bound,
      window_high_v = upper_bound,
      p_span_global = NA_real_,
      s_span_global = NA_real_,
      p_sd_bar_p95 = NA_real_,
      s_sd_mm_p95 = NA_real_,
      p_mad_bar_median = NA_real_,
      s_mad_mm_median = NA_real_,
      p_repeatability_pct = NA_real_,
      s_repeatability_pct = NA_real_,
      settle_assessment = "insufficient usable data"
    ))
  }

  p_span_global <- diff(range(df_use$p, na.rm = TRUE))
  s_span_global <- diff(range(df_use$s, na.rm = TRUE))

  p_sd_bar_p95 <- as.numeric(quantile(repeatability_use$p_sd_bar, probs = 0.95, na.rm = TRUE, names = FALSE))
  s_sd_mm_p95 <- as.numeric(quantile(repeatability_use$s_sd_mm, probs = 0.95, na.rm = TRUE, names = FALSE))
  p_mad_bar_median <- median(repeatability_use$p_mad_bar, na.rm = TRUE)
  s_mad_mm_median <- median(repeatability_use$s_mad_mm, na.rm = TRUE)

  p_repeatability_pct <- if (is.finite(p_span_global) && p_span_global > 0) {
    100 * p_sd_bar_p95 / p_span_global
  } else {
    NA_real_
  }

  s_repeatability_pct <- if (is.finite(s_span_global) && s_span_global > 0) {
    100 * s_sd_mm_p95 / s_span_global
  } else {
    NA_real_
  }

  assessment <- case_when(
    is.na(p_repeatability_pct) || is.na(s_repeatability_pct) ~ "insufficient usable data",
    p_repeatability_pct <= 1 && s_repeatability_pct <= 3 ~ "settling likely sufficient",
    p_repeatability_pct <= 2 && s_repeatability_pct <= 5 ~ "settling possibly sufficient",
    TRUE ~ "settling may be too short or run repeatability is poor"
  )

  tibble(
    n_usable = nrow(df_use),
    window_low_v = lower_bound,
    window_high_v = upper_bound,
    p_span_global = p_span_global,
    s_span_global = s_span_global,
    p_sd_bar_p95 = p_sd_bar_p95,
    s_sd_mm_p95 = s_sd_mm_p95,
    p_mad_bar_median = p_mad_bar_median,
    s_mad_mm_median = s_mad_mm_median,
    p_repeatability_pct = p_repeatability_pct,
    s_repeatability_pct = s_repeatability_pct,
    settle_assessment = assessment
  )
}

# Summarize each full experiment run at a high level.
create_run_summary <- function(raw_data) {
  raw_data %>%
    group_by(run_id, file, step_v, settle_ms) %>%
    summarise(
      n_rows = n(),
      cycles = n_distinct(cycle_id),
      directions = n_distinct(direction),
      u_min = min(u),
      u_max = max(u),
      p_min = min(p),
      p_max = max(p),
      s_min = min(s),
      s_max = max(s),
      p_mean = mean(p),
      s_mean = mean(s),
      .groups = "drop"
    )
}

# Summarize each detected cycle and direction segment.
create_cycle_summary <- function(raw_data) {
  raw_data %>%
    group_by(run_id, cycle_id, direction) %>%
    summarise(
      n_rows = n(),
      u_min = min(u),
      u_max = max(u),
      p_max = max(p),
      s_max = max(s),
      .groups = "drop"
    )
}

# Intersect valid up/down cutoffs to get one usable window per run.
create_run_overlap_summary <- function(cutoff_summary) {
  cutoff_summary %>%
    group_by(run_id, file, step_v, settle_ms) %>%
    summarise(
      successful_directions = sum(fit_valid, na.rm = TRUE),
      stable_directions = sum(bootstrap_stable %in% TRUE, na.rm = TRUE),
      overlap_low_v = if (all(c("up", "down") %in% as.character(direction[fit_valid]))) {
        max(cutoff_low_v[fit_valid], na.rm = TRUE)
      } else {
        NA_real_
      },
      overlap_high_v = if (all(c("up", "down") %in% as.character(direction[fit_valid]))) {
        min(cutoff_high_v[fit_valid], na.rm = TRUE)
      } else {
        NA_real_
      },
      overlap_exists = !is.na(overlap_low_v) &&
        !is.na(overlap_high_v) &&
        overlap_low_v < overlap_high_v,
      overlap_width_v = if (overlap_exists) overlap_high_v - overlap_low_v else NA_real_,
      any_direction_unstable = any(fit_valid & bootstrap_enabled %in% TRUE & !bootstrap_stable),
      any_stroke_endpoint_flag = any(any_stroke_endpoint_flag %in% TRUE),
      .groups = "drop"
    )
}

# Intersect valid run windows to get one conservative global control window.
create_recommended_control_window <- function(run_overlap_summary) {
  valid_run_windows <- run_overlap_summary %>%
    filter(overlap_exists)

  if (nrow(valid_run_windows) == 0) {
    stop("No valid run-level overlap available to define a global recommended control window.")
  }

  global_window_low <- max(valid_run_windows$overlap_low_v, na.rm = TRUE)
  global_window_high <- min(valid_run_windows$overlap_high_v, na.rm = TRUE)

  if (
    !is.finite(global_window_low) ||
      !is.finite(global_window_high) ||
      global_window_low >= global_window_high
  ) {
    stop("Run-level overlaps do not share a common global control window.")
  }

  bind_rows(
    run_overlap_summary %>%
      transmute(
        level = "run",
        run_id = as.character(run_id),
        file,
        step_v,
        settle_ms,
        window_low_v = overlap_low_v,
        window_high_v = overlap_high_v,
        window_width_v = overlap_width_v,
        window_exists = overlap_exists,
        runs_used = 1L,
        any_direction_unstable,
        any_stroke_endpoint_flag
      ),
    tibble(
      level = "global",
      run_id = "global_conservative",
      file = NA_character_,
      step_v = NA_real_,
      settle_ms = NA_integer_,
      window_low_v = global_window_low,
      window_high_v = global_window_high,
      window_width_v = global_window_high - global_window_low,
      window_exists = TRUE,
      runs_used = nrow(valid_run_windows),
      any_direction_unstable = any(valid_run_windows$any_direction_unstable),
      any_stroke_endpoint_flag = any(valid_run_windows$any_stroke_endpoint_flag)
    )
  )
}

# Re-run the settling diagnostic inside the final global control window.
create_settling_summary <- function(raw_data, recommended_control_window) {
  global_window <- recommended_control_window %>%
    filter(level == "global") %>%
    slice(1)

  if (nrow(global_window) == 0) {
    return(
      raw_data %>%
        distinct(run_id, file, step_v, settle_ms) %>%
        mutate(
          n_usable = NA_integer_,
          window_low_v = NA_real_,
          window_high_v = NA_real_,
          p_span_global = NA_real_,
          s_span_global = NA_real_,
          p_sd_bar_p95 = NA_real_,
          s_sd_mm_p95 = NA_real_,
          p_mad_bar_median = NA_real_,
          s_mad_mm_median = NA_real_,
          p_repeatability_pct = NA_real_,
          s_repeatability_pct = NA_real_,
          settle_assessment = "insufficient usable data"
        )
    )
  }

  raw_data %>%
    group_split(run_id, .keep = TRUE) %>%
    map_dfr(function(df) {
      bind_cols(
        df %>% distinct(run_id, file, step_v, settle_ms),
        assess_settling(df, global_window)
      )
    })
}

# Compare valid cutoff estimates across runs and directions.
create_cutoff_validation <- function(cutoff_summary) {
  bind_rows(
    cutoff_summary %>%
      filter(fit_valid) %>%
      mutate(scope = "all"),
    cutoff_summary %>%
      filter(fit_valid) %>%
      mutate(scope = as.character(direction))
  ) %>%
    group_by(scope) %>%
    summarise(
      n_fits = n(),
      unstable_fits = sum(!bootstrap_stable, na.rm = TRUE),
      cutoff_low_min = min(cutoff_low_v, na.rm = TRUE),
      cutoff_low_max = max(cutoff_low_v, na.rm = TRUE),
      cutoff_high_min = min(cutoff_high_v, na.rm = TRUE),
      cutoff_high_max = max(cutoff_high_v, na.rm = TRUE),
      low_cutoff_range_v = cutoff_low_max - cutoff_low_min,
      high_cutoff_range_v = cutoff_high_max - cutoff_high_min,
      low_ci_overlap = if (all(is.finite(cutoff_low_ci_low)) && all(is.finite(cutoff_low_ci_high))) {
        max(cutoff_low_ci_low) <= min(cutoff_low_ci_high)
      } else {
        NA
      },
      high_ci_overlap = if (all(is.finite(cutoff_high_ci_low)) && all(is.finite(cutoff_high_ci_high))) {
        max(cutoff_high_ci_low) <= min(cutoff_high_ci_high)
      } else {
        NA
      },
      mean_rmse_p_bar = mean(rmse_p_bar, na.rm = TRUE),
      .groups = "drop"
    )
}

# Collect failed or warning-producing fits for quick inspection.
create_displacement_validation_flags <- function(cutoff_summary) {
  cutoff_summary %>%
    filter(
      any_stroke_endpoint_flag %in% TRUE |
        !fit_valid |
        (bootstrap_enabled %in% TRUE & !bootstrap_stable)
    ) %>%
    dplyr::select(
      run_id,
      direction,
      success,
      fit_valid,
      bootstrap_stable,
      cutoff_low_v,
      cutoff_high_v,
      s_cutoff_low_pct_span,
      s_cutoff_high_pct_span,
      any_stroke_endpoint_flag,
      fit_message
    )
}
