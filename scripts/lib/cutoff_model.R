# Average settled pressure and stroke at each voltage for validation summaries.
aggregate_direction_curve <- function(df) {
  df %>%
    group_by(u) %>%
    summarise(
      p_mean = mean(p, na.rm = TRUE),
      s_mean = mean(s, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(u)
}

# Choose reasonable initial breakpoints for the segmented fit.
initialize_breakpoints <- function(u_values, endpoint_margin, min_width) {
  psi_init <- as.numeric(stats::quantile(u_values, probs = c(0.25, 0.75), names = FALSE))
  psi_init[1] <- max(psi_init[1], min(u_values) + endpoint_margin)
  psi_init[2] <- min(psi_init[2], max(u_values) - endpoint_margin)

  if (!is.finite(psi_init[1]) || !is.finite(psi_init[2]) || psi_init[1] >= psi_init[2]) {
    idx <- round(c(length(u_values) / 3, 2 * length(u_values) / 3))
    idx <- pmin(pmax(idx, 2L), length(u_values) - 1L)
    psi_init <- u_values[idx]
  }

  if ((psi_init[2] - psi_init[1]) < min_width) {
    center <- mean(psi_init)
    psi_init <- c(center - min_width / 2, center + min_width / 2)
    psi_init[1] <- max(psi_init[1], min(u_values) + endpoint_margin)
    psi_init[2] <- min(psi_init[2], max(u_values) - endpoint_margin)
  }

  sort(psi_init)
}

# Fit a segmented piecewise pressure model and derive the usable cutoffs.
fit_piecewise_cutoffs <- function(
  df,
  grid_n = 400,
  min_segment_points = 3L,
  min_middle_slope_fraction = 0.10,
  min_width_fraction = 0.10,
  min_width_abs = 0.50,
  endpoint_margin_fraction = 0.05
) {
  empty_summary <- tibble(
    success = FALSE,
    fit_message = NA_character_,
    fit_method = "segmented_piecewise",
    fit_basis = "raw_pressure_samples",
    cutoff_low_v = NA_real_,
    cutoff_high_v = NA_real_,
    usable_width_v = NA_real_,
    slope_1 = NA_real_,
    slope_2 = NA_real_,
    slope_3 = NA_real_,
    rmse_p_bar = NA_real_,
    fit_valid = FALSE,
    s_cutoff_low_mm = NA_real_,
    s_cutoff_high_mm = NA_real_,
    s_min_mm = NA_real_,
    s_max_mm = NA_real_,
    s_span_mm = NA_real_,
    s_cutoff_low_pct_span = NA_real_,
    s_cutoff_high_pct_span = NA_real_,
    cutoff_low_near_stroke_endpoint = NA,
    cutoff_high_near_stroke_endpoint = NA,
    any_stroke_endpoint_flag = NA
  )

  if (nrow(df) < 6 || n_distinct(df$u) < 6) {
    return(list(
      summary = mutate(empty_summary, fit_message = "insufficient unique voltage levels"),
      curve = tibble(
        u = seq(min(df$u), max(df$u), length.out = grid_n),
        p_fit = NA_real_,
        p_fit_ci_low = NA_real_,
        p_fit_ci_high = NA_real_
      )
    ))
  }

  stroke_curve <- aggregate_direction_curve(df)
  u_values <- sort(unique(df$u))
  u_range <- range(u_values, na.rm = TRUE)
  u_span <- diff(u_range)
  min_step <- min(diff(u_values), na.rm = TRUE)
  if (!is.finite(min_step) || min_step <= 0) {
    min_step <- u_span / max(1, length(u_values) - 1)
  }
  pressure_span <- diff(range(df$p, na.rm = TRUE))
  middle_slope_min <- min_middle_slope_fraction * pressure_span /
    max(u_span, .Machine$double.eps)
  endpoint_margin <- max(endpoint_margin_fraction * u_span, 2 * min_step)
  min_width <- max(min_width_abs, min_width_fraction * u_span, 4 * min_step)

  grid <- tibble(
    u = seq(u_range[1], u_range[2], length.out = grid_n),
    p_fit = NA_real_,
    p_fit_ci_low = NA_real_,
    p_fit_ci_high = NA_real_
  )

  psi_init <- initialize_breakpoints(
    u_values = u_values,
    endpoint_margin = endpoint_margin,
    min_width = min_width
  )

  model_fit <- tryCatch(
    segmented::segmented(
      obj = lm(p ~ u, data = df),
      seg.Z = ~u,
      npsi = 2,
      psi = psi_init,
      control = segmented::seg.control(
        n.boot = 0,
        display = FALSE,
        it.max = 30,
        fix.npsi = TRUE
      )
    ),
    error = function(e) e
  )

  if (inherits(model_fit, "error")) {
    return(list(
      summary = mutate(empty_summary, fit_message = paste("segmented fit failed:", conditionMessage(model_fit))),
      curve = grid
    ))
  }

  breakpoint_est <- sort(as.numeric(model_fit$psi[, "Est."]))
  if (length(breakpoint_est) != 2 || any(!is.finite(breakpoint_est))) {
    return(list(
      summary = mutate(empty_summary, fit_message = "segmented fit did not return two finite breakpoints"),
      curve = grid
    ))
  }

  n_left <- sum(u_values <= breakpoint_est[1])
  n_mid <- sum(u_values > breakpoint_est[1] & u_values <= breakpoint_est[2])
  n_right <- sum(u_values > breakpoint_est[2])

  if (
    breakpoint_est[1] <= u_range[1] + endpoint_margin ||
      breakpoint_est[2] >= u_range[2] - endpoint_margin ||
      (breakpoint_est[2] - breakpoint_est[1]) < min_width ||
      n_left < min_segment_points ||
      n_mid < min_segment_points ||
      n_right < min_segment_points
  ) {
    return(list(
      summary = mutate(empty_summary, fit_message = "segmented breakpoints failed feasibility checks"),
      curve = grid
    ))
  }

  slope_table <- segmented::slope(model_fit)$u
  slope_est <- as.numeric(slope_table[, "Est."])

  slope_pattern_ok <- slope_est[2] > middle_slope_min &&
    slope_est[2] > slope_est[1] &&
    slope_est[2] > slope_est[3]

  if (!isTRUE(slope_pattern_ok)) {
    return(list(
      summary = mutate(empty_summary, fit_message = "segmented fit did not satisfy slope-pattern constraints"),
      curve = grid
    ))
  }

  obs_pred <- as.numeric(predict(model_fit, newdata = tibble(u = df$u)))
  rmse <- sqrt(mean((df$p - obs_pred)^2))
  grid_pred <- predict(
    model_fit,
    newdata = tibble(
      u = grid$u
    ),
    se.fit = TRUE,
    interval = "confidence",
    level = 0.95
  )
  grid$p_fit <- as.numeric(grid_pred$fit[, "fit"])
  grid$p_fit_ci_low <- as.numeric(grid_pred$fit[, "lwr"])
  grid$p_fit_ci_high <- as.numeric(grid_pred$fit[, "upr"])

  s_interp <- approx(
    x = stroke_curve$u,
    y = stroke_curve$s_mean,
    xout = breakpoint_est,
    ties = mean,
    rule = 2
  )$y
  s_min <- min(stroke_curve$s_mean, na.rm = TRUE)
  s_max <- max(stroke_curve$s_mean, na.rm = TRUE)
  s_span <- s_max - s_min
  s_pct <- if (is.finite(s_span) && s_span > 0) {
    100 * (s_interp - s_min) / s_span
  } else {
    c(NA_real_, NA_real_)
  }
  endpoint_flag <- ifelse(
    is.na(s_pct),
    NA,
    pmin(s_pct, 100 - s_pct) <= 5
  )

  list(
    summary = tibble(
      success = TRUE,
      fit_message = NA_character_,
      fit_method = "segmented_piecewise",
      fit_basis = "raw_pressure_samples",
      cutoff_low_v = breakpoint_est[1],
      cutoff_high_v = breakpoint_est[2],
      usable_width_v = breakpoint_est[2] - breakpoint_est[1],
      slope_1 = slope_est[1],
      slope_2 = slope_est[2],
      slope_3 = slope_est[3],
      rmse_p_bar = rmse,
      fit_valid = TRUE,
      s_cutoff_low_mm = s_interp[1],
      s_cutoff_high_mm = s_interp[2],
      s_min_mm = s_min,
      s_max_mm = s_max,
      s_span_mm = s_span,
      s_cutoff_low_pct_span = s_pct[1],
      s_cutoff_high_pct_span = s_pct[2],
      cutoff_low_near_stroke_endpoint = endpoint_flag[1],
      cutoff_high_near_stroke_endpoint = endpoint_flag[2],
      any_stroke_endpoint_flag = any(endpoint_flag %in% TRUE)
    ),
    curve = grid
  )
}

# Optionally bootstrap whole cycles to quantify breakpoint uncertainty.
bootstrap_piecewise_cutoffs <- function(
  df,
  B = 200,
  bootstrap_seed = 123L
) {
  cycle_ids <- sort(unique(df$cycle_id))
  if (B <= 0) {
    return(tibble(
      bootstrap_enabled = FALSE,
      bootstrap_reps = 0L,
      bootstrap_successes = NA_integer_,
      bootstrap_stable = NA,
      cutoff_low_ci_low = NA_real_,
      cutoff_low_ci_high = NA_real_,
      cutoff_high_ci_low = NA_real_,
      cutoff_high_ci_high = NA_real_
    ))
  }

  if (length(cycle_ids) == 0) {
    return(tibble(
      bootstrap_enabled = TRUE,
      bootstrap_reps = B,
      bootstrap_successes = 0L,
      bootstrap_stable = FALSE,
      cutoff_low_ci_low = NA_real_,
      cutoff_low_ci_high = NA_real_,
      cutoff_high_ci_low = NA_real_,
      cutoff_high_ci_high = NA_real_
    ))
  }

  set.seed(bootstrap_seed)
  bootstrap_results <- map_dfr(seq_len(B), function(i) {
    sampled_cycles <- sample(cycle_ids, size = length(cycle_ids), replace = TRUE)
    sampled_df <- map2_dfr(
      sampled_cycles,
      seq_along(sampled_cycles),
      function(cycle_id, draw_id) {
        df %>%
          filter(cycle_id == !!cycle_id) %>%
          mutate(cycle_id = paste0(cycle_id, "_boot_", draw_id))
      }
    )

    fit <- fit_piecewise_cutoffs(
      sampled_df,
      grid_n = 200
    )
    tibble(
      success = fit$summary$success && fit$summary$fit_valid,
      cutoff_low_v = fit$summary$cutoff_low_v,
      cutoff_high_v = fit$summary$cutoff_high_v
    )
  })

  success_n <- sum(bootstrap_results$success, na.rm = TRUE)
  valid_results <- bootstrap_results %>% filter(success)

  tibble(
    bootstrap_enabled = TRUE,
    bootstrap_reps = B,
    bootstrap_successes = success_n,
    bootstrap_stable = success_n >= 100,
    cutoff_low_ci_low = if (success_n > 0) quantile(valid_results$cutoff_low_v, 0.025) else NA_real_,
    cutoff_low_ci_high = if (success_n > 0) quantile(valid_results$cutoff_low_v, 0.975) else NA_real_,
    cutoff_high_ci_low = if (success_n > 0) quantile(valid_results$cutoff_high_v, 0.025) else NA_real_,
    cutoff_high_ci_high = if (success_n > 0) quantile(valid_results$cutoff_high_v, 0.975) else NA_real_
  )
}

# Combine the fitted cutoff summary, optional bootstrap, and fitted curve.
estimate_piecewise_direction <- function(
  df,
  bootstrap_B = 0,
  grid_n = 400,
  bootstrap_seed = 123L
) {
  fit <- fit_piecewise_cutoffs(
    df,
    grid_n = grid_n
  )
  bootstrap_summary <- if (fit$summary$success && fit$summary$fit_valid) {
    bootstrap_piecewise_cutoffs(
      df,
      B = bootstrap_B,
      bootstrap_seed = bootstrap_seed
    )
  } else {
    tibble(
      bootstrap_enabled = bootstrap_B > 0,
      bootstrap_reps = bootstrap_B,
      bootstrap_successes = if (bootstrap_B > 0) 0L else NA_integer_,
      bootstrap_stable = if (bootstrap_B > 0) FALSE else NA,
      cutoff_low_ci_low = NA_real_,
      cutoff_low_ci_high = NA_real_,
      cutoff_high_ci_low = NA_real_,
      cutoff_high_ci_high = NA_real_
    )
  }

  list(
    summary = bind_cols(fit$summary, bootstrap_summary),
    curve = fit$curve,
    mean_curve = aggregate_direction_curve(df)
  )
}

# Run the direction-level cutoff estimation for every run and branch.
estimate_direction_results <- function(
  raw_data,
  bootstrap_B,
  bootstrap_seed
) {
  direction_results <- raw_data %>%
    group_split(run_id, direction, .keep = TRUE) %>%
    map(function(df) {
      est <- estimate_piecewise_direction(
        df,
        bootstrap_B = bootstrap_B,
        bootstrap_seed = bootstrap_seed
      )
      meta <- df %>%
        distinct(run_id, file, step_v, settle_ms, direction) %>%
        slice(1)

      list(
        summary = bind_cols(meta, est$summary),
        curve = bind_cols(meta, est$curve),
        mean_curve = bind_cols(meta, est$mean_curve)
      )
    })

  list(
    cutoff_summary = bind_rows(map(direction_results, "summary")),
    cutoff_curves = bind_rows(map(direction_results, "curve")),
    cutoff_mean_curves = bind_rows(map(direction_results, "mean_curve"))
  )
}
