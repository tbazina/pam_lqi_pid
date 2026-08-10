#!/usr/bin/env Rscript

# Identify a branchwise transient pressure reference for LQI feedback.  This
# captures repeatable valve-plus-pressure dynamics after a command step; it is
# not an isolated pneumatic or material parameter model.

find_entry_script_dir <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) return(dirname(normalizePath(sub("^--file=", "", script_match[[1]]))))
  candidate <- file.path(getwd(), "scripts")
  if (dir.exists(candidate)) normalizePath(candidate) else normalizePath(getwd())
}

parse_settings <- function(defaults) {
  settings <- defaults
  for (arg in commandArgs(trailingOnly = TRUE)) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) < 2) stop(sprintf("Arguments must use --name=value: %s", arg))
    key <- pieces[[1]]
    value <- paste(pieces[-1], collapse = "=")
    if (!key %in% names(settings)) stop(sprintf("Unknown setting override: %s", key))
    if (is.integer(settings[[key]])) settings[[key]] <- as.integer(value) else if (is.numeric(settings[[key]])) settings[[key]] <- as.numeric(value) else settings[[key]] <- value
  }
  settings
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else stats::sd(x)
}

safe_slope <- function(y, time_s) {
  keep <- is.finite(y) & is.finite(time_s)
  if (sum(keep) < 2 || length(unique(y[keep])) < 2) return(0)
  as.numeric(stats::coef(stats::lm(y[keep] ~ time_s[keep]))[[2]])
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else sqrt(mean(x^2))
}

load_excitation <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing excitation file: %s", path))
  df <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(df) <- trimws(names(df))
  required <- c("u[V]", "s[mm]", "p[bar]", "t[ms]")
  if (!identical(names(df), required)) stop(sprintf("Expected columns %s", paste(required, collapse = ", ")))
  df <- df %>% dplyr::transmute(
    source_sample_id = dplyr::row_number(),
    u_cmd_v = as.numeric(`u[V]`),
    s_meas_mm = as.numeric(`s[mm]`),
    p_meas_bar = as.numeric(`p[bar]`),
    loop_dt_ms = as.numeric(`t[ms]`)
  )
  if (any(!is.finite(unlist(df)))) stop("Excitation input contains non-finite values.")
  df
}

infer_sampling_time <- function(df) {
  valid <- df$loop_dt_ms[df$source_sample_id > 2L & df$loop_dt_ms > 0 & df$loop_dt_ms < 1000]
  if (length(valid) == 0) stop("Could not infer sample time from valid loop periods.")
  stats::median(valid) / 1000
}

segment_holds <- function(df, ts_s) {
  df <- df %>% dplyr::mutate(
    hold_id = cumsum(dplyr::row_number() == 1L | abs(u_cmd_v - dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v))) > 1e-12)
  )
  hold_summary <- df %>%
    dplyr::group_by(hold_id) %>%
    dplyr::summarise(n_samples = dplyr::n(), u_cmd_v = dplyr::first(u_cmd_v), .groups = "drop") %>%
    dplyr::mutate(
      u_prev_v = dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v)),
      delta_u_v = u_cmd_v - u_prev_v,
      branch = dplyr::case_when(delta_u_v > 1e-12 ~ "up", delta_u_v < -1e-12 ~ "down", TRUE ~ "start")
    )
  full_modal_n <- hold_summary %>%
    dplyr::filter(branch %in% c("up", "down"), n_samples > 1L) %>%
    dplyr::count(n_samples, sort = TRUE) %>% dplyr::slice(1) %>% dplyr::pull(n_samples)
  if (length(full_modal_n) != 1) stop("Could not determine modal full-hold length.")
  df %>%
    dplyr::left_join(hold_summary %>% dplyr::select(-u_cmd_v), by = "hold_id") %>%
    dplyr::group_by(hold_id) %>%
    dplyr::mutate(sample_in_hold = dplyr::row_number(), time_since_step_s = (sample_in_hold - 1L) * ts_s) %>%
    dplyr::ungroup() %>% dplyr::mutate(full_modal_n = as.integer(full_modal_n))
}

summarise_holds <- function(segmented_df, settings) {
  split_holds <- split(segmented_df, segmented_df$hold_id)
  dplyr::bind_rows(lapply(split_holds, function(hold_df) {
    n_samples <- nrow(hold_df)
    tail_n <- max(5L, ceiling(settings$terminal_tail_fraction * n_samples))
    tail_df <- utils::tail(hold_df, tail_n)
    s_start <- hold_df$s_meas_mm[[1]]
    p_start <- hold_df$p_meas_bar[[1]]
    s_terminal <- mean(tail_df$s_meas_mm)
    p_terminal <- mean(tail_df$p_meas_bar)
    is_full <- n_samples >= settings$min_full_hold_samples && abs(n_samples - hold_df$full_modal_n[[1]]) <= settings$full_hold_count_tolerance
    terminal_stable <-
      safe_sd(tail_df$s_meas_mm) <= settings$max_terminal_s_sd_mm &&
      abs(safe_slope(tail_df$s_meas_mm, tail_df$time_since_step_s)) <= settings$max_terminal_s_slope_mm_s &&
      safe_sd(tail_df$p_meas_bar) <= settings$max_terminal_p_sd_bar &&
      abs(safe_slope(tail_df$p_meas_bar, tail_df$time_since_step_s)) <= settings$max_terminal_p_slope_bar_s
    displacement_step <- s_terminal - s_start
    pressure_step <- p_terminal - p_start
    reason <- dplyr::case_when(
      hold_df$branch[[1]] == "start" ~ "startup_or_zero_command_hold",
      !is_full ~ "incomplete_hold",
      !terminal_stable ~ "unstable_terminal_tail",
      abs(hold_df$delta_u_v[[1]]) < settings$min_abs_delta_u_v ~ "small_command_step",
      abs(displacement_step) < settings$min_abs_displacement_step_mm ~ "small_displacement_step",
      abs(pressure_step) < settings$min_abs_pressure_step_bar ~ "small_pressure_step",
      TRUE ~ ""
    )
    tibble::tibble(
      hold_id = hold_df$hold_id[[1]], branch = hold_df$branch[[1]], n_samples = n_samples,
      full_modal_n = hold_df$full_modal_n[[1]], is_full_hold = is_full,
      u_prev_v = hold_df$u_prev_v[[1]], u_cmd_v = hold_df$u_cmd_v[[1]], delta_u_v = hold_df$delta_u_v[[1]],
      duration_s = n_samples * settings$ts_s, tail_n_samples = tail_n,
      s_start_mm = s_start, s_terminal_mm = s_terminal, p_start_bar = p_start, p_terminal_bar = p_terminal,
      displacement_step_mm = displacement_step, pressure_step_bar = pressure_step,
      s_terminal_sd_mm = safe_sd(tail_df$s_meas_mm), p_terminal_sd_bar = safe_sd(tail_df$p_meas_bar),
      s_terminal_slope_mm_s = safe_slope(tail_df$s_meas_mm, tail_df$time_since_step_s),
      p_terminal_slope_bar_s = safe_slope(tail_df$p_meas_bar, tail_df$time_since_step_s),
      accepted_for_pressure_reference = identical(reason, ""), exclusion_reason = dplyr::na_if(reason, "")
    )
  }))
}

attach_hold_metrics <- function(segmented_df, hold_metrics) {
  segmented_df %>%
    dplyr::left_join(
      hold_metrics %>% dplyr::select(hold_id, branch, accepted_for_pressure_reference, exclusion_reason, p_start_bar, p_terminal_bar, pressure_step_bar),
      by = c("hold_id", "branch")
    ) %>%
    dplyr::mutate(
      terminal_relative_p_bar = p_meas_bar - p_terminal_bar,
      normalized_pressure_progress = dplyr::if_else(
        accepted_for_pressure_reference,
        (p_meas_bar - p_start_bar) / pressure_step_bar,
        NA_real_
      )
    )
}

profile_by_time <- function(cleaned_samples) {
  cleaned_samples %>%
    dplyr::filter(accepted_for_pressure_reference, branch %in% c("up", "down"), is.finite(normalized_pressure_progress)) %>%
    dplyr::group_by(branch, time_since_step_s) %>%
    dplyr::summarise(
      n_samples = dplyr::n(), n_holds = dplyr::n_distinct(hold_id),
      pressure_progress_median = stats::median(normalized_pressure_progress),
      pressure_progress_q10 = as.numeric(stats::quantile(normalized_pressure_progress, 0.10, names = FALSE)),
      pressure_progress_q90 = as.numeric(stats::quantile(normalized_pressure_progress, 0.90, names = FALSE)),
      terminal_relative_p_median_bar = stats::median(terminal_relative_p_bar),
      terminal_relative_p_q10_bar = as.numeric(stats::quantile(terminal_relative_p_bar, 0.10, names = FALSE)),
      terminal_relative_p_q90_bar = as.numeric(stats::quantile(terminal_relative_p_bar, 0.90, names = FALSE)),
      .groups = "drop"
    )
}

predict_pressure_model <- function(model, parameters, time_s) {
  t <- pmax(as.numeric(time_s), 0)
  if (model == "static_immediate") return(rep(1, length(t)))
  if (model == "first_order") return(1 - exp(-t / parameters[["tau_s"]]))
  if (model == "fopdt") {
    elapsed <- pmax(t - parameters[["theta_s"]], 0)
    return(ifelse(t < parameters[["theta_s"]], 0, 1 - exp(-elapsed / parameters[["tau_s"]])))
  }
  if (model == "fopdt_double_exponential") {
    elapsed <- pmax(t - parameters[["theta_s"]], 0)
    response <- 1 - parameters[["weight_fast"]] * exp(-elapsed / parameters[["tau_fast_s"]]) -
      (1 - parameters[["weight_fast"]]) * exp(-elapsed / parameters[["tau_slow_s"]])
    return(ifelse(t < parameters[["theta_s"]], 0, response))
  }
  stop(sprintf("Unknown pressure model: %s", model))
}

fit_by_optim <- function(model, profile_df, settings) {
  x <- profile_df$time_since_step_s
  y <- profile_df$pressure_progress_median
  objective <- function(par) sum((y - predict_pressure_model(model, par, x))^2)
  if (model == "static_immediate") return(list(model = model, parameters = c(theta_s = NA_real_, tau_s = NA_real_, tau_fast_s = NA_real_, tau_slow_s = NA_real_, weight_fast = NA_real_)))
  if (model == "first_order") {
    fit <- stats::optim(c(tau_s = 0.15), objective, method = "L-BFGS-B", lower = c(tau_s = 0.01), upper = c(tau_s = settings$max_tau_s))
    return(list(model = model, parameters = c(theta_s = 0, tau_s = fit$par[["tau_s"]], tau_fast_s = NA_real_, tau_slow_s = NA_real_, weight_fast = NA_real_)))
  }
  if (model == "fopdt") {
    fit <- stats::optim(c(theta_s = 0.08, tau_s = 0.12), objective, method = "L-BFGS-B",
      lower = c(theta_s = 0, tau_s = 0.01), upper = c(theta_s = settings$max_dead_time_s, tau_s = settings$max_tau_s))
    return(list(model = model, parameters = c(theta_s = fit$par[["theta_s"]], tau_s = fit$par[["tau_s"]], tau_fast_s = NA_real_, tau_slow_s = NA_real_, weight_fast = NA_real_)))
  }
  if (model == "fopdt_double_exponential") {
    fit <- stats::optim(c(theta_s = 0.08, tau_fast_s = 0.06, tau_slow_s = 0.24, weight_fast = 0.65), objective, method = "L-BFGS-B",
      lower = c(theta_s = 0, tau_fast_s = 0.01, tau_slow_s = 0.01, weight_fast = 0),
      upper = c(theta_s = settings$max_dead_time_s, tau_fast_s = settings$max_tau_s, tau_slow_s = settings$max_tau_s, weight_fast = 1))
    return(list(model = model, parameters = c(theta_s = fit$par[["theta_s"]], tau_s = NA_real_, tau_fast_s = fit$par[["tau_fast_s"]], tau_slow_s = fit$par[["tau_slow_s"]], weight_fast = fit$par[["weight_fast"]])))
  }
  NULL
}

model_validity <- function(fit, settings) {
  time_grid <- seq(0, settings$runtime_horizon_s, by = 0.001)
  predicted <- predict_pressure_model(fit$model, fit$parameters, time_grid)
  finite <- all(is.finite(predicted))
  bounded <- finite && all(predicted >= -1e-9 & predicted <= 1 + 1e-9)
  nondecreasing <- finite && all(diff(predicted) >= -1e-7)
  list(valid = finite && bounded && nondecreasing, finite = finite, bounded = bounded, nondecreasing = nondecreasing)
}

make_folds <- function(hold_metrics, k_folds) {
  hold_metrics %>%
    dplyr::filter(accepted_for_pressure_reference) %>%
    dplyr::group_by(branch) %>%
    dplyr::arrange(abs(delta_u_v), .by_group = TRUE) %>%
    dplyr::mutate(cv_fold = rep(seq_len(k_folds), length.out = dplyr::n())) %>%
    dplyr::ungroup() %>% dplyr::select(hold_id, cv_fold)
}

cross_validate_models <- function(cleaned_samples, folds, settings) {
  model_names <- c("static_immediate", "first_order", "fopdt", "fopdt_double_exponential")
  complexity <- c(static_immediate = 0L, first_order = 1L, fopdt = 2L, fopdt_double_exponential = 4L)
  dplyr::bind_rows(lapply(c("up", "down"), function(branch) {
    samples <- cleaned_samples %>% dplyr::filter(accepted_for_pressure_reference, branch == !!branch) %>% dplyr::inner_join(folds, by = "hold_id")
    dplyr::bind_rows(lapply(model_names, function(model) {
      fold_errors <- vapply(seq_len(settings$cv_folds), function(fold) {
        train <- samples %>% dplyr::filter(cv_fold != fold)
        test <- samples %>% dplyr::filter(cv_fold == fold)
        train_profile <- train %>% dplyr::group_by(time_since_step_s) %>% dplyr::summarise(pressure_progress_median = stats::median(normalized_pressure_progress), .groups = "drop")
        fit <- fit_by_optim(model, train_profile, settings)
        if (is.null(fit) || !model_validity(fit, settings)$valid) return(NA_real_)
        safe_rmse(test$normalized_pressure_progress - predict_pressure_model(model, fit$parameters, test$time_since_step_s))
      }, numeric(1))
      full_profile <- samples %>% dplyr::group_by(time_since_step_s) %>% dplyr::summarise(pressure_progress_median = stats::median(normalized_pressure_progress), .groups = "drop")
      full_fit <- fit_by_optim(model, full_profile, settings)
      validity <- model_validity(full_fit, settings)
      tibble::tibble(
        branch = branch, model = model, model_complexity = unname(complexity[[model]]),
        cv_rmse = if (any(is.finite(fold_errors))) mean(fold_errors, na.rm = TRUE) else NA_real_,
        cv_rmse_sd = if (sum(is.finite(fold_errors)) > 1) stats::sd(fold_errors, na.rm = TRUE) else NA_real_,
        cv_folds_completed = sum(is.finite(fold_errors)),
        theta_s = full_fit$parameters[["theta_s"]], tau_s = full_fit$parameters[["tau_s"]],
        tau_fast_s = full_fit$parameters[["tau_fast_s"]], tau_slow_s = full_fit$parameters[["tau_slow_s"]], weight_fast = full_fit$parameters[["weight_fast"]],
        valid_runtime_profile = validity$valid, runtime_profile_finite = validity$finite,
        runtime_profile_bounded = validity$bounded, runtime_profile_nondecreasing = validity$nondecreasing
      )
    }))
  }))
}

select_models <- function(model_comparison, simplicity_margin) {
  model_comparison %>%
    dplyr::mutate(eligible_for_runtime = valid_runtime_profile & is.finite(cv_rmse) & cv_folds_completed >= 3L) %>%
    dplyr::group_by(branch) %>%
    dplyr::mutate(
      best_eligible_cv_rmse = min(cv_rmse[eligible_for_runtime], na.rm = TRUE),
      within_simplicity_margin = eligible_for_runtime & cv_rmse <= (1 + simplicity_margin) * best_eligible_cv_rmse,
      selection_order = dplyr::if_else(within_simplicity_margin, model_complexity * 1e6 + cv_rmse, Inf),
      selected_model = dplyr::row_number() == which.min(selection_order)
    ) %>% dplyr::ungroup() %>% dplyr::select(-selection_order)
}

build_runtime_profiles <- function(selected_models, settings, calibrated_horizon_s) {
  time_grid <- seq(0, settings$runtime_horizon_s, by = settings$profile_dt_s)
  dplyr::bind_rows(lapply(c("up", "down"), function(branch) {
    row <- selected_models %>% dplyr::filter(branch == !!branch, selected_model)
    if (nrow(row) != 1) stop(sprintf("Missing selected pressure model for %s.", branch))
    parameters <- c(theta_s = row$theta_s[[1]], tau_s = row$tau_s[[1]], tau_fast_s = row$tau_fast_s[[1]], tau_slow_s = row$tau_slow_s[[1]], weight_fast = row$weight_fast[[1]])
    tibble::tibble(
      branch = branch, time_since_step_s = time_grid,
      pressure_progress_fraction = predict_pressure_model(row$model[[1]], parameters, time_grid),
      selected_model = row$model[[1]], theta_s = row$theta_s[[1]], tau_s = row$tau_s[[1]],
      tau_fast_s = row$tau_fast_s[[1]], tau_slow_s = row$tau_slow_s[[1]], weight_fast = row$weight_fast[[1]],
      calibrated_horizon_s = calibrated_horizon_s, runtime_horizon_s = settings$runtime_horizon_s,
      interval = dplyr::if_else(time_grid <= calibrated_horizon_s + 1e-12, "calibrated_open_loop_hold", "model_extension")
    )
  }))
}

runtime_formula_for_model <- function(model) {
  switch(model,
    static_immediate = "g=1;_p_ff_dynamic=p_ff_static;_p_tilde=p_meas-p_ff_dynamic",
    first_order = "g=1-exp(-t/tau);_p_ff_dynamic=p_step+(p_ff_static-p_step)*g;_p_tilde=p_meas-p_ff_dynamic",
    fopdt = "g=0_for_t<theta;_else_1-exp(-(t-theta)/tau);_p_ff_dynamic=p_step+(p_ff_static-p_step)*g;_p_tilde=p_meas-p_ff_dynamic",
    fopdt_double_exponential = "g=0_for_t<theta;_else_1-w*exp(-(t-theta)/tau_fast)-(1-w)*exp(-(t-theta)/tau_slow);_p_ff_dynamic=p_step+(p_ff_static-p_step)*g;_p_tilde=p_meas-p_ff_dynamic",
    stop(sprintf("No runtime formula defined for model: %s", model))
  )
}

create_plots <- function(profile_df, hold_metrics, model_comparison, runtime_profiles, settings, plot_theme) {
  profile_fit <- ggplot2::ggplot(profile_df, ggplot2::aes(time_since_step_s, pressure_progress_median, colour = branch, fill = branch)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pressure_progress_q10, ymax = pressure_progress_q90), alpha = 0.16, colour = NA) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_line(data = runtime_profiles, ggplot2::aes(time_since_step_s, pressure_progress_fraction, colour = branch, group = branch), linetype = "dashed", linewidth = 0.50, inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = max(profile_df$time_since_step_s), linetype = "dotdash") +
    ggsci::scale_color_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::labs(x = "Time since voltage step [s]", y = "Normalized pressure progress [-]", colour = "Branch", fill = "Branch") + plot_theme
  comparison <- model_comparison %>% dplyr::filter(is.finite(cv_rmse)) %>%
    ggplot2::ggplot(ggplot2::aes(model, cv_rmse, fill = branch)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_point(data = model_comparison %>% dplyr::filter(selected_model), ggplot2::aes(model, cv_rmse), inherit.aes = FALSE, shape = 21, fill = "white", size = 1.8, stroke = 0.25) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::scale_x_discrete(labels = c(static_immediate = "Static target", first_order = "First-order", fopdt = "FOPDT", fopdt_double_exponential = "Double exponential FOPDT")) +
    ggplot2::labs(x = "Candidate pressure-reference model", y = "Normalized-pressure CV RMSE", fill = "Branch") + plot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 9))
  quality <- hold_metrics %>%
    dplyr::mutate(status = dplyr::if_else(accepted_for_pressure_reference, "accepted", dplyr::coalesce(exclusion_reason, "excluded"))) %>%
    dplyr::count(branch, status) %>% dplyr::filter(branch %in% c("up", "down")) %>%
    ggplot2::ggplot(ggplot2::aes(status, n, fill = branch)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::scale_x_discrete(labels = c(accepted = "Accepted", incomplete_hold = "Incomplete hold", small_command_step = "Small command step", small_displacement_step = "Small displacement step", startup_or_zero_command_hold = "Startup or zero-command hold", insufficient_pressure_change = "Insufficient pressure change", excluded = "Excluded")) +
    ggplot2::labs(x = "Hold outcome", y = "Number of holds", fill = "Branch") + plot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  selected_parameters <- model_comparison %>%
    dplyr::filter(selected_model) %>%
    dplyr::select(branch, model, theta_s, tau_s, tau_fast_s, tau_slow_s, weight_fast)
  residuals <- profile_df %>%
    dplyr::left_join(selected_parameters, by = "branch") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pressure_progress_model = predict_pressure_model(
        model,
        c(theta_s = theta_s, tau_s = tau_s, tau_fast_s = tau_fast_s, tau_slow_s = tau_slow_s, weight_fast = weight_fast),
        time_since_step_s
      )[[1]],
      residual = pressure_progress_median - pressure_progress_model
    ) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(time_since_step_s, residual, colour = branch)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.2) + ggplot2::geom_line(linewidth = 0.45) +
    ggsci::scale_color_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::labs(x = "Time since voltage step [s]", y = "Observed minus model progress [-]", colour = "Branch") + plot_theme
  list(profile_fit = profile_fit, comparison = comparison, quality = quality, residuals = residuals)
}

main <- function() {
  script_dir <- find_entry_script_dir()
  source(file.path(script_dir, "lib", "setup.R"), local = FALSE)
  load_analysis_packages()
  context <- initialize_analysis_context()
  settings <- parse_settings(list(
    input_csv = "experiment/excitation_experiment_100hz_15min_4s_hold.csv",
    terminal_tail_fraction = 0.05, min_full_hold_samples = 380L, full_hold_count_tolerance = 2L,
    min_abs_delta_u_v = 0.35, min_abs_displacement_step_mm = 1.0, min_abs_pressure_step_bar = 0.10,
    max_terminal_s_sd_mm = 0.60, max_terminal_s_slope_mm_s = 1.20,
    max_terminal_p_sd_bar = 0.12, max_terminal_p_slope_bar_s = 0.25,
    cv_folds = 5L, simplicity_margin = 0.05, max_dead_time_s = 0.50, max_tau_s = 3.99,
    # A 4 s hold sampled at 100 Hz has measured samples from 0 to 3.99 s.
    runtime_horizon_s = 3.99, profile_dt_s = 0.01
  ))
  input_path <- file.path(context$repo_dir, settings$input_csv)
  raw_df <- load_excitation(input_path)
  settings$ts_s <- infer_sampling_time(raw_df)
  segmented_df <- segment_holds(raw_df, settings$ts_s)
  hold_metrics <- summarise_holds(segmented_df, settings)
  cleaned_samples <- attach_hold_metrics(segmented_df, hold_metrics)
  profile_df <- profile_by_time(cleaned_samples)
  if (nrow(profile_df) == 0) stop("No accepted pressure-transition holds remain after filtering.")
  folds <- make_folds(hold_metrics, settings$cv_folds)
  model_comparison <- cross_validate_models(cleaned_samples, folds, settings) %>% select_models(settings$simplicity_margin)
  if (sum(model_comparison$selected_model) != 2L) stop("Exactly one pressure model must be selected per branch.")
  calibrated_horizon_s <- max(profile_df$time_since_step_s)
  runtime_profiles <- build_runtime_profiles(model_comparison, settings, calibrated_horizon_s)
  if (any(!is.finite(runtime_profiles$pressure_progress_fraction)) || any(runtime_profiles$pressure_progress_fraction < -1e-9) || any(runtime_profiles$pressure_progress_fraction > 1 + 1e-9)) stop("Invalid runtime pressure profile.")
  if (any(runtime_profiles %>% dplyr::group_by(branch) %>% dplyr::summarise(bad = any(diff(pressure_progress_fraction) < -1e-7), .groups = "drop") %>% dplyr::pull(bad))) stop("Runtime pressure profile must be nondecreasing.")

  summary_df <- model_comparison %>% dplyr::filter(selected_model) %>%
    dplyr::left_join(hold_metrics %>% dplyr::group_by(branch) %>% dplyr::summarise(
      n_holds_total = dplyr::n(), n_holds_full = sum(is_full_hold), n_holds_accepted = sum(accepted_for_pressure_reference), n_holds_excluded = sum(!accepted_for_pressure_reference), .groups = "drop"
    ), by = "branch") %>%
    dplyr::mutate(
      input_dataset = settings$input_csv, sampling_time_s = settings$ts_s, full_hold_samples = unique(segmented_df$full_modal_n)[[1]],
      calibrated_horizon_s = calibrated_horizon_s, runtime_horizon_s = settings$runtime_horizon_s,
      controller_scope = "lookup_plus_branchwise_minimal_LQI|lookup_plus_branchwise_minimal_LQI_VEL",
      implementation_policy = "direct_selected_model_equation;_profiles_are_verification_only"
    )
  static_rows <- model_comparison %>% dplyr::filter(model == "static_immediate") %>% dplyr::select(branch, static_immediate_cv_rmse = cv_rmse)
  summary_df <- summary_df %>% dplyr::left_join(static_rows, by = "branch") %>% dplyr::mutate(dynamic_improves_static = cv_rmse < static_immediate_cv_rmse)
  if (any(!summary_df$dynamic_improves_static)) stop("Selected dynamic model does not improve on the static-immediate baseline.")

  runtime_contract <- summary_df %>% dplyr::transmute(
    branch, controller_scope, selected_model = model,
    required_runtime_inputs = "conditioned_p_meas_bar_at_step,conditioned_p_meas_bar,s_ref_mm,s_ref_prev_mm,selected_branch,time_since_material_reference_step_s",
    material_step_rule = "reset_and_cache_when_abs(s_ref_mm-s_ref_prev_mm)_exceeds_existing_branch_deadband_mm",
    formula = vapply(model, runtime_formula_for_model, character(1)),
    static_lookup_policy = "interpolate_and_cache_open_loop_steady_state_p_lookup(active_branch,s_ref_mm)_at_material_step;_post_correction_deferred",
    branch_policy = "apply_only_when_material_step_direction_matches_selected_branch; otherwise use_static_p_ff",
    horizon_policy = "use_g=1_after_runtime_horizon_s; do_not_extend_model_beyond_exported_horizon",
    interaction_with_creep = "independent_of_u_ff_creep_layer; modifies_only_p_tilde",
    theta_s, tau_s, tau_fast_s, tau_slow_s, weight_fast,
    calibrated_horizon_s, runtime_horizon_s
  )

  readr::write_csv(cleaned_samples, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_cleaned_samples.csv"))
  readr::write_csv(hold_metrics, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_hold_metrics.csv"))
  readr::write_csv(profile_df, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_profile_points.csv"))
  readr::write_csv(model_comparison, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_model_comparison.csv"))
  readr::write_csv(summary_df, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_summary.csv"))
  readr::write_csv(runtime_contract, file.path(context$output_dir, "open_loop_dynamic_pressure_reference_runtime_contract.csv"))
  for (branch in c("up", "down")) readr::write_csv(runtime_profiles %>% dplyr::filter(branch == !!branch), file.path(context$output_dir, sprintf("open_loop_dynamic_pressure_reference_%s.csv", branch)))

  plot_theme <- create_analysis_plot_theme()
  plots <- create_plots(profile_df, hold_metrics, model_comparison, runtime_profiles, settings, plot_theme)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_profile.png"), plots$profile_fit, width = 15, height = 8.4375, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_model_comparison.png"), plots$comparison, width = 15, height = 8.4375, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_hold_quality.png"), plots$quality, width = 15, height = 8.4375, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_residuals.png"), plots$residuals, width = 15, height = 7.5, units = "cm", dpi = 600)

  cat("Open-loop dynamic pressure-reference analysis complete\n")
  cat(sprintf("Input: %s\n", input_path))
  cat(sprintf("Sampling time: %.6f s; modal full hold: %d samples; calibrated horizon: %.3f s\n", settings$ts_s, unique(segmented_df$full_modal_n)[[1]], calibrated_horizon_s))
  print(summary_df %>% dplyr::select(branch, n_holds_accepted, model, cv_rmse, static_immediate_cv_rmse, theta_s, tau_s, dynamic_improves_static))
  cat(sprintf("Outputs written to: %s\n", context$output_dir))
}

main()
