#!/usr/bin/env Rscript

# Estimate repeatable branchwise within-hold relaxation from the open-loop
# excitation experiment. The result is a dynamic feedforward layer, not a
# material-only viscoelastic parameter estimate: valve and pressure dynamics
# remain part of the observed response.

find_entry_script_dir <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_match[[1]]))))
  }
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

safe_slope <- function(y, time_s) {
  keep <- is.finite(y) & is.finite(time_s)
  if (sum(keep) < 2 || length(unique(y[keep])) < 2) return(0)
  as.numeric(stats::coef(stats::lm(y[keep] ~ time_s[keep]))[[2]])
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else sqrt(mean(x^2))
}

read_lookup <- function(path, expected_branch) {
  if (!file.exists(path)) stop(sprintf("Missing lookup table: %s", path))
  lookup <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  required <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")
  if (!identical(names(lookup), required)) stop(sprintf("Unexpected lookup schema: %s", path))
  if (any(lookup$branch != expected_branch) || any(!is.finite(unlist(lookup[required[-1]])))) {
    stop(sprintf("Invalid %s lookup table: %s", expected_branch, path))
  }
  lookup <- lookup %>% dplyr::arrange(s_ref_mm)
  if (anyDuplicated(lookup$s_ref_mm) || any(diff(lookup$s_ref_mm) <= 0)) {
    stop(sprintf("Lookup references must be unique and ascending: %s", path))
  }
  lookup
}

interpolate_lookup <- function(lookup, s_ref_mm, value_col = "u_ff_v") {
  stats::approx(lookup$s_ref_mm, lookup[[value_col]], xout = s_ref_mm, method = "linear", rule = 1)$y
}

load_excitation <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing excitation experiment: %s", path))
  df <- readr::read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(df) <- trimws(names(df))
  required <- c("u[V]", "s[mm]", "p[bar]", "t[ms]")
  if (!identical(names(df), required)) {
    stop(sprintf("Expected columns %s in %s", paste(required, collapse = ", "), path))
  }
  df <- df %>%
    dplyr::transmute(
      source_sample_id = dplyr::row_number(),
      u_cmd_v = as.numeric(`u[V]`),
      s_meas_mm = as.numeric(`s[mm]`),
      p_meas_bar = as.numeric(`p[bar]`),
      loop_dt_ms = as.numeric(`t[ms]`)
    )
  if (any(!is.finite(unlist(df[c("u_cmd_v", "s_meas_mm", "p_meas_bar", "loop_dt_ms")]))) ) {
    stop("The excitation file contains non-finite values.")
  }
  df
}

infer_sampling_time <- function(df) {
  valid <- df$loop_dt_ms[df$source_sample_id > 2L & df$loop_dt_ms > 0 & df$loop_dt_ms < 1000]
  if (length(valid) == 0) stop("Could not infer a sampling period from valid loop times.")
  stats::median(valid) / 1000
}

segment_holds <- function(df, ts_s) {
  df <- df %>%
    dplyr::mutate(hold_id = cumsum(dplyr::row_number() == 1L | abs(u_cmd_v - dplyr::lag(u_cmd_v, default = dplyr::first(u_cmd_v))) > 1e-12))

  hold_summary <- df %>%
    dplyr::group_by(hold_id) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      u_cmd_v = dplyr::first(u_cmd_v),
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

  full_modal_n <- hold_summary %>%
    dplyr::filter(branch %in% c("up", "down"), n_samples > 1L) %>%
    dplyr::count(n_samples, sort = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::pull(n_samples)
  if (length(full_modal_n) != 1) stop("Could not determine the modal full-hold sample count.")

  df %>%
    # The sample table already contains u_cmd_v. Add only derived hold fields
    # to avoid a suffixed duplicate command column in audit outputs.
    dplyr::left_join(hold_summary %>% dplyr::select(-u_cmd_v), by = "hold_id") %>%
    dplyr::group_by(hold_id) %>%
    dplyr::mutate(
      sample_in_hold = dplyr::row_number(),
      time_since_step_s = (sample_in_hold - 1L) * ts_s
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(full_modal_n = as.integer(full_modal_n))
}

summarise_holds <- function(segmented_df, settings) {
  split_holds <- split(segmented_df, segmented_df$hold_id)
  dplyr::bind_rows(lapply(split_holds, function(hold_df) {
    n_samples <- nrow(hold_df)
    tail_n <- max(5L, ceiling(settings$terminal_tail_fraction * n_samples))
    tail_df <- utils::tail(hold_df, tail_n)
    duration_s <- n_samples * settings$ts_s
    s_start <- hold_df$s_meas_mm[[1]]
    s_terminal <- mean(tail_df$s_meas_mm)
    p_terminal <- mean(tail_df$p_meas_bar)
    is_full <- n_samples >= settings$min_full_hold_samples && abs(n_samples - hold_df$full_modal_n[[1]]) <= settings$full_hold_count_tolerance
    terminal_stable <-
      safe_sd(tail_df$s_meas_mm) <= settings$max_terminal_s_sd_mm &&
      abs(safe_slope(tail_df$s_meas_mm, tail_df$time_since_step_s)) <= settings$max_terminal_s_slope_mm_s &&
      safe_sd(tail_df$p_meas_bar) <= settings$max_terminal_p_sd_bar &&
      abs(safe_slope(tail_df$p_meas_bar, tail_df$time_since_step_s)) <= settings$max_terminal_p_slope_bar_s
    movement_mm <- s_terminal - s_start
    reason <- dplyr::case_when(
      hold_df$branch[[1]] == "start" ~ "startup_or_zero_command_hold",
      !is_full ~ "incomplete_hold",
      !terminal_stable ~ "unstable_terminal_tail",
      abs(hold_df$delta_u_v[[1]]) < settings$min_abs_delta_u_v ~ "small_command_step",
      abs(movement_mm) < settings$min_abs_displacement_step_mm ~ "small_displacement_step",
      TRUE ~ ""
    )
    tibble::tibble(
      hold_id = hold_df$hold_id[[1]],
      branch = hold_df$branch[[1]],
      n_samples = n_samples,
      full_modal_n = hold_df$full_modal_n[[1]],
      is_full_hold = is_full,
      u_prev_v = hold_df$u_prev_v[[1]],
      u_cmd_v = hold_df$u_cmd_v[[1]],
      delta_u_v = hold_df$delta_u_v[[1]],
      duration_s = duration_s,
      tail_n_samples = tail_n,
      s_start_mm = s_start,
      s_terminal_mm = s_terminal,
      p_start_bar = hold_df$p_meas_bar[[1]],
      p_terminal_bar = p_terminal,
      displacement_step_mm = movement_mm,
      s_terminal_sd_mm = safe_sd(tail_df$s_meas_mm),
      p_terminal_sd_bar = safe_sd(tail_df$p_meas_bar),
      s_terminal_slope_mm_s = safe_slope(tail_df$s_meas_mm, tail_df$time_since_step_s),
      p_terminal_slope_bar_s = safe_slope(tail_df$p_meas_bar, tail_df$time_since_step_s),
      accepted_for_creep = identical(reason, ""),
      exclusion_reason = dplyr::na_if(reason, "")
    )
  }))
}

attach_hold_metrics <- function(segmented_df, hold_metrics) {
  segmented_df %>%
    dplyr::left_join(
      hold_metrics %>% dplyr::select(hold_id, branch, accepted_for_creep, exclusion_reason, s_start_mm, s_terminal_mm, p_terminal_bar, displacement_step_mm),
      by = c("hold_id", "branch")
    ) %>%
    dplyr::mutate(
      terminal_relative_s_mm = s_meas_mm - s_terminal_mm,
      terminal_relative_p_bar = p_meas_bar - p_terminal_bar,
      normalized_remaining_displacement = dplyr::if_else(
        accepted_for_creep,
        (s_terminal_mm - s_meas_mm) / displacement_step_mm,
        NA_real_
      )
    )
}

profile_by_time <- function(cleaned_samples) {
  cleaned_samples %>%
    dplyr::filter(accepted_for_creep, branch %in% c("up", "down")) %>%
    dplyr::group_by(branch, time_since_step_s) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_holds = dplyr::n_distinct(hold_id),
      normalized_remaining_median = stats::median(normalized_remaining_displacement),
      normalized_remaining_q10 = as.numeric(stats::quantile(normalized_remaining_displacement, 0.10, names = FALSE)),
      normalized_remaining_q90 = as.numeric(stats::quantile(normalized_remaining_displacement, 0.90, names = FALSE)),
      terminal_relative_s_median_mm = stats::median(terminal_relative_s_mm),
      terminal_relative_p_median_bar = stats::median(terminal_relative_p_bar),
      terminal_relative_p_q10_bar = as.numeric(stats::quantile(terminal_relative_p_bar, 0.10, names = FALSE)),
      terminal_relative_p_q90_bar = as.numeric(stats::quantile(terminal_relative_p_bar, 0.90, names = FALSE)),
      .groups = "drop"
    )
}

find_pressure_plateau <- function(profile_df, settings) {
  branches <- c("up", "down")
  dplyr::bind_rows(lapply(branches, function(branch) {
    branch_profile <- profile_df %>% dplyr::filter(branch == !!branch) %>% dplyr::arrange(time_since_step_s)
    candidate_idx <- which(branch_profile$time_since_step_s >= settings$minimum_creep_start_s)
    plateau_idx <- NA_integer_
    if (length(candidate_idx) > 0) {
      for (idx in candidate_idx) {
        later <- branch_profile$terminal_relative_p_median_bar[idx:nrow(branch_profile)]
        if (all(abs(later) <= settings$pressure_plateau_band_bar)) {
          plateau_idx <- idx
          break
        }
      }
    }
    fallback <- is.na(plateau_idx)
    start_s <- if (fallback) settings$fallback_creep_start_s else branch_profile$time_since_step_s[[plateau_idx]]
    tibble::tibble(
      branch = branch,
      pressure_plateau_start_s = start_s,
      pressure_plateau_band_bar = settings$pressure_plateau_band_bar,
      pressure_plateau_fallback_used = fallback
    )
  }))
}

predict_creep_model <- function(model, parameters, elapsed_after_start_s) {
  x <- pmax(as.numeric(elapsed_after_start_s), 0)
  if (model == "zero") return(rep(0, length(x)))
  if (model == "linear_to_zero") return(pmax(0, parameters[["a"]] - parameters[["m"]] * x))
  if (model == "single_exponential") return(pmax(0, parameters[["a"]] * exp(-x / parameters[["tau_s"]])))
  if (model == "linear_then_exponential") {
    join_s <- parameters[["join_s"]]
    linear <- pmax(0, parameters[["a"]] - parameters[["m"]] * x)
    at_join <- pmax(0, parameters[["a"]] - parameters[["m"]] * join_s)
    tail <- at_join * exp(-(x - join_s) / parameters[["tau_s"]])
    return(ifelse(x <= join_s, linear, pmax(0, tail)))
  }
  stop(sprintf("Unknown creep model: %s", model))
}

fit_linear_to_zero <- function(profile_df) {
  fit <- stats::lm(normalized_remaining_median ~ elapsed_after_start_s, data = profile_df)
  coeff <- stats::coef(fit)
  list(model = "linear_to_zero", parameters = c(a = max(0, coeff[[1]]), m = max(0, -coeff[[2]]), tau_s = NA_real_, join_s = NA_real_))
}

fit_single_exponential <- function(profile_df) {
  y <- profile_df$normalized_remaining_median
  x <- profile_df$elapsed_after_start_s
  start_a <- max(0.001, y[[1]])
  fit <- tryCatch(
    stats::nls(
      y ~ a * exp(-x / tau_s),
      start = list(a = start_a, tau_s = 0.15),
      algorithm = "port",
      lower = c(a = 0, tau_s = 0.01),
      upper = c(a = 2, tau_s = 5),
      control = stats::nls.control(maxiter = 500, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  coeff <- stats::coef(fit)
  list(model = "single_exponential", parameters = c(a = coeff[["a"]], m = NA_real_, tau_s = coeff[["tau_s"]], join_s = NA_real_))
}

fit_linear_then_exponential <- function(profile_df) {
  x <- profile_df$elapsed_after_start_s
  y <- profile_df$normalized_remaining_median
  if (length(x) < 12) return(NULL)
  join_candidates <- unique(stats::quantile(x, probs = c(0.20, 0.30, 0.40, 0.50, 0.60), names = FALSE))
  # Search relaxation constants over the complete measured 4 s hold horizon.
  tau_candidates <- seq(0.03, 3.99, length.out = 80)
  best <- NULL
  for (join_s in join_candidates) {
    early <- x <= join_s
    if (sum(early) < 3) next
    linear_fit <- stats::lm(y[early] ~ x[early])
    linear_coeff <- stats::coef(linear_fit)
    a <- max(0, linear_coeff[[1]])
    m <- max(0, -linear_coeff[[2]])
    at_join <- max(0, a - m * join_s)
    if (at_join <= 0) next
    for (tau_s in tau_candidates) {
      predicted <- ifelse(x <= join_s, pmax(0, a - m * x), at_join * exp(-(x - join_s) / tau_s))
      rss <- sum((y - predicted)^2)
      if (is.null(best) || rss < best$rss) best <- list(rss = rss, a = a, m = m, tau_s = tau_s, join_s = join_s)
    }
  }
  if (is.null(best)) return(NULL)
  list(model = "linear_then_exponential", parameters = unlist(best[c("a", "m", "tau_s", "join_s")]))
}

fit_model <- function(model, profile_df) {
  if (model == "zero") return(list(model = "zero", parameters = c(a = 0, m = 0, tau_s = NA_real_, join_s = NA_real_)))
  if (model == "linear_to_zero") return(fit_linear_to_zero(profile_df))
  if (model == "single_exponential") return(fit_single_exponential(profile_df))
  if (model == "linear_then_exponential") return(fit_linear_then_exponential(profile_df))
  stop(sprintf("Unknown model %s", model))
}

model_validity <- function(fit, runtime_horizon_s) {
  # Validate the selected profile through the last measured sample of the 4 s hold.
  horizon_s <- max(0, runtime_horizon_s)
  grid <- if (horizon_s <= 0) 0 else seq(0, horizon_s, by = 0.001)
  predicted <- predict_creep_model(fit$model, fit$parameters, grid)
  finite <- all(is.finite(predicted))
  nonnegative <- finite && all(predicted >= -1e-9)
  nonincreasing <- finite && all(diff(predicted) <= 1e-7)
  list(valid = finite && nonnegative && nonincreasing, finite = finite, nonnegative = nonnegative, nonincreasing = nonincreasing)
}

make_folds <- function(hold_metrics, k_folds) {
  hold_metrics %>%
    dplyr::filter(accepted_for_creep) %>%
    dplyr::group_by(branch) %>%
    dplyr::arrange(abs(delta_u_v), .by_group = TRUE) %>%
    dplyr::mutate(cv_fold = rep(seq_len(k_folds), length.out = dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::select(hold_id, cv_fold)
}

profile_from_samples <- function(samples, plateau_df) {
  samples %>%
    dplyr::inner_join(plateau_df %>% dplyr::select(branch, pressure_plateau_start_s), by = "branch") %>%
    dplyr::filter(time_since_step_s >= pressure_plateau_start_s) %>%
    dplyr::mutate(elapsed_after_start_s = time_since_step_s - pressure_plateau_start_s) %>%
    dplyr::group_by(branch, time_since_step_s, elapsed_after_start_s) %>%
    dplyr::summarise(normalized_remaining_median = stats::median(normalized_remaining_displacement), .groups = "drop")
}

cross_validate_models <- function(cleaned_samples, folds, plateau_df, settings) {
  model_names <- c("zero", "linear_to_zero", "single_exponential", "linear_then_exponential")
  branch_rows <- lapply(c("up", "down"), function(branch) {
    branch_samples <- cleaned_samples %>%
      dplyr::filter(accepted_for_creep, branch == !!branch) %>%
      dplyr::inner_join(folds, by = "hold_id") %>%
      dplyr::inner_join(plateau_df %>% dplyr::filter(branch == !!branch), by = "branch") %>%
      dplyr::filter(time_since_step_s >= pressure_plateau_start_s) %>%
      dplyr::mutate(elapsed_after_start_s = time_since_step_s - pressure_plateau_start_s)
    lapply(model_names, function(model) {
      fold_errors <- lapply(seq_len(settings$cv_folds), function(fold) {
        train <- branch_samples %>% dplyr::filter(cv_fold != fold)
        test <- branch_samples %>% dplyr::filter(cv_fold == fold)
        train_profile <- train %>%
          dplyr::group_by(time_since_step_s, elapsed_after_start_s) %>%
          dplyr::summarise(normalized_remaining_median = stats::median(normalized_remaining_displacement), .groups = "drop")
        fit <- fit_model(model, train_profile)
        if (is.null(fit)) return(NA_real_)
        validity <- model_validity(fit, settings$runtime_horizon_s - branch_samples$pressure_plateau_start_s[[1]])
        if (!validity$valid) return(NA_real_)
        prediction <- predict_creep_model(fit$model, fit$parameters, test$elapsed_after_start_s)
        safe_rmse(test$normalized_remaining_displacement - prediction)
      })
      full_profile <- branch_samples %>%
        dplyr::group_by(time_since_step_s, elapsed_after_start_s) %>%
        dplyr::summarise(normalized_remaining_median = stats::median(normalized_remaining_displacement), .groups = "drop")
      full_fit <- fit_model(model, full_profile)
      if (is.null(full_fit)) {
        return(tibble::tibble(branch = branch, model = model, cv_rmse = NA_real_, cv_rmse_sd = NA_real_, cv_folds_completed = 0L, a = NA_real_, m = NA_real_, tau_s = NA_real_, join_s = NA_real_, valid_runtime_profile = FALSE, runtime_profile_finite = FALSE, runtime_profile_nonnegative = FALSE, runtime_profile_nonincreasing = FALSE))
      }
      validity <- model_validity(full_fit, settings$runtime_horizon_s - branch_samples$pressure_plateau_start_s[[1]])
      tibble::tibble(
        branch = branch,
        model = model,
        cv_rmse = mean(unlist(fold_errors), na.rm = TRUE),
        cv_rmse_sd = stats::sd(unlist(fold_errors), na.rm = TRUE),
        cv_folds_completed = sum(is.finite(unlist(fold_errors))),
        a = full_fit$parameters[["a"]],
        m = full_fit$parameters[["m"]],
        tau_s = full_fit$parameters[["tau_s"]],
        join_s = full_fit$parameters[["join_s"]],
        valid_runtime_profile = validity$valid,
        runtime_profile_finite = validity$finite,
        runtime_profile_nonnegative = validity$nonnegative,
        runtime_profile_nonincreasing = validity$nonincreasing
      )
    })
  })
  dplyr::bind_rows(unlist(branch_rows, recursive = FALSE))
}

select_models <- function(model_comparison) {
  runtime_models <- c("zero", "single_exponential", "linear_then_exponential")
  complexity <- c(zero = 0L, linear_to_zero = 2L, single_exponential = 2L, linear_then_exponential = 4L)
  model_comparison %>%
    dplyr::mutate(model_complexity = unname(complexity[model]), eligible_for_runtime = model %in% runtime_models & valid_runtime_profile & is.finite(cv_rmse) & cv_folds_completed >= 3L) %>%
    dplyr::group_by(branch) %>%
    dplyr::mutate(
      best_eligible_cv_rmse = min(cv_rmse[eligible_for_runtime], na.rm = TRUE),
      within_simplicity_margin = eligible_for_runtime & cv_rmse <= 1.05 * best_eligible_cv_rmse,
      selected_model = within_simplicity_margin & model_complexity == min(model_complexity[within_simplicity_margin], na.rm = TRUE)
    ) %>%
    dplyr::arrange(branch, model_complexity, cv_rmse, .by_group = TRUE) %>%
    dplyr::mutate(selected_model = selected_model & dplyr::row_number() == which(selected_model)[1]) %>%
    dplyr::ungroup()
}

build_runtime_profiles <- function(selected_models, plateau_df, settings, calibrated_horizon_s) {
  time_grid <- seq(0, settings$runtime_horizon_s, by = settings$profile_dt_s)
  dplyr::bind_rows(lapply(c("up", "down"), function(branch) {
    model_row <- selected_models %>% dplyr::filter(branch == !!branch, selected_model)
    plateau_row <- plateau_df %>% dplyr::filter(branch == !!branch)
    if (nrow(model_row) != 1 || nrow(plateau_row) != 1) stop(sprintf("Missing selected model or plateau for %s branch.", branch))
    parameters <- c(a = model_row$a[[1]], m = model_row$m[[1]], tau_s = model_row$tau_s[[1]], join_s = model_row$join_s[[1]])
    elapsed <- pmax(0, time_grid - plateau_row$pressure_plateau_start_s[[1]])
    fraction <- predict_creep_model(model_row$model[[1]], parameters, elapsed)
    fraction[time_grid < plateau_row$pressure_plateau_start_s[[1]]] <- 0
    tibble::tibble(
      branch = branch,
      time_since_step_s = time_grid,
      elapsed_after_creep_start_s = elapsed,
      creep_fraction = fraction,
      selected_model = model_row$model[[1]],
      pressure_plateau_start_s = plateau_row$pressure_plateau_start_s[[1]],
      calibrated_horizon_s = calibrated_horizon_s,
      runtime_horizon_s = settings$runtime_horizon_s,
      interval = dplyr::case_when(
        time_grid < plateau_row$pressure_plateau_start_s[[1]] ~ "pre_pressure_plateau_no_correction",
        time_grid <= calibrated_horizon_s + 1e-12 ~ "calibrated_open_loop_hold",
        TRUE ~ "exponential_model_extension"
      )
    )
  }))
}

build_voltage_examples <- function(runtime_profiles, lookup_specs, control_window, settings) {
  dplyr::bind_rows(lapply(seq_len(nrow(lookup_specs)), function(idx) {
    spec <- lookup_specs[idx, ]
    lookup <- read_lookup(spec$path, spec$branch)
    span <- max(lookup$s_ref_mm) - min(lookup$s_ref_mm)
    ref_values <- stats::setNames(
      as.numeric(stats::quantile(lookup$s_ref_mm, probs = c(0.20, 0.50, 0.80), names = FALSE)),
      c("low", "mid", "high")
    )
    step_values <- span * c(0.10, 0.30, 0.50) * ifelse(spec$branch == "up", 1, -1)
    profile <- runtime_profiles %>% dplyr::filter(branch == spec$branch)
    dplyr::bind_rows(lapply(names(ref_values), function(reference_level) {
      s_ref_mm <- ref_values[[reference_level]]
      dplyr::bind_rows(lapply(step_values, function(delta_s_ref_mm) {
        s_virtual_raw <- s_ref_mm + delta_s_ref_mm * profile$creep_fraction
        s_virtual_mm <- pmin(max(lookup$s_ref_mm), pmax(min(lookup$s_ref_mm), s_virtual_raw))
        u_base_v <- interpolate_lookup(lookup, rep(s_ref_mm, nrow(profile)))
        u_creep_v <- interpolate_lookup(lookup, s_virtual_mm)
        tibble::tibble(
          controller = spec$controller,
          branch = spec$branch,
          lookup_file = basename(spec$path),
          reference_level = reference_level,
          s_ref_mm = s_ref_mm,
          delta_s_ref_mm = delta_s_ref_mm,
          time_since_step_s = profile$time_since_step_s,
          interval = profile$interval,
          creep_fraction = profile$creep_fraction,
          s_virtual_mm = s_virtual_mm,
          virtual_reference_clipped = abs(s_virtual_mm - s_virtual_raw) > 1e-12,
          u_ff_base_v = u_base_v,
          u_ff_creep_v = u_creep_v,
          delta_u_creep_v = u_creep_v - u_base_v,
          u_creep_within_control_window = dplyr::between(u_creep_v, control_window$window_low_v[[1]], control_window$window_high_v[[1]])
        )
      }))
    }))
  }))
}

create_plots <- function(profile_df, hold_metrics, model_comparison, runtime_profiles, voltage_examples, plateau_df, settings, plot_theme) {
  displacement <- ggplot2::ggplot(profile_df, ggplot2::aes(time_since_step_s, normalized_remaining_median, color = branch, fill = branch)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = normalized_remaining_q10, ymax = normalized_remaining_q90), alpha = 0.16, colour = NA) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_line(data = runtime_profiles, ggplot2::aes(time_since_step_s, creep_fraction, colour = branch), linewidth = 0.50, linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::geom_vline(data = plateau_df, ggplot2::aes(xintercept = pressure_plateau_start_s, colour = branch), linetype = "dotted", show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = max(profile_df$time_since_step_s), linetype = "dotdash") +
    ggsci::scale_color_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::labs(x = "Time since voltage step [s]", y = "Normalized remaining displacement [-]", colour = "Branch", fill = "Branch") +
    plot_theme
  pressure <- ggplot2::ggplot(profile_df, ggplot2::aes(time_since_step_s, terminal_relative_p_median_bar, colour = branch, fill = branch)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = terminal_relative_p_q10_bar, ymax = terminal_relative_p_q90_bar), alpha = 0.16, colour = NA) +
    ggplot2::geom_hline(yintercept = c(-settings$pressure_plateau_band_bar, settings$pressure_plateau_band_bar), linetype = "dashed", colour = "grey35") +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_vline(data = plateau_df, ggplot2::aes(xintercept = pressure_plateau_start_s, colour = branch), linetype = "dotted", show.legend = FALSE) +
    ggsci::scale_color_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::labs(x = "Time since voltage step [s]", y = expression(p - p[terminal] ~ "[bar]"), colour = "Branch", fill = "Branch") +
    plot_theme
  comparison <- model_comparison %>%
    dplyr::filter(is.finite(cv_rmse)) %>%
    ggplot2::ggplot(ggplot2::aes(model, cv_rmse, fill = branch)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_point(data = model_comparison %>% dplyr::filter(selected_model, is.finite(cv_rmse)), ggplot2::aes(model, cv_rmse), inherit.aes = FALSE, shape = 21, fill = "white", size = 1.8, stroke = 0.25) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::scale_x_discrete(labels = c(zero = "No correction", linear_to_zero = "Linear-to-zero", single_exponential = "Single exponential", linear_then_exponential = "Linear + exponential")) +
    ggplot2::labs(x = "Candidate creep model", y = "Normalized-displacement CV RMSE", fill = "Branch") +
    plot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 9))
  voltage <- voltage_examples %>%
    dplyr::group_by(controller, branch) %>%
    dplyr::filter(abs(delta_s_ref_mm) == max(abs(delta_s_ref_mm))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(reference_level = factor(reference_level, levels = c("low", "mid", "high"))) %>%
    ggplot2::ggplot(ggplot2::aes(time_since_step_s, delta_u_creep_v, colour = controller, linetype = interval, group = interaction(controller, interval))) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(branch ~ reference_level, scales = "free_y", labeller = ggplot2::labeller(branch = c(down = "Downward", up = "Upward"), reference_level = c(low = "Low reference", mid = "Mid reference", high = "High reference"))) +
    ggsci::scale_color_nejm(labels = c(open_loop_lookup = "Raw open-loop lookup")) +
    ggplot2::scale_linetype_discrete(labels = c(pre_pressure_plateau_no_correction = "Before pressure plateau", calibrated_open_loop_hold = "Calibrated hold", exponential_model_extension = "Exponential extension")) +
    ggplot2::labs(x = "Time since reference step [s]", y = expression(Delta * u[ff] ~ "[V]"), colour = "Lookup source", linetype = "Interval") +
    plot_theme
  quality <- hold_metrics %>%
    dplyr::mutate(status = dplyr::if_else(accepted_for_creep, "accepted", dplyr::coalesce(exclusion_reason, "excluded"))) %>%
    dplyr::count(branch, status) %>%
    dplyr::filter(branch %in% c("up", "down")) %>%
    ggplot2::ggplot(ggplot2::aes(status, n, fill = branch)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggsci::scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
    ggplot2::scale_x_discrete(labels = c(accepted = "Accepted", incomplete_hold = "Incomplete hold", small_command_step = "Small command step", small_displacement_step = "Small displacement step", startup_or_zero_command_hold = "Startup or zero-command hold", excluded = "Excluded")) +
    ggplot2::labs(x = "Hold outcome", y = "Number of holds", fill = "Branch") +
    plot_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  list(displacement = displacement, pressure = pressure, comparison = comparison, voltage = voltage, quality = quality)
}

main <- function() {
  script_dir <- find_entry_script_dir()
  source(file.path(script_dir, "lib", "setup.R"), local = FALSE)
  load_analysis_packages()
  context <- initialize_analysis_context()
  settings <- parse_settings(list(
    input_csv = "experiment/excitation_experiment_100hz_15min_4s_hold.csv",
    terminal_tail_fraction = 0.05,
    min_full_hold_samples = 380L,
    full_hold_count_tolerance = 2L,
    min_abs_delta_u_v = 0.35,
    min_abs_displacement_step_mm = 1.0,
    max_terminal_s_sd_mm = 0.60,
    max_terminal_s_slope_mm_s = 1.20,
    max_terminal_p_sd_bar = 0.12,
    max_terminal_p_slope_bar_s = 0.25,
    minimum_creep_start_s = 0.15,
    pressure_plateau_band_bar = 0.05,
    fallback_creep_start_s = 0.50,
    cv_folds = 5L,
    # A 4 s hold sampled at 100 Hz has measured samples from 0 to 3.99 s.
    runtime_horizon_s = 3.99,
    profile_dt_s = 0.01
  ))

  input_path <- file.path(context$repo_dir, settings$input_csv)
  raw_df <- load_excitation(input_path)
  settings$ts_s <- infer_sampling_time(raw_df)
  segmented_df <- segment_holds(raw_df, settings$ts_s)
  hold_metrics <- summarise_holds(segmented_df, settings)
  cleaned_samples <- attach_hold_metrics(segmented_df, hold_metrics)
  profile_df <- profile_by_time(cleaned_samples)
  plateau_df <- find_pressure_plateau(profile_df, settings)
  folds <- make_folds(hold_metrics, settings$cv_folds)
  model_comparison <- cross_validate_models(cleaned_samples, folds, plateau_df, settings)
  model_comparison <- select_models(model_comparison)
  if (sum(model_comparison$selected_model) != 2L) stop("Exactly one creep model must be selected per branch.")

  calibrated_horizon_s <- max(profile_df$time_since_step_s)
  runtime_profiles <- build_runtime_profiles(model_comparison, plateau_df, settings, calibrated_horizon_s)

  control_window <- readr::read_csv(file.path(context$output_dir, "recommended_control_window.csv"), show_col_types = FALSE) %>%
    dplyr::filter(level == "global") %>% dplyr::slice(1)
  if (nrow(control_window) != 1) stop("Missing global control window for voltage-example validation.")
  lookup_specs <- tibble::tribble(
    ~controller, ~branch, ~path,
    "open_loop_lookup", "up", file.path(context$output_dir, "lqr_steady_state_lookup_up.csv"),
    "open_loop_lookup", "down", file.path(context$output_dir, "lqr_steady_state_lookup_down.csv")
  )
  voltage_examples <- build_voltage_examples(runtime_profiles, lookup_specs, control_window, settings)

  summary_df <- model_comparison %>%
    dplyr::filter(selected_model) %>%
    dplyr::left_join(plateau_df, by = "branch") %>%
    dplyr::left_join(
      hold_metrics %>% dplyr::group_by(branch) %>% dplyr::summarise(
        n_holds_total = dplyr::n(),
        n_holds_full = sum(is_full_hold),
        n_holds_accepted = sum(accepted_for_creep),
        n_holds_excluded = sum(!accepted_for_creep),
        .groups = "drop"
      ), by = "branch"
    ) %>%
    dplyr::mutate(
      input_dataset = settings$input_csv,
      sampling_time_s = settings$ts_s,
      full_hold_samples = unique(segmented_df$full_modal_n)[[1]],
      calibrated_horizon_s = calibrated_horizon_s,
      runtime_horizon_s = settings$runtime_horizon_s,
      extension_policy = dplyr::if_else(
        settings$runtime_horizon_s <= calibrated_horizon_s,
        "no_extension_runtime_horizon_within_calibrated_hold",
        "selected_exponential_model_extension_after_calibrated_horizon"
      ),
      correction_interpretation = "repeatable_branchwise_open_loop_relaxation_not_isolated_material_parameter",
      post_correction_policy = "post_correction_is_deferred_until_after_real_system_experiment"
    )

  runtime_contract <- summary_df %>%
    dplyr::transmute(
      branch,
      profile_file = sprintf("open_loop_creep_compensation_%s.csv", branch),
      static_lookup_layer = "open_loop_steady_state_lookup",
      required_runtime_inputs = "s_ref_mm,s_ref_prev_mm,time_since_material_reference_step_s,selected_branch",
      material_step_rule = "reset time when abs(s_ref_mm-s_ref_prev_mm)>0.001",
      formula = "delta_s_ref=s_ref_mm-s_ref_prev_mm; s_virtual=clip(s_ref_mm+delta_s_ref*creep_fraction,loaded_lookup_domain); u_ff_creep=U_branch(s_virtual); delta_u_creep=u_ff_creep-U_branch(s_ref_mm)",
      p_ff_policy = "unchanged",
      branch_policy = "select_existing_branch_before_creep_correction; apply_only_when_step_direction_matches_selected_branch",
      output_policy = "use_existing_voltage_saturation_after_feedforward_and_feedback",
      domain_policy = "reject_or_constrain_out_of_lookup_reference; do_not_extrapolate_static_lookup",
      calibrated_horizon_s,
      runtime_horizon_s,
      extension_policy,
      selected_model = model,
      pressure_plateau_start_s,
      coefficient_a = a,
      coefficient_m = m,
      coefficient_tau_s = tau_s,
      coefficient_join_s = join_s
    )

  readr::write_csv(cleaned_samples, file.path(context$output_dir, "open_loop_creep_cleaned_samples.csv"))
  readr::write_csv(hold_metrics, file.path(context$output_dir, "open_loop_creep_hold_metrics.csv"))
  readr::write_csv(profile_df, file.path(context$output_dir, "open_loop_creep_profile_points.csv"))
  readr::write_csv(model_comparison, file.path(context$output_dir, "open_loop_creep_model_comparison.csv"))
  readr::write_csv(summary_df, file.path(context$output_dir, "open_loop_creep_summary.csv"))
  readr::write_csv(runtime_contract, file.path(context$output_dir, "open_loop_creep_runtime_contract.csv"))
  readr::write_csv(voltage_examples, file.path(context$output_dir, "open_loop_creep_voltage_examples.csv"))
  for (branch in c("up", "down")) {
    readr::write_csv(runtime_profiles %>% dplyr::filter(branch == !!branch), file.path(context$output_dir, sprintf("open_loop_creep_compensation_%s.csv", branch)))
  }

  plot_theme <- create_analysis_plot_theme()
  plots <- create_plots(profile_df, hold_metrics, model_comparison, runtime_profiles, voltage_examples, plateau_df, settings, plot_theme)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_creep_normalized_displacement.png"), plots$displacement, width = 15, height = 8.4375, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_creep_pressure_profile.png"), plots$pressure, width = 15, height = 3.75, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_creep_model_comparison.png"), plots$comparison, width = 15, height = 8.4375, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_creep_voltage_examples.png"), plots$voltage, width = 15, height = 7.5, units = "cm", dpi = 600)
  ggplot2::ggsave(file.path(context$output_dir, "open_loop_creep_hold_quality.png"), plots$quality, width = 15, height = 8.4375, units = "cm", dpi = 600)

  cat("Open-loop creep compensation analysis complete\n")
  cat(sprintf("Input: %s\n", input_path))
  cat(sprintf("Sampling time: %.6f s; modal full hold: %d samples; calibrated horizon: %.3f s\n", settings$ts_s, unique(segmented_df$full_modal_n)[[1]], calibrated_horizon_s))
  print(summary_df %>% dplyr::select(branch, n_holds_accepted, pressure_plateau_start_s, model, cv_rmse, a, m, tau_s, join_s))
  cat(sprintf("Outputs written to: %s\n", context$output_dir))
}

main()
