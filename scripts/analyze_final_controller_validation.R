#!/usr/bin/env Rscript

# Analyze active LabVIEW/PAM final-validation logs against the 2 s reference
# stream. These are measured closed-loop results, not controller-selection
# simulations.

find_entry_script_dir <- function() {
  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", script_match[1]))))
  }
  candidate <- file.path(getwd(), "scripts")
  if (dir.exists(candidate)) normalizePath(candidate) else normalizePath(getwd())
}

parse_settings <- function(defaults) {
  settings <- defaults
  for (arg in commandArgs(trailingOnly = TRUE)) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) < 2) stop(sprintf("Arguments must use --name=value: %s", arg))
    key <- pieces[1]
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

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else sum(x)
}

read_required_csv <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", label, path))
  read_csv(path, show_col_types = FALSE)
}

read_lookup <- function(path, expected_branch) {
  lookup <- read_required_csv(path, "raw open-loop lookup")
  required <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")
  if (!identical(names(lookup), required)) stop(sprintf("Unexpected lookup schema in %s", path))
  numeric_values <- unlist(lookup[c("s_ref_mm", "u_ff_v", "p_ff_bar")])
  if (!all(lookup$branch == expected_branch) || any(!is.finite(numeric_values))) {
    stop(sprintf("Invalid %s lookup: %s", expected_branch, path))
  }
  if (anyDuplicated(lookup$s_ref_mm) || any(diff(lookup$s_ref_mm) <= 0)) {
    stop(sprintf("Lookup references must be unique and ascending: %s", path))
  }
  lookup
}

read_runtime_profile <- function(path, value_column, expected_branch) {
  profile <- read_required_csv(path, sprintf("%s runtime profile", expected_branch))
  required <- c("branch", "time_since_step_s", value_column)
  if (!all(required %in% names(profile))) stop(sprintf("Unexpected profile schema in %s", path))
  profile <- profile %>% dplyr::select(all_of(required)) %>% arrange(time_since_step_s)
  if (!all(profile$branch == expected_branch) || any(!is.finite(profile$time_since_step_s)) ||
      any(!is.finite(profile[[value_column]])) || any(diff(profile$time_since_step_s) <= 0)) {
    stop(sprintf("Invalid %s runtime profile: %s", expected_branch, path))
  }
  profile
}

profile_value <- function(profile, value_column, elapsed_s, terminal_value) {
  horizon_s <- max(profile$time_since_step_s)
  if (elapsed_s > horizon_s + 1e-12) return(terminal_value)
  approx(profile$time_since_step_s, profile[[value_column]], xout = elapsed_s, method = "linear", rule = 2)$y
}

read_reference_inputs <- function(output_dir, settings) {
  samples <- read_required_csv(file.path(output_dir, settings$reference_samples_file), "validation reference samples")
  schedule <- read_required_csv(file.path(output_dir, settings$reference_schedule_file), "validation reference schedule")
  sample_required <- c("sample_id", "time_s", "hold_id", "s_ref_mm")
  schedule_required <- c("hold_id", "t_start_s", "t_end_s", "s_ref_mm", "delta_s_mm", "delta_s_pct_span", "direction", "step_bin")
  if (!identical(names(samples), sample_required)) stop("validation_reference_samples.csv has an unexpected schema.")
  if (!identical(names(schedule), schedule_required)) stop("validation_reference_schedule.csv has an unexpected schema.")
  if (any(!is.finite(unlist(samples))) || any(diff(samples$time_s) <= 0)) {
    stop("Reference samples must be finite with strictly increasing time.")
  }
  hold_reference <- samples %>%
    group_by(hold_id, s_ref_mm) %>%
    summarise(expected_samples = n(), .groups = "drop") %>%
    left_join(schedule, by = c("hold_id", "s_ref_mm")) %>%
    arrange(hold_id)
  if (nrow(hold_reference) != nrow(schedule) || any(is.na(hold_reference$direction)) ||
      any(hold_reference$expected_samples != settings$samples_per_hold)) {
    stop("Reference samples and schedule are not a complete fixed-hold validation stream.")
  }
  samples <- samples %>%
    left_join(schedule %>% dplyr::select(hold_id, direction, step_bin, delta_s_mm, delta_s_pct_span), by = "hold_id")
  list(samples = samples, schedule = schedule, holds = hold_reference)
}

read_validation_log <- function(path) {
  raw <- read_required_csv(path, "LabVIEW validation log")
  names(raw) <- trimws(names(raw))
  required <- c("xref[mm]", "s[mm]", "p[bar]", "u[V]", "t[ms]")
  if (!identical(names(raw), required)) {
    stop(sprintf("Unexpected validation-log schema in %s. Expected: %s", path, paste(required, collapse = ", ")))
  }
  clean <- raw %>% transmute(
    source_row_id = row_number(),
    xref_logged_mm = as.numeric(`xref[mm]`),
    s_meas_mm = as.numeric(`s[mm]`),
    p_meas_bar = as.numeric(`p[bar]`),
    u_cmd_v = as.numeric(`u[V]`),
    loop_dt_ms = as.numeric(`t[ms]`)
  )
  if (any(!is.finite(unlist(clean[-1])))) stop(sprintf("Non-finite validation values in %s", path))
  clean
}

segment_logged_holds <- function(log_df, settings) {
  positive_rows <- which(log_df$xref_logged_mm > settings$alignment_tolerance_mm)
  if (length(positive_rows) == 0) stop("Validation log has no positive scheduled references.")
  first_scheduled <- min(positive_rows)
  last_scheduled <- max(positive_rows)
  if (first_scheduled != 2L || any(log_df$xref_logged_mm[seq_len(first_scheduled - 1L)] > settings$alignment_tolerance_mm)) {
    stop("Expected exactly one non-schedule startup row before the first positive reference.")
  }
  if (any(log_df$xref_logged_mm[(last_scheduled + 1L):nrow(log_df)] > settings$alignment_tolerance_mm)) {
    stop("Unexpected positive-reference row after the scheduled validation stream.")
  }
  scheduled <- log_df %>% slice(first_scheduled:last_scheduled)
  run_change <- c(TRUE, abs(diff(scheduled$xref_logged_mm)) > settings$alignment_tolerance_mm)
  scheduled$logged_run_id <- cumsum(run_change)
  runs <- scheduled %>%
    group_by(logged_run_id) %>%
    summarise(
      logged_xref_mm = first(xref_logged_mm),
      observed_samples = n(),
      source_row_start = min(source_row_id),
      source_row_end = max(source_row_id),
      .groups = "drop"
    )
  if (any(runs$observed_samples != settings$samples_per_hold)) {
    stop("Every retained logged reference plateau must have the expected complete-hold sample count.")
  }
  list(
    samples = scheduled,
    runs = runs,
    startup_rows_dropped = first_scheduled - 1L,
    terminal_rows_dropped = nrow(log_df) - last_scheduled
  )
}

align_logged_holds <- function(logged_runs, reference_holds, settings) {
  expected_xref <- round(reference_holds$s_ref_mm, 3)
  next_expected <- 1L
  matched_rows <- vector("list", nrow(logged_runs))
  missing_rows <- list()
  for (idx in seq_len(nrow(logged_runs))) {
    candidates <- which(abs(expected_xref[next_expected:length(expected_xref)] - logged_runs$logged_xref_mm[idx]) <= settings$alignment_tolerance_mm)
    if (length(candidates) == 0) {
      stop(sprintf("Cannot monotonically align logged hold %d at %.3f mm.", idx, logged_runs$logged_xref_mm[idx]))
    }
    expected_idx <- next_expected + candidates[1] - 1L
    if (expected_idx > next_expected) {
      missing_rows[[length(missing_rows) + 1L]] <- reference_holds[next_expected:(expected_idx - 1L), ] %>%
        mutate(alignment_status = "missing_logged_hold", logged_run_id = NA_integer_, logged_xref_mm = NA_real_, match_distance_mm = NA_real_)
    }
    matched_rows[[idx]] <- bind_cols(
      logged_runs[idx, ],
      reference_holds[expected_idx, ] %>% dplyr::select(-expected_samples)
    ) %>% mutate(alignment_status = "matched", match_distance_mm = abs(logged_xref_mm - round(s_ref_mm, 3)))
    next_expected <- expected_idx + 1L
  }
  if (next_expected <= nrow(reference_holds)) {
    missing_rows[[length(missing_rows) + 1L]] <- reference_holds[next_expected:nrow(reference_holds), ] %>%
      mutate(alignment_status = "missing_logged_hold", logged_run_id = NA_integer_, logged_xref_mm = NA_real_, match_distance_mm = NA_real_)
  }
  audit <- bind_rows(bind_rows(matched_rows), bind_rows(missing_rows)) %>% arrange(hold_id)
  if (any(audit$alignment_status == "matched" & audit$match_distance_mm > settings$alignment_tolerance_mm)) {
    stop("Reference hold alignment exceeded tolerance.")
  }
  audit
}

build_aligned_samples <- function(logged_samples, alignment_audit, reference_samples) {
  matched <- alignment_audit %>% filter(alignment_status == "matched") %>% arrange(logged_run_id)
  sample_blocks <- map(seq_len(nrow(matched)), function(idx) {
    run <- matched[idx, ]
    observed <- logged_samples %>% filter(logged_run_id == run$logged_run_id) %>% arrange(source_row_id)
    canonical <- reference_samples %>% filter(hold_id == run$hold_id) %>% arrange(sample_id)
    if (nrow(observed) != nrow(canonical)) stop("Aligned run and canonical reference block have different lengths.")
    bind_cols(observed %>% dplyr::select(-logged_run_id), canonical, run %>% dplyr::select(logged_run_id, source_row_start, source_row_end, match_distance_mm))
  })
  bind_rows(sample_blocks) %>% arrange(hold_id, sample_id)
}

first_order_lowpass <- function(x, ts_s, tau_s) {
  if (tau_s <= 0 || length(x) == 0) return(x)
  alpha <- exp(-ts_s / tau_s)
  y <- numeric(length(x))
  y[1] <- x[1]
  for (idx in 2:length(x)) y[idx] <- alpha * y[idx - 1L] + (1 - alpha) * x[idx]
  y
}

sliding_linear_slope <- function(s_mm, ts_s, window_samples) {
  n <- length(s_mm)
  slopes <- numeric(n)
  if (n < 2L) return(slopes)
  window_samples <- min(as.integer(window_samples), n)
  time_index <- seq_len(window_samples) - 1L
  centered_time <- time_index - mean(time_index)
  denominator <- sum(centered_time^2)
  # stats::filter is causal for sides = 1; reverse the oldest-to-newest slope weights.
  slopes <- as.numeric(stats::filter(s_mm, filter = rev(centered_time / denominator / ts_s), sides = 1L))
  slopes[is.na(slopes)] <- 0.0
  warmup_end <- min(window_samples - 1L, n)
  if (warmup_end >= 2L) {
    for (idx in 2:warmup_end) {
      segment <- s_mm[seq_len(idx)]
      prefix_time <- (seq_len(idx) - 1L) * ts_s
      centered_prefix_time <- prefix_time - mean(prefix_time)
      prefix_denominator <- sum(centered_prefix_time^2)
      slopes[idx] <- sum(centered_prefix_time * (segment - mean(segment))) / prefix_denominator
    }
  }
  slopes
}

lookup_linear <- function(lookup, column, reference_mm) {
  if (reference_mm < min(lookup$s_ref_mm) - 1e-9 || reference_mm > max(lookup$s_ref_mm) + 1e-9) {
    stop(sprintf("Reference %.6f mm is outside the raw %s lookup domain.", reference_mm, lookup$branch[1]))
  }
  approx(lookup$s_ref_mm, lookup[[column]], xout = reference_mm, method = "linear", rule = 1)$y
}

reconstruct_runtime_layers <- function(sample_df, up_lookup, down_lookup, creep_profiles, pressure_profiles, config, settings) {
  n <- nrow(sample_df)
  reference <- sample_df$s_ref_mm
  measured <- sample_df$s_meas_mm
  hold_ids <- sample_df$hold_id
  active_code <- integer(n)
  current_branch <- 1L
  for (idx in seq_len(n)) {
    if (reference[idx] > measured[idx] + config$branch_deadband_mm) {
      current_branch <- 1L
    } else if (reference[idx] < measured[idx] - config$branch_deadband_mm) {
      current_branch <- 2L
    }
    active_code[idx] <- current_branch
  }
  active_branch <- if_else(active_code == 1L, "up", "down")
  u_static_up <- approx(up_lookup$s_ref_mm, up_lookup$u_ff_v, xout = reference, method = "linear", rule = 1)$y
  u_static_down <- approx(down_lookup$s_ref_mm, down_lookup$u_ff_v, xout = reference, method = "linear", rule = 1)$y
  p_static_up <- approx(up_lookup$s_ref_mm, up_lookup$p_ff_bar, xout = reference, method = "linear", rule = 1)$y
  p_static_down <- approx(down_lookup$s_ref_mm, down_lookup$p_ff_bar, xout = reference, method = "linear", rule = 1)$y
  u_static <- if_else(active_code == 1L, u_static_up, u_static_down)
  p_static <- if_else(active_code == 1L, p_static_up, p_static_down)

  hold_start <- !duplicated(hold_ids)
  hold_index <- match(hold_ids, unique(hold_ids))
  hold_references <- reference[hold_start]
  hold_delta <- c(0.0, diff(hold_references))
  hold_branch <- if_else(hold_delta > 0, "up", if_else(hold_delta < 0, "down", "start"))
  elapsed_s <- (ave(seq_len(n), hold_ids, FUN = seq_along) - 1L) * settings$ts_s
  step_delta <- hold_delta[hold_index]
  step_branch <- hold_branch[hold_index]
  branch_at_hold_start <- active_branch[hold_start]
  direction_matches_hold <- hold_branch == branch_at_hold_start & abs(hold_delta) > settings$alignment_tolerance_mm
  direction_matches <- direction_matches_hold[hold_index]

  creep_fraction_up <- approx(creep_profiles$up$time_since_step_s, creep_profiles$up$creep_fraction, xout = elapsed_s, method = "linear", rule = 2)$y
  creep_fraction_down <- approx(creep_profiles$down$time_since_step_s, creep_profiles$down$creep_fraction, xout = elapsed_s, method = "linear", rule = 2)$y
  creep_fraction <- if_else(active_code == 1L, creep_fraction_up, creep_fraction_down)
  creep_fraction[!config$uses_creep_compensation | !direction_matches] <- 0.0
  virtual_reference <- reference + step_delta * creep_fraction
  virtual_reference[active_code == 1L] <- pmin(pmax(virtual_reference[active_code == 1L], min(up_lookup$s_ref_mm)), max(up_lookup$s_ref_mm))
  virtual_reference[active_code == 2L] <- pmin(pmax(virtual_reference[active_code == 2L], min(down_lookup$s_ref_mm)), max(down_lookup$s_ref_mm))
  u_creep_up <- approx(up_lookup$s_ref_mm, up_lookup$u_ff_v, xout = virtual_reference, method = "linear", rule = 1)$y
  u_creep_down <- approx(down_lookup$s_ref_mm, down_lookup$u_ff_v, xout = virtual_reference, method = "linear", rule = 1)$y
  u_creep <- if_else(active_code == 1L, u_creep_up, u_creep_down)

  dynamic_active_hold <- config$uses_dynamic_pressure_reference & direction_matches_hold
  p_step_hold <- sample_df$p_meas_bar[hold_start]
  p_static_hold <- p_static[hold_start]
  pressure_progress_up <- approx(pressure_profiles$up$time_since_step_s, pressure_profiles$up$pressure_progress_fraction, xout = elapsed_s, method = "linear", rule = 2)$y
  pressure_progress_down <- approx(pressure_profiles$down$time_since_step_s, pressure_profiles$down$pressure_progress_fraction, xout = elapsed_s, method = "linear", rule = 2)$y
  pressure_progress <- if_else(branch_at_hold_start[hold_index] == "up", pressure_progress_up, pressure_progress_down)
  pressure_progress[!dynamic_active_hold[hold_index]] <- 1.0
  p_dynamic <- p_static
  active_dynamic <- dynamic_active_hold[hold_index]
  p_dynamic[active_dynamic] <- p_step_hold[hold_index][active_dynamic] + (p_static_hold[hold_index][active_dynamic] - p_step_hold[hold_index][active_dynamic]) * pressure_progress[active_dynamic]

  tibble(
    active_branch = active_branch,
    time_since_reference_step_s = elapsed_s,
    reference_step_delta_mm = step_delta,
    reference_step_branch = step_branch,
    u_ff_static_v = u_static,
    u_ff_v = u_creep,
    creep_fraction = creep_fraction,
    p_ff_static_bar = p_static,
    p_ff_dynamic_bar = p_dynamic,
    dynamic_pressure_progress = pressure_progress,
    p_tilde_static_bar = sample_df$p_meas_bar - p_static,
    p_tilde_dynamic_bar = sample_df$p_meas_bar - p_dynamic
  )
}

classify_step_bin <- function(reference_delta_mm, reference_span_mm) {
  percentage <- 100 * abs(reference_delta_mm) / reference_span_mm
  if (!is.finite(percentage) || percentage < 5) return("below_minimum")
  if (percentage < 10) return("xsmall")
  if (percentage < 25) return("small")
  if (percentage < 40) return("medium")
  if (percentage < 60) return("large")
  if (percentage <= 85) return("xlarge")
  "above_design_range"
}

compute_hold_metrics <- function(sample_df, settings, reference_span_mm) {
  hold_ids <- unique(sample_df$hold_id)
  rows <- vector("list", length(hold_ids))
  previous_end_s <- NA_real_
  previous_ref <- NA_real_
  for (idx in seq_along(hold_ids)) {
    hold_df <- sample_df %>% filter(hold_id == hold_ids[idx]) %>% arrange(sample_id)
    is_initial <- idx == 1L
    tail_n <- max(5L, ceiling(settings$settled_window_fraction * nrow(hold_df)))
    error <- hold_df$error_mm
    tail_error <- tail(error, tail_n)
    tail_command <- tail(hold_df$u_cmd_v, tail_n)
    tail_velocity_common <- tail(hold_df$velocity_common_mm_s, tail_n)
    tail_velocity_native <- tail(hold_df$velocity_native_mm_s, tail_n)
    ref <- hold_df$s_ref_mm[1]
    initial_s <- if (is_initial) NA_real_ else previous_end_s
    reference_delta <- if (is_initial) NA_real_ else ref - previous_ref
    step_mm <- if (is_initial) NA_real_ else ref - initial_s
    direction <- if (is_initial) "start" else if (reference_delta >= 0) "up" else "down"
    step_bin <- if (is_initial) "start" else classify_step_bin(reference_delta, reference_span_mm)
    transition_eligible <- !is_initial
    step_eligible <- transition_eligible && abs(step_mm) >= settings$step_metric_min_mm
    local_time <- (seq_len(nrow(hold_df)) - 1L) * settings$ts_s
    row <- tibble(
      controller = hold_df$controller[1],
      hold_id = hold_ids[idx],
      direction = direction,
      step_bin = step_bin,
      scheduled_direction = hold_df$direction[1],
      scheduled_step_bin = hold_df$step_bin[1],
      n_samples = nrow(hold_df),
      complete_hold = TRUE,
      is_initial_hold = is_initial,
      transition_metric_eligible = transition_eligible,
      step_metric_eligible = step_eligible,
      s_ref_mm = ref,
      initial_s_mm = initial_s,
      reference_delta_mm = reference_delta,
      step_mm = step_mm,
      settled_error_rms_mm = sqrt(mean(tail_error^2)),
      settled_error_p2p_mm = max(tail_error) - min(tail_error),
      settled_command_p2p_v = max(tail_command) - min(tail_command),
      settled_velocity_rms_mm_s = sqrt(mean(tail_velocity_common^2)),
      settled_velocity_rms_native_mm_s = if (all(is.finite(tail_velocity_native))) sqrt(mean(tail_velocity_native^2)) else NA_real_,
      steady_state_bias_mm = mean(tail_error),
      abs_steady_state_bias_mm = abs(mean(tail_error)),
      IAE_mm_s = sum(abs(error)) * settings$ts_s,
      ISE_mm2_s = sum(error^2) * settings$ts_s,
      ITAE_mm_s2 = sum(local_time * abs(error)) * settings$ts_s,
      control_effort_mean_abs_v = mean(abs(hold_df$u_cmd_v - hold_df$u_ff_v)),
      control_residual_rms_v = sqrt(mean((hold_df$u_cmd_v - hold_df$u_ff_v)^2)),
      command_total_variation_v = if (nrow(hold_df) > 1) sum(abs(diff(hold_df$u_cmd_v))) else 0.0,
      saturation_fraction = mean(hold_df$command_saturated),
      rise_time_s_10_90 = NA_real_,
      settling_time_s_2pct = NA_real_,
      peak_time_s = NA_real_,
      overshoot_mm = NA_real_,
      overshoot_pct = NA_real_
    )
    if (step_eligible) {
      step_sign <- ifelse(step_mm >= 0, 1.0, -1.0)
      progress <- step_sign * (hold_df$s_meas_mm - initial_s) / abs(step_mm)
      rise_low <- which(progress >= settings$rise_low)[1]
      rise_high <- which(progress >= settings$rise_high)[1]
      settling_band <- max(settings$settling_band_fraction * abs(step_mm), settings$settling_band_floor_mm)
      deviation <- abs(hold_df$s_meas_mm - ref)
      settling_candidates <- which(vapply(seq_along(deviation), function(position) all(deviation[position:length(deviation)] <= settling_band), logical(1)))
      directional_response <- step_sign * (hold_df$s_meas_mm - initial_s)
      peak_idx <- which.max(directional_response)
      row$rise_time_s_10_90 <- if (is.na(rise_low) || is.na(rise_high) || rise_high < rise_low) NA_real_ else (rise_high - rise_low) * settings$ts_s
      row$settling_time_s_2pct <- if (length(settling_candidates) == 0) NA_real_ else (settling_candidates[1] - 1L) * settings$ts_s
      row$peak_time_s <- (peak_idx - 1L) * settings$ts_s
      row$overshoot_mm <- max(step_sign * (hold_df$s_meas_mm - ref), 0.0)
      row$overshoot_pct <- 100 * row$overshoot_mm / abs(step_mm)
    }
    rows[[idx]] <- row
    previous_end_s <- tail(hold_df$s_meas_mm, 1)
    previous_ref <- ref
  }
  bind_rows(rows)
}

summarize_controller <- function(sample_df, hold_metrics, quality_row) {
  transition_metrics <- hold_metrics %>% filter(transition_metric_eligible)
  step_metrics <- hold_metrics %>% filter(step_metric_eligible)
  tibble(
    controller = quality_row$controller,
    analysis_scope = "available_complete_aligned_holds",
    expected_holds = quality_row$expected_holds,
    complete_holds = quality_row$complete_holds,
    missing_holds = quality_row$missing_holds,
    missing_hold_ids = quality_row$missing_hold_ids,
    completed_duration_min = quality_row$completed_duration_min,
    full_reference_run_completed = quality_row$full_reference_run_completed,
    n_tracking_samples = nrow(sample_df),
    n_transition_holds = nrow(transition_metrics),
    n_step_metric_holds = nrow(step_metrics),
    n_small_step_holds_excluded = nrow(transition_metrics) - nrow(step_metrics),
    tracking_rms_mm = sqrt(mean(sample_df$error_mm^2)),
    tracking_overshoot_mm = max(abs(sample_df$error_mm)),
    control_effort_mean_abs_v = mean(abs(sample_df$u_cmd_v - sample_df$u_ff_v)),
    settled_error_rms_mm = safe_mean(transition_metrics$settled_error_rms_mm),
    settled_error_p2p_mm = safe_mean(transition_metrics$settled_error_p2p_mm),
    settled_command_p2p_v = safe_mean(transition_metrics$settled_command_p2p_v),
    settled_velocity_rms_mm_s = safe_mean(transition_metrics$settled_velocity_rms_mm_s),
    settled_velocity_rms_native_mm_s = safe_mean(transition_metrics$settled_velocity_rms_native_mm_s),
    rise_time_s_10_90 = safe_mean(step_metrics$rise_time_s_10_90),
    settling_time_s_2pct = safe_mean(step_metrics$settling_time_s_2pct),
    peak_time_s = safe_mean(step_metrics$peak_time_s),
    overshoot_mm = safe_mean(step_metrics$overshoot_mm),
    overshoot_pct = safe_mean(step_metrics$overshoot_pct),
    steady_state_bias_mm = safe_mean(transition_metrics$steady_state_bias_mm),
    abs_steady_state_bias_mm = safe_mean(transition_metrics$abs_steady_state_bias_mm),
    IAE_mm_s = safe_sum(transition_metrics$IAE_mm_s),
    ISE_mm2_s = safe_sum(transition_metrics$ISE_mm2_s),
    ITAE_mm_s2 = safe_sum(transition_metrics$ITAE_mm_s2),
    control_residual_rms_v = sqrt(mean((sample_df$u_cmd_v - sample_df$u_ff_v)^2)),
    command_total_variation_v = sum(abs(diff(sample_df$u_cmd_v))),
    saturation_fraction = mean(sample_df$command_saturated)
  )
}

extract_lqi_gains <- function(path, variant, terms, require_variant) {
  gains <- read_required_csv(path, "LQI feedback gains")
  if (require_variant) gains <- gains %>% filter(variant == !!variant)
  required <- if (require_variant) c("matrix", "row", "col", "value") else c("matrix", "row", "col", "value")
  if (!all(required %in% names(gains))) stop(sprintf("Unexpected LQI gain schema in %s", path))
  gains %>%
    filter(row == 1L, matrix %in% c("K_aug_up", "K_aug_down")) %>%
    mutate(branch = sub("^K_aug_", "", matrix), term = terms[col]) %>%
    dplyr::select(branch, term, simulated_value = value) %>%
    arrange(branch, match(term, terms))
}

read_deployed_gain_file <- function(path, controller, terms) {
  raw <- read_csv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(raw) != 2L || nrow(raw) != length(terms)) {
    stop(sprintf("%s must contain exactly %d rows and two branch columns.", path, length(terms)))
  }
  values <- as.matrix(raw)
  storage.mode(values) <- "numeric"
  if (any(!is.finite(values))) stop(sprintf("Non-finite deployed gains in %s", path))
  map_dfr(seq_along(terms), function(index) {
    tibble(controller = controller, term = terms[index], branch = c("up", "down"), deployed_value = values[index, ])
  })
}

read_gain_audit <- function(output_dir, input_dir) {
  base_simulated <- extract_lqi_gains(file.path(output_dir, "lqr_feedback_gains.csv"), "", c("K_s", "K_p", "K_i"), FALSE) %>% mutate(controller = "LQI_no_vel")
  vel_simulated <- extract_lqi_gains(file.path(output_dir, "lqr_feedback_gains_nonselected.csv"), "lookup_plus_branchwise_minimal_LQI_VEL", c("K_s", "K_p", "K_v", "K_i"), TRUE) %>% mutate(controller = "LQI_vel")
  pid_simulated <- read_required_csv(file.path(output_dir, "pid_feedback_gains.csv"), "PID feedback gains") %>%
    transmute(controller = "PID", branch, term = recode(term, kp_v_per_mm = "K_p", ki_v_per_mm_s = "K_i", kd_v_s_per_mm = "K_d"), simulated_value = value)
  deployed <- bind_rows(
    read_deployed_gain_file(file.path(input_dir, "PID_validation_experiment_gains.csv"), "PID", c("K_p", "K_i", "K_d")),
    read_deployed_gain_file(file.path(input_dir, "LQI_no_vel_validation_experiment_gains.csv"), "LQI_no_vel", c("K_s", "K_p", "K_i")),
    read_deployed_gain_file(file.path(input_dir, "LQI_vel_validation_experiment_gains.csv"), "LQI_vel", c("K_s", "K_p", "K_v", "K_i"))
  )
  audit <- deployed %>%
    left_join(bind_rows(base_simulated, vel_simulated, pid_simulated), by = c("controller", "term", "branch")) %>%
    mutate(
      absolute_change = deployed_value - simulated_value,
      relative_change_pct = 100 * absolute_change / abs(simulated_value),
      change_status = if_else(abs(absolute_change) <= 1e-9, "unchanged", "manually_adjusted")
    ) %>% arrange(match(controller, c("PID", "LQI_no_vel", "LQI_vel")), match(term, c("K_s", "K_p", "K_v", "K_i", "K_d")), branch)
  if (any(!is.finite(audit$simulated_value))) stop("Could not match every deployed gain to its simulation export.")
  bind_rows(
    audit,
    tibble(controller = "Feedforward_only", term = "not_applicable", branch = NA_character_, simulated_value = NA_real_, deployed_value = NA_real_, absolute_change = NA_real_, relative_change_pct = NA_real_, change_status = "no_feedback_gains")
  )
}

read_controller_configs <- function(output_dir, settings) {
  lqr_summary <- read_required_csv(file.path(output_dir, "lqr_controller_design_summary.csv"), "LQI design summary")
  pid_summary <- read_required_csv(file.path(output_dir, "pid_controller_design_summary.csv"), "PID design summary")
  feedforward_summary <- read_required_csv(file.path(output_dir, "feedforward_only_controller_summary.csv"), "feedforward-only summary")
  base_row <- lqr_summary %>% filter(variant == "lookup_plus_branchwise_minimal_LQI") %>% slice(1)
  vel_row <- lqr_summary %>% filter(variant == "lookup_plus_branchwise_minimal_LQI_VEL") %>% slice(1)
  pid_row <- pid_summary %>% filter(chosen_variant) %>% slice(1)
  feedforward_row <- feedforward_summary %>% slice(1)
  if (nrow(base_row) != 1 || nrow(vel_row) != 1 || nrow(pid_row) != 1 || nrow(feedforward_row) != 1) {
    stop("Could not resolve active controller configurations.")
  }
  tribble(
    ~controller, ~experiment_file, ~branch_deadband_mm, ~native_slope_window_samples, ~native_position_prefilter_tau_s, ~native_velocity_filter_tau_s, ~uses_creep_compensation, ~uses_dynamic_pressure_reference,
    "Feedforward_only", "feedforward_only_validation_experiment.csv", as.numeric(feedforward_row$branch_deadband_mm), NA_integer_, 0.0, 0.0, TRUE, FALSE,
    "PID", "PID_validation_experiment.csv", as.numeric(pid_row$branch_deadband_mm), as.integer(pid_row$slope_window_samples), as.numeric(pid_row$position_prefilter_tau_s), as.numeric(pid_row$tau_d_s), TRUE, FALSE,
    "LQI_no_vel", "LQI_no_vel_validation_experiment.csv", as.numeric(base_row$branch_deadband_mm), NA_integer_, 0.0, 0.0, TRUE, TRUE,
    "LQI_vel", "LQI_vel_validation_experiment.csv", as.numeric(vel_row$branch_deadband_mm), as.integer(vel_row$velocity_slope_window_samples), as.numeric(vel_row$velocity_position_prefilter_tau_s), as.numeric(vel_row$velocity_filter_tau_s), TRUE, TRUE
  ) %>% mutate(u_min_v = settings$deployed_u_min_v, u_max_v = settings$deployed_u_max_v)
}

read_nominal_simulation_metrics <- function(output_dir) {
  lqr <- read_required_csv(file.path(output_dir, "lqr_controller_design_summary.csv"), "LQI design summary") %>%
    transmute(controller = recode(variant, lookup_plus_branchwise_minimal_LQI = "LQI_no_vel", lookup_plus_branchwise_minimal_LQI_VEL = "LQI_vel"), across(where(is.numeric)))
  pid <- read_required_csv(file.path(output_dir, "pid_controller_design_summary.csv"), "PID design summary") %>%
    filter(chosen_variant) %>% transmute(controller = "PID", across(where(is.numeric)))
  feedforward <- read_required_csv(file.path(output_dir, "feedforward_only_controller_summary.csv"), "feedforward-only summary") %>%
    transmute(controller = "Feedforward_only", across(where(is.numeric)))
  bind_rows(feedforward, pid, lqr)
}

metric_units <- function(metric) {
  recode(metric,
    tracking_rms_mm = "mm", tracking_overshoot_mm = "mm", settled_error_rms_mm = "mm", settled_error_p2p_mm = "mm",
    settled_command_p2p_v = "V", settled_velocity_rms_mm_s = "mm/s", rise_time_s_10_90 = "s", settling_time_s_2pct = "s",
    peak_time_s = "s", overshoot_mm = "mm", overshoot_pct = "%", steady_state_bias_mm = "mm", abs_steady_state_bias_mm = "mm",
    IAE_mm_s = "mm s", ISE_mm2_s = "mm2 s", ITAE_mm_s2 = "mm s2", control_residual_rms_v = "V",
    command_total_variation_v = "V", saturation_fraction = "fraction", .default = "")
}

build_simulation_comparison <- function(summary_df, nominal_df) {
  metrics <- c("tracking_rms_mm", "tracking_overshoot_mm", "settled_error_rms_mm", "settled_error_p2p_mm", "settled_command_p2p_v", "settled_velocity_rms_mm_s", "rise_time_s_10_90", "settling_time_s_2pct", "peak_time_s", "overshoot_mm", "overshoot_pct", "steady_state_bias_mm", "abs_steady_state_bias_mm", "IAE_mm_s", "ISE_mm2_s", "ITAE_mm_s2", "control_residual_rms_v", "command_total_variation_v", "saturation_fraction")
  measured_long <- summary_df %>% dplyr::select(controller, all_of(metrics)) %>% pivot_longer(-controller, names_to = "metric", values_to = "measured_value")
  nominal_long <- nominal_df %>% dplyr::select(controller, any_of(metrics)) %>% pivot_longer(-controller, names_to = "metric", values_to = "nominal_simulated_value")
  measured_long %>% left_join(nominal_long, by = c("controller", "metric")) %>%
    mutate(
      unit = metric_units(metric),
      absolute_difference = measured_value - nominal_simulated_value,
      comparison_basis = "measured_deployed_manual_gains_vs_nominal_original_gain_single_dataset_simulation",
      interpretation = "context_only_not_matched_gain_validation"
    ) %>% arrange(match(controller, c("Feedforward_only", "PID", "LQI_no_vel", "LQI_vel")), metric)
}

create_plots <- function(cleaned_df, hold_metrics, summary_df, settings, plot_theme) {
  display_controller <- function(x) {
    factor(recode(x,
      Feedforward_only = "Feedforward",
      PID = "PID",
      LQI_no_vel = "Base LQI",
      LQI_vel = "LQI-VEL"
    ), levels = c("Feedforward", "PID", "Base LQI", "LQI-VEL"))
  }
  overview_df <- cleaned_df %>%
    filter(sample_id %% 10L == 1L) %>% dplyr::select(controller, time_s, s_ref_mm, s_meas_mm) %>%
    pivot_longer(c(s_ref_mm, s_meas_mm), names_to = "signal", values_to = "value_mm") %>%
    mutate(controller = display_controller(controller), signal = recode(signal, s_ref_mm = "Reference", s_meas_mm = "Measured displacement"))
  zoom_df <- cleaned_df %>%
    filter(time_s <= settings$zoom_s) %>% dplyr::select(controller, time_s, s_ref_mm, s_meas_mm) %>%
    pivot_longer(c(s_ref_mm, s_meas_mm), names_to = "signal", values_to = "value_mm") %>%
    mutate(controller = display_controller(controller), signal = recode(signal, s_ref_mm = "Reference", s_meas_mm = "Measured displacement"))
  metric_df <- summary_df %>%
    dplyr::select(controller, tracking_rms_mm, settled_error_rms_mm, settled_error_p2p_mm, rise_time_s_10_90, settling_time_s_2pct, overshoot_pct, saturation_fraction) %>% mutate(controller = display_controller(controller)) %>%
    pivot_longer(-controller, names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric, tracking_rms_mm = "Tracking RMS [mm]", settled_error_rms_mm = "Settled RMS [mm]", settled_error_p2p_mm = "Settled p2p [mm]", rise_time_s_10_90 = "Rise time [s]", settling_time_s_2pct = "Settling time [s]", overshoot_pct = "Overshoot [%]", saturation_fraction = "Saturation fraction"))
  step_df <- hold_metrics %>%
    filter(step_metric_eligible) %>% dplyr::select(controller, direction, rise_time_s_10_90, settling_time_s_2pct, overshoot_pct) %>% mutate(controller = display_controller(controller)) %>%
    pivot_longer(-c(controller, direction), names_to = "metric", values_to = "value") %>% filter(is.finite(value)) %>%
    mutate(metric = recode(metric, rise_time_s_10_90 = "Rise time [s]", settling_time_s_2pct = "Settling time [s]", overshoot_pct = "Overshoot [%]"))
  wide_plot_theme <- create_wide_analysis_plot_theme()
  list(
    overview = ggplot(overview_df, aes(time_s / 60, value_mm, color = signal)) + geom_step(linewidth = 0.36) + facet_wrap(~controller, ncol = 1) + scale_color_nejm() + labs(x = "Time [min]", y = "Displacement [mm]", color = NULL) + wide_plot_theme,
    zoom = ggplot(zoom_df, aes(time_s, value_mm, color = signal)) + geom_step(linewidth = 0.45) + facet_wrap(~controller, ncol = 1) + scale_color_nejm() + labs(x = "Time [s]", y = "Displacement [mm]", color = NULL) + wide_plot_theme,
    error = ggplot(cleaned_df %>% filter(sample_id %% 5L == 1L) %>% mutate(controller = display_controller(controller)), aes(controller, error_mm, fill = controller)) + geom_boxplot(outlier.shape = NA, linewidth = 0.25) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + scale_fill_nejm() + labs(x = "Controller", y = expression(s - s[ref] ~ "[mm]")) + plot_theme + theme(legend.position = "none"),
    comparison = ggplot(metric_df, aes(controller, value, fill = controller)) + geom_col(width = 0.7) + facet_wrap(~metric, scales = "free_y", ncol = 4) + scale_fill_nejm() + scale_x_discrete(labels = c(`Feedforward` = "FF-only", `PID` = "PID", `Base LQI` = "Base LQI", `LQI-VEL` = "LQI-VEL")) + labs(x = "Controller", y = NULL) + wide_plot_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1, size = 8)),
    step_response = ggplot(step_df, aes(controller, value, fill = direction)) + geom_boxplot(outlier.shape = NA, linewidth = 0.25, position = position_dodge(width = 0.75)) + scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) + scale_x_discrete(labels = c(`Feedforward` = "FF-only", `PID` = "PID", `Base LQI` = "Base LQI", `LQI-VEL` = "LQI-VEL")) + facet_wrap(~metric, scales = "free_y", ncol = 3) + labs(x = "Controller", y = NULL, fill = "Direction") + plot_theme + theme(axis.text.x = element_text(size = 8)),
    command = ggplot(cleaned_df %>% filter(sample_id %% 10L == 1L) %>% mutate(controller = display_controller(controller)), aes(time_s / 60, u_cmd_v, color = command_saturated)) + geom_step(linewidth = 0.36) + facet_wrap(~controller, ncol = 1) + scale_color_manual(values = c(`FALSE` = "#3C5488", `TRUE` = "#B24745"), labels = c(`FALSE` = "Within deployed limits", `TRUE` = "At deployed limit")) + labs(x = "Time [min]", y = "Valve command [V]", color = "Saturation status") + wide_plot_theme
  )
}

save_plots <- function(plots, output_dir) {
  filenames <- c(overview = "final_validation_tracking_overview.png", zoom = "final_validation_tracking_zoom.png", error = "final_validation_error_distribution.png", comparison = "final_validation_metric_comparison.png", step_response = "final_validation_step_response_by_direction.png", command = "final_validation_command_saturation.png")
  for (name in names(plots)) {
    dimensions <- switch(name,
      overview = c(17, 9.5625),
      zoom = c(17, 9.5625),
      comparison = c(17, 8.5),
      step_response = c(15, 5),
      command = c(17, 9.5625),
      error = c(15, 8.4375)
    )
    ggsave(file.path(output_dir, filenames[[name]]), plots[[name]], width = dimensions[[1]], height = dimensions[[2]], units = "cm", dpi = 600)
  }
}

script_dir <- find_entry_script_dir()
source(file.path(script_dir, "lib", "setup.R"), local = FALSE)
load_analysis_packages()
context <- initialize_analysis_context()
plot_theme <- create_analysis_plot_theme()

defaults <- list(
  reference_samples_file = "validation_reference_samples.csv",
  reference_schedule_file = "validation_reference_schedule.csv",
  alignment_tolerance_mm = 0.0005,
  saturation_tolerance_v = 0.001,
  deployed_u_min_v = 0.9479660643688704,
  deployed_u_max_v = 7.2,
  settled_window_fraction = 0.30,
  settling_band_fraction = 0.02,
  settling_band_floor_mm = 0.15,
  step_metric_min_mm = 0.50,
  rise_low = 0.10,
  rise_high = 0.90,
  standardized_velocity_window_samples = 5L,
  samples_per_hold = 200L,
  zoom_s = 60.0,
  ts_s = 0.01
)
settings <- parse_settings(defaults)

reference <- read_reference_inputs(context$output_dir, settings)
controller_configs <- read_controller_configs(context$output_dir, settings)
gain_audit_df <- read_gain_audit(context$output_dir, context$input_dir)
up_lookup <- read_lookup(file.path(context$output_dir, "lqr_steady_state_lookup_up.csv"), "up")
down_lookup <- read_lookup(file.path(context$output_dir, "lqr_steady_state_lookup_down.csv"), "down")
creep_profiles <- list(
  up = read_runtime_profile(file.path(context$output_dir, "open_loop_creep_compensation_up.csv"), "creep_fraction", "up"),
  down = read_runtime_profile(file.path(context$output_dir, "open_loop_creep_compensation_down.csv"), "creep_fraction", "down")
)
pressure_profiles <- list(
  up = read_runtime_profile(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_up.csv"), "pressure_progress_fraction", "up"),
  down = read_runtime_profile(file.path(context$output_dir, "open_loop_dynamic_pressure_reference_down.csv"), "pressure_progress_fraction", "down")
)
reference_span_mm <- max(reference$samples$s_ref_mm) - min(reference$samples$s_ref_mm)

all_cleaned <- list()
all_hold_metrics <- list()
quality_rows <- list()
alignment_rows <- list()
for (idx in seq_len(nrow(controller_configs))) {
  config <- controller_configs[idx, ]
  log_df <- read_validation_log(file.path(context$input_dir, config$experiment_file))
  segmented <- segment_logged_holds(log_df, settings)
  alignment_audit <- align_logged_holds(segmented$runs, reference$holds, settings) %>% mutate(controller = config$controller)
  expected_missing <- if (config$controller == "Feedforward_only") 446L else integer()
  missing_ids <- alignment_audit %>% filter(alignment_status == "missing_logged_hold") %>% pull(hold_id)
  if (!identical(as.integer(missing_ids), expected_missing)) {
    stop(sprintf("Unexpected missing canonical holds for %s: %s", config$controller, paste(missing_ids, collapse = ",")))
  }
  sample_df <- build_aligned_samples(segmented$samples, alignment_audit, reference$samples) %>%
    mutate(
      controller = config$controller,
      error_mm = s_meas_mm - s_ref_mm,
      command_saturated = u_cmd_v <= config$u_min_v + settings$saturation_tolerance_v | u_cmd_v >= config$u_max_v - settings$saturation_tolerance_v
    )
  runtime_df <- reconstruct_runtime_layers(sample_df, up_lookup, down_lookup, creep_profiles, pressure_profiles, config, settings)
  common_velocity <- sliding_linear_slope(sample_df$s_meas_mm, settings$ts_s, settings$standardized_velocity_window_samples)
  native_velocity <- rep(NA_real_, nrow(sample_df))
  if (!is.na(config$native_slope_window_samples)) {
    native_position <- first_order_lowpass(sample_df$s_meas_mm, settings$ts_s, config$native_position_prefilter_tau_s)
    native_velocity <- first_order_lowpass(sliding_linear_slope(native_position, settings$ts_s, config$native_slope_window_samples), settings$ts_s, config$native_velocity_filter_tau_s)
  }
  sample_df <- bind_cols(sample_df, runtime_df) %>% mutate(velocity_common_mm_s = common_velocity, velocity_native_mm_s = native_velocity)
  hold_metrics <- compute_hold_metrics(sample_df, settings, reference_span_mm)
  quality_rows[[config$controller]] <- tibble(
    controller = config$controller,
    source_file = config$experiment_file,
    raw_rows = nrow(log_df),
    startup_rows_dropped = segmented$startup_rows_dropped,
    terminal_rows_dropped = segmented$terminal_rows_dropped,
    expected_holds = nrow(reference$holds),
    complete_holds = sum(alignment_audit$alignment_status == "matched"),
    missing_holds = sum(alignment_audit$alignment_status == "missing_logged_hold"),
    missing_hold_ids = if (length(missing_ids) == 0) "" else paste(missing_ids, collapse = ";"),
    completed_duration_min = sum(alignment_audit$alignment_status == "matched") * settings$samples_per_hold * settings$ts_s / 60,
    full_reference_run_completed = length(missing_ids) == 0,
    reference_match_fraction = 1.0,
    schedule_coverage_fraction = mean(alignment_audit$alignment_status == "matched"),
    pressure_negative_samples = sum(sample_df$p_meas_bar < 0),
    command_lower_limit_samples = sum(sample_df$u_cmd_v <= config$u_min_v + settings$saturation_tolerance_v),
    command_upper_limit_samples = sum(sample_df$u_cmd_v >= config$u_max_v - settings$saturation_tolerance_v),
    u_min_v = config$u_min_v,
    u_max_v = config$u_max_v,
    branch_deadband_mm = config$branch_deadband_mm,
    native_slope_window_samples = config$native_slope_window_samples,
    native_position_prefilter_tau_s = config$native_position_prefilter_tau_s,
    native_velocity_filter_tau_s = config$native_velocity_filter_tau_s,
    uses_creep_compensation = config$uses_creep_compensation,
    uses_dynamic_pressure_reference = config$uses_dynamic_pressure_reference
  )
  all_cleaned[[config$controller]] <- sample_df
  all_hold_metrics[[config$controller]] <- hold_metrics
  alignment_rows[[config$controller]] <- alignment_audit
}

cleaned_df <- bind_rows(all_cleaned)
hold_metrics_df <- bind_rows(all_hold_metrics)
quality_df <- bind_rows(quality_rows)
alignment_df <- bind_rows(alignment_rows)
summary_df <- bind_rows(map(quality_rows, function(row) summarize_controller(cleaned_df %>% filter(controller == row$controller), hold_metrics_df %>% filter(controller == row$controller), row))) %>%
  left_join(quality_df, by = c("controller", "expected_holds", "complete_holds", "missing_holds", "missing_hold_ids", "completed_duration_min", "full_reference_run_completed"))
step_summary_df <- hold_metrics_df %>%
  filter(step_metric_eligible) %>%
  group_by(controller, direction, step_bin) %>%
  summarise(n_holds = n(), rise_time_s_10_90 = safe_mean(rise_time_s_10_90), settling_time_s_2pct = safe_mean(settling_time_s_2pct), peak_time_s = safe_mean(peak_time_s), overshoot_mm = safe_mean(overshoot_mm), overshoot_pct = safe_mean(overshoot_pct), settled_error_rms_mm = safe_mean(settled_error_rms_mm), settled_error_p2p_mm = safe_mean(settled_error_p2p_mm), .groups = "drop")
simulation_comparison_df <- build_simulation_comparison(summary_df, read_nominal_simulation_metrics(context$output_dir))

write_csv(cleaned_df, file.path(context$output_dir, "final_validation_cleaned_samples.csv"))
write_csv(hold_metrics_df, file.path(context$output_dir, "final_validation_hold_metrics.csv"))
write_csv(summary_df, file.path(context$output_dir, "final_validation_summary.csv"))
write_csv(quality_df, file.path(context$output_dir, "final_validation_data_quality.csv"))
write_csv(alignment_df, file.path(context$output_dir, "final_validation_alignment_audit.csv"))
write_csv(step_summary_df, file.path(context$output_dir, "final_validation_step_metrics_by_direction.csv"))
write_csv(gain_audit_df, file.path(context$output_dir, "final_validation_deployed_gains.csv"))
write_csv(simulation_comparison_df, file.path(context$output_dir, "final_validation_simulation_comparison.csv"))
save_plots(create_plots(cleaned_df, hold_metrics_df, summary_df, settings, plot_theme), context$output_dir)

cat("\nMeasured final-controller validation summary\n")
print(summary_df)
cat(sprintf("\nSaved validation analysis outputs to: %s\n", context$output_dir))
