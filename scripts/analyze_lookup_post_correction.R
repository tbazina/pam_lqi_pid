#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

TAIL_SAMPLES <- 100L
MIN_HOLD_SAMPLES <- 380L
REFERENCE_MATCH_TOLERANCE_MM <- 0.001
WARMUP_REFERENCES_MM <- c(5, 10, 15, 20, 25, 28)

MAX_S_SD_MM <- 0.03
MAX_S_SLOPE_MM_S <- 0.10
MAX_P_SD_BAR <- 0.03
MAX_P_SLOPE_BAR_S <- 0.10
MAX_U_SD_V <- 0.03
MAX_U_SLOPE_V_S <- 0.10

EXPECTED_OUTPUT_ROWS <- c(
  "LQI_no_vel_up" = 194L,
  "LQI_no_vel_down" = 197L,
  "LQI_vel_up" = 186L,
  "LQI_vel_down" = 197L,
  "PID_up" = 186L,
  "PID_down" = 197L
)

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

near_any <- function(x, values, tolerance = REFERENCE_MATCH_TOLERANCE_MM) {
  vapply(x, function(value) any(abs(value - values) <= tolerance), logical(1))
}

safe_sd <- function(x) {
  if (length(x) < 2 || all(!is.finite(x))) {
    return(NA_real_)
  }
  stats::sd(x, na.rm = TRUE)
}

safe_slope <- function(y, dt_ms) {
  keep <- is.finite(y) & is.finite(dt_ms) & dt_ms > 0
  y <- y[keep]
  dt_ms <- dt_ms[keep]
  if (length(y) < 2 || length(unique(y)) < 2) {
    return(0)
  }

  time_s <- c(0, cumsum(dt_ms[-1])) / 1000
  as.numeric(stats::coef(stats::lm(y ~ time_s))[[2]])
}

safe_rmse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else sqrt(mean(x^2))
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = probability, names = FALSE, type = 7))
}

parse_experiment_filename <- function(path) {
  file_name <- basename(path)
  match <- regexec(
    "^(LQI_no_vel|LQI_vel|PID)_lookup_(up|down)_post_correction_experiment\\.csv$",
    file_name
  )
  captures <- regmatches(file_name, match)[[1]]
  if (length(captures) != 3) {
    stop(sprintf("Unsupported post-correction experiment filename: %s", file_name))
  }

  list(
    source_file = file_name,
    controller = captures[[2]],
    filename_branch = captures[[3]]
  )
}

load_experiment <- function(path) {
  raw_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  expected_headers <- c("xref[mm]", "s[mm]", "p[bar]", "u[V]", "t[ms]")
  if (ncol(raw_df) != length(expected_headers) || !identical(names(raw_df), expected_headers)) {
    stop(sprintf(
      "%s must contain the five recorded columns in this order: %s",
      basename(path),
      paste(expected_headers, collapse = ", ")
    ))
  }

  # The recorded headers after xref are shifted. Physical ranges establish the
  # actual order as command voltage, displacement, pressure, and loop period.
  names(raw_df) <- c("xref_mm", "u_cmd_v", "s_meas_mm", "p_meas_bar", "dt_ms")
  experiment_df <- raw_df %>%
    transmute(
      row_id = row_number(),
      xref_mm = as.numeric(xref_mm),
      u_cmd_v = as.numeric(u_cmd_v),
      s_meas_mm = as.numeric(s_meas_mm),
      p_meas_bar = as.numeric(p_meas_bar),
      dt_ms = as.numeric(dt_ms)
    )

  if (any(!is.finite(as.matrix(experiment_df %>% select(-row_id))))) {
    stop(sprintf("%s contains non-finite measurements.", basename(path)))
  }

  range_checks <- c(
    all(dplyr::between(experiment_df$xref_mm, 0, 35)),
    all(dplyr::between(experiment_df$u_cmd_v, -0.5, 10.5)),
    all(dplyr::between(experiment_df$s_meas_mm, -2, 40)),
    all(dplyr::between(experiment_df$p_meas_bar, -1.1, 8.5)),
    all(dplyr::between(experiment_df$dt_ms, 0, 100))
  )
  if (!all(range_checks)) {
    stop(sprintf(
      "%s fails physical-range validation for the xref/u/s/p/dt column interpretation.",
      basename(path)
    ))
  }

  experiment_df %>%
    mutate(
      segment_id = cumsum(row_number() == 1L | xref_mm != lag(xref_mm, default = first(xref_mm)))
    )
}

infer_branch <- function(experiment_df, source_file) {
  segment_refs <- experiment_df %>%
    group_by(segment_id) %>%
    summarise(xref_mm = first(xref_mm), .groups = "drop") %>%
    filter(xref_mm > 0, !near_any(xref_mm, WARMUP_REFERENCES_MM)) %>%
    pull(xref_mm)

  deltas <- diff(segment_refs)
  deltas <- deltas[abs(deltas) > REFERENCE_MATCH_TOLERANCE_MM]
  if (length(deltas) == 0) {
    stop(sprintf("Cannot infer lookup direction from %s.", source_file))
  }

  up_fraction <- mean(deltas > 0)
  down_fraction <- mean(deltas < 0)
  if (up_fraction >= 0.95) {
    return("up")
  }
  if (down_fraction >= 0.95) {
    return("down")
  }

  stop(sprintf(
    "%s is not sufficiently monotonic to infer a lookup branch (up=%.3f, down=%.3f).",
    source_file,
    up_fraction,
    down_fraction
  ))
}

load_original_lookup <- function(output_dir, branch) {
  path <- file.path(output_dir, sprintf("lqr_steady_state_lookup_%s.csv", branch))
  if (!file.exists(path)) {
    stop(sprintf("Missing original lookup table: %s", path))
  }

  lookup_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(lookup_df) <- trimws(names(lookup_df))
  required_cols <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")
  if (!identical(names(lookup_df), required_cols)) {
    stop(sprintf("%s must contain exactly: %s", basename(path), paste(required_cols, collapse = ", ")))
  }

  lookup_df <- lookup_df %>%
    transmute(
      branch = trimws(as.character(branch)),
      s_ref_mm = as.numeric(s_ref_mm),
      u_ff_v = as.numeric(u_ff_v),
      p_ff_bar = as.numeric(p_ff_bar)
    ) %>%
    arrange(s_ref_mm)

  if (nrow(lookup_df) == 0 || any(lookup_df$branch != branch)) {
    stop(sprintf("Original %s lookup is empty or has incorrect branch labels.", branch))
  }
  if (any(!is.finite(as.matrix(lookup_df %>% select(-branch)))) || anyDuplicated(lookup_df$s_ref_mm)) {
    stop(sprintf("Original %s lookup contains invalid or duplicate numeric rows.", branch))
  }

  lookup_df
}

summarize_segments <- function(experiment_df, metadata, inferred_branch, original_lookup) {
  segment_list <- split(experiment_df, experiment_df$segment_id)

  rows <- lapply(segment_list, function(segment_df) {
    n_samples <- nrow(segment_df)
    tail_n_samples <- min(TAIL_SAMPLES, n_samples)
    tail_df <- tail(segment_df, tail_n_samples)
    xref_mm <- first(segment_df$xref_mm)

    nearest_index <- which.min(abs(original_lookup$s_ref_mm - xref_mm))
    matched_s_ref_mm <- original_lookup$s_ref_mm[[nearest_index]]
    match_distance_mm <- abs(matched_s_ref_mm - xref_mm)

    tibble(
      source_file = metadata$source_file,
      controller = metadata$controller,
      filename_branch = metadata$filename_branch,
      inferred_branch = inferred_branch,
      filename_branch_mismatch = metadata$filename_branch != inferred_branch,
      segment_id = first(segment_df$segment_id),
      xref_mm = xref_mm,
      n_samples = n_samples,
      duration_s = sum(segment_df$dt_ms[segment_df$dt_ms > 0]) / 1000,
      tail_n_samples = tail_n_samples,
      median_dt_ms = stats::median(segment_df$dt_ms),
      u_ss_v = stats::median(tail_df$u_cmd_v),
      s_ss_mm = stats::median(tail_df$s_meas_mm),
      p_ss_bar = stats::median(tail_df$p_meas_bar),
      tracking_error_mm = s_ss_mm - xref_mm,
      abs_tracking_error_mm = abs(tracking_error_mm),
      u_tail_sd_v = safe_sd(tail_df$u_cmd_v),
      s_tail_sd_mm = safe_sd(tail_df$s_meas_mm),
      p_tail_sd_bar = safe_sd(tail_df$p_meas_bar),
      u_tail_slope_v_s = safe_slope(tail_df$u_cmd_v, tail_df$dt_ms),
      s_tail_slope_mm_s = safe_slope(tail_df$s_meas_mm, tail_df$dt_ms),
      p_tail_slope_bar_s = safe_slope(tail_df$p_meas_bar, tail_df$dt_ms),
      matched_s_ref_mm = matched_s_ref_mm,
      match_distance_mm = match_distance_mm,
      original_u_ff_v = original_lookup$u_ff_v[[nearest_index]],
      original_p_ff_bar = original_lookup$p_ff_bar[[nearest_index]]
    )
  })

  bind_rows(rows) %>%
    mutate(
      exclusion_reason = case_when(
        xref_mm <= 0 ~ "startup_or_zero_reference",
        inferred_branch == "down" & near_any(xref_mm, WARMUP_REFERENCES_MM) ~ "down_branch_warmup",
        n_samples < MIN_HOLD_SAMPLES ~ "incomplete_hold",
        !is.finite(median_dt_ms) | median_dt_ms <= 0 ~ "invalid_loop_period",
        !is.finite(s_tail_sd_mm) | abs(s_tail_slope_mm_s) > MAX_S_SLOPE_MM_S |
          s_tail_sd_mm > MAX_S_SD_MM ~ "unstable_displacement_tail",
        !is.finite(p_tail_sd_bar) | abs(p_tail_slope_bar_s) > MAX_P_SLOPE_BAR_S |
          p_tail_sd_bar > MAX_P_SD_BAR ~ "unstable_pressure_tail",
        !is.finite(u_tail_sd_v) | abs(u_tail_slope_v_s) > MAX_U_SLOPE_V_S |
          u_tail_sd_v > MAX_U_SD_V ~ "unstable_command_tail",
        match_distance_mm > REFERENCE_MATCH_TOLERANCE_MM ~ "no_original_lookup_match",
        TRUE ~ NA_character_
      ),
      keep_for_corrected_lookup = is.na(exclusion_reason),
      raw_u_correction_v = u_ss_v - original_u_ff_v,
      raw_p_correction_bar = p_ss_bar - original_p_ff_bar
    )
}

build_corrected_lookup <- function(points_df, controller, branch) {
  accepted_df <- points_df %>%
    filter(keep_for_corrected_lookup) %>%
    arrange(matched_s_ref_mm)

  if (nrow(accepted_df) == 0) {
    stop(sprintf("No accepted post-correction points remain for %s %s.", controller, branch))
  }
  if (anyDuplicated(accepted_df$matched_s_ref_mm)) {
    stop(sprintf("%s %s maps multiple experiment holds to one original reference.", controller, branch))
  }

  corrected_df <- accepted_df %>%
    transmute(
      branch = branch,
      s_ref_mm = matched_s_ref_mm,
      u_ff_v = as.numeric(stats::isoreg(matched_s_ref_mm, u_ss_v)$yf),
      p_ff_bar = as.numeric(stats::isoreg(matched_s_ref_mm, p_ss_bar)$yf)
    )

  if (any(diff(corrected_df$s_ref_mm) <= 0) ||
      any(diff(corrected_df$u_ff_v) < -1e-12) ||
      any(diff(corrected_df$p_ff_bar) < -1e-12)) {
    stop(sprintf("Monotonic cleanup failed for %s %s.", controller, branch))
  }
  if (any(!is.finite(as.matrix(corrected_df %>% select(-branch))))) {
    stop(sprintf("Corrected lookup contains non-finite values for %s %s.", controller, branch))
  }
  if (any(!dplyr::between(corrected_df$u_ff_v, 0, 10)) ||
      any(!dplyr::between(corrected_df$p_ff_bar, 0, 8))) {
    stop(sprintf("Corrected lookup exceeds the valve or pressure-sensor range for %s %s.", controller, branch))
  }

  corrected_df
}

build_summary_row <- function(points_df, corrected_df, output_file) {
  accepted_df <- points_df %>% filter(keep_for_corrected_lookup)
  raw_u_violations <- sum(diff(accepted_df$u_ss_v[order(accepted_df$matched_s_ref_mm)]) < -1e-12)
  raw_p_violations <- sum(diff(accepted_df$p_ss_bar[order(accepted_df$matched_s_ref_mm)]) < -1e-12)

  tibble(
    controller = first(points_df$controller),
    branch = first(points_df$inferred_branch),
    source_file = first(points_df$source_file),
    output_file = output_file,
    filename_branch = first(points_df$filename_branch),
    inferred_branch = first(points_df$inferred_branch),
    filename_branch_mismatch = first(points_df$filename_branch_mismatch),
    n_segments = nrow(points_df),
    n_startup_excluded = sum(points_df$exclusion_reason == "startup_or_zero_reference", na.rm = TRUE),
    n_warmup_excluded = sum(points_df$exclusion_reason == "down_branch_warmup", na.rm = TRUE),
    n_incomplete_excluded = sum(points_df$exclusion_reason == "incomplete_hold", na.rm = TRUE),
    n_unstable_excluded = sum(grepl("^unstable_", points_df$exclusion_reason), na.rm = TRUE),
    n_unmatched_excluded = sum(points_df$exclusion_reason == "no_original_lookup_match", na.rm = TRUE),
    n_accepted = nrow(accepted_df),
    output_rows = nrow(corrected_df),
    corrected_s_ref_min_mm = min(corrected_df$s_ref_mm),
    corrected_s_ref_max_mm = max(corrected_df$s_ref_mm),
    tracking_error_bias_mm = mean(accepted_df$tracking_error_mm),
    tracking_error_median_mm = stats::median(accepted_df$tracking_error_mm),
    tracking_error_rmse_mm = safe_rmse(accepted_df$tracking_error_mm),
    tracking_error_p95_abs_mm = safe_quantile(accepted_df$abs_tracking_error_mm, 0.95),
    tracking_error_max_abs_mm = max(accepted_df$abs_tracking_error_mm),
    u_correction_median_v = stats::median(accepted_df$raw_u_correction_v),
    u_correction_rmse_v = safe_rmse(accepted_df$raw_u_correction_v),
    u_correction_max_abs_v = max(abs(accepted_df$raw_u_correction_v)),
    p_correction_median_bar = stats::median(accepted_df$raw_p_correction_bar),
    p_correction_rmse_bar = safe_rmse(accepted_df$raw_p_correction_bar),
    p_correction_max_abs_bar = max(abs(accepted_df$raw_p_correction_bar)),
    raw_u_monotonic_violations = raw_u_violations,
    raw_p_monotonic_violations = raw_p_violations,
    final_u_monotonic_violations = sum(diff(corrected_df$u_ff_v) < -1e-12),
    final_p_monotonic_violations = sum(diff(corrected_df$p_ff_bar) < -1e-12),
    reference_match_tolerance_mm = REFERENCE_MATCH_TOLERANCE_MM,
    min_hold_samples = MIN_HOLD_SAMPLES,
    tail_samples = TAIL_SAMPLES
  )
}

output_file_name <- function(controller, branch) {
  sprintf("%s_lookup_%s_post_corrected.csv", tolower(controller), branch)
}

main <- function() {
  script_dir <- find_entry_script_dir()
  repo_dir <- dirname(script_dir)
  experiment_dir <- file.path(repo_dir, "experiment")
  output_dir <- file.path(repo_dir, "analysis_outputs")

  experiment_paths <- list.files(
    experiment_dir,
    pattern = "^(LQI_no_vel|LQI_vel|PID)_lookup_(up|down)_post_correction_experiment\\.csv$",
    full.names = TRUE
  )
  if (length(experiment_paths) != 6) {
    stop(sprintf("Expected exactly six post-correction experiments; found %d.", length(experiment_paths)))
  }

  original_lookups <- list(
    up = load_original_lookup(output_dir, "up"),
    down = load_original_lookup(output_dir, "down")
  )

  result_list <- lapply(sort(experiment_paths), function(path) {
    metadata <- parse_experiment_filename(path)
    experiment_df <- load_experiment(path)
    inferred_branch <- infer_branch(experiment_df, metadata$source_file)
    original_lookup <- original_lookups[[inferred_branch]]
    points_df <- summarize_segments(experiment_df, metadata, inferred_branch, original_lookup)
    corrected_df <- build_corrected_lookup(points_df, metadata$controller, inferred_branch)
    file_name <- output_file_name(metadata$controller, inferred_branch)

    list(
      key = sprintf("%s_%s", metadata$controller, inferred_branch),
      points = points_df,
      corrected = corrected_df,
      output_file = file_name,
      summary = build_summary_row(points_df, corrected_df, file_name)
    )
  })

  result_keys <- vapply(result_list, `[[`, character(1), "key")
  if (anyDuplicated(result_keys) || !setequal(result_keys, names(EXPECTED_OUTPUT_ROWS))) {
    stop(sprintf(
      "Post-correction experiments do not resolve to six unique expected controller/branch pairs: %s",
      paste(result_keys, collapse = ", ")
    ))
  }

  actual_counts <- setNames(
    vapply(result_list, function(result) nrow(result$corrected), integer(1)),
    result_keys
  )
  count_mismatches <- names(EXPECTED_OUTPUT_ROWS)[
    actual_counts[names(EXPECTED_OUTPUT_ROWS)] != EXPECTED_OUTPUT_ROWS
  ]
  if (length(count_mismatches) > 0) {
    details <- vapply(count_mismatches, function(key) {
      sprintf("%s=%d (expected %d)", key, actual_counts[[key]], EXPECTED_OUTPUT_ROWS[[key]])
    }, character(1))
    stop(sprintf("Unexpected corrected lookup row counts: %s", paste(details, collapse = "; ")))
  }

  diagnostics_df <- bind_rows(lapply(result_list, `[[`, "points")) %>%
    arrange(controller, inferred_branch, segment_id)
  summary_df <- bind_rows(lapply(result_list, `[[`, "summary")) %>%
    arrange(controller, branch)

  for (result in result_list) {
    write_csv(result$corrected, file.path(output_dir, result$output_file))
  }
  write_csv(diagnostics_df, file.path(output_dir, "lookup_post_correction_cleaned_points.csv"))
  write_csv(summary_df, file.path(output_dir, "lookup_post_correction_summary.csv"))

  message("Lookup post-correction analysis complete.")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    message(sprintf(
      "%s %s: %d corrected rows, %.3f--%.3f mm, tracking RMSE %.3f mm%s",
      row$controller,
      row$branch,
      row$output_rows,
      row$corrected_s_ref_min_mm,
      row$corrected_s_ref_max_mm,
      row$tracking_error_rmse_mm,
      if (row$filename_branch_mismatch) " (filename branch mismatch corrected)" else ""
    ))
  }
  message(sprintf("Outputs written to: %s", output_dir))
}

main()
