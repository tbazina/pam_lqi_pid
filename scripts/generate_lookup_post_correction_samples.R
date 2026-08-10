#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

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

# Read one branchwise lookup table and enforce the expected schema.
load_lookup_table <- function(path, expected_branch) {
  if (!file.exists(path)) {
    stop(sprintf("Missing lookup table: %s", path))
  }

  raw_df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE)
  names(raw_df) <- trimws(names(raw_df))

  required_cols <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")
  missing_cols <- setdiff(required_cols, names(raw_df))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Lookup table %s must contain columns: %s",
      basename(path),
      paste(required_cols, collapse = ", ")
    ))
  }

  lookup_df <- raw_df %>%
    transmute(
      branch = trimws(as.character(branch)),
      s_ref_mm = as.numeric(s_ref_mm),
      u_ff_v = as.numeric(u_ff_v),
      p_ff_bar = as.numeric(p_ff_bar)
    )

  if (nrow(lookup_df) == 0) {
    stop(sprintf("Lookup table %s is empty.", basename(path)))
  }

  if (!all(lookup_df$branch == expected_branch)) {
    stop(sprintf(
      "Lookup table %s contains unexpected branch values. Expected only '%s'.",
      basename(path),
      expected_branch
    ))
  }

  non_finite_rows <- which(!is.finite(lookup_df$s_ref_mm))
  if (length(non_finite_rows) > 0) {
    stop(sprintf(
      "Lookup table %s contains non-finite s_ref_mm values at rows: %s",
      basename(path),
      paste(non_finite_rows, collapse = ", ")
    ))
  }

  lookup_df
}

# Reorder one lookup table so branch playback follows the intended movement direction.
orient_lookup_table <- function(lookup_df, expected_branch) {
  if (expected_branch == "up") {
    ordered_df <- lookup_df %>%
      arrange(s_ref_mm)
  } else if (expected_branch == "down") {
    ordered_df <- lookup_df %>%
      arrange(desc(s_ref_mm))
  } else {
    stop(sprintf("Unsupported branch orientation request: %s", expected_branch))
  }

  ordered_df
}

# Expand one lookup row per hold into an exact 10 ms sampled reference stream.
expand_lookup_to_samples <- function(lookup_df, sample_period_s = 0.01, hold_duration_s = 4.0) {
  samples_per_hold <- as.integer(round(hold_duration_s / sample_period_s))

  if (!isTRUE(all.equal(samples_per_hold * sample_period_s, hold_duration_s, tolerance = 1e-12))) {
    stop("hold_duration_s must be an integer multiple of sample_period_s.")
  }

  hold_template <- tibble(
    sample_offset = seq.int(0L, samples_per_hold - 1L),
    time_offset_s = seq.int(0L, samples_per_hold - 1L) * sample_period_s
  )

  sample_df <- bind_rows(lapply(seq_len(nrow(lookup_df)), function(i) {
    hold_start_s <- (i - 1L) * hold_duration_s
    s_ref_value <- lookup_df$s_ref_mm[[i]]

    tibble(
      time_s = hold_start_s + hold_template$time_offset_s,
      hold_id = i,
      s_ref_mm = s_ref_value
    )
  })) %>%
    mutate(sample_id = row_number()) %>%
    relocate(sample_id)

  sample_df
}

# Emit a compact console summary after writing the LabVIEW-ready outputs.
print_generation_summary <- function(branch_name, lookup_df, sample_df) {
  hold_counts <- sample_df %>%
    count(hold_id, name = "n_samples")

  message(sprintf(
    "%s: %d holds expanded to %d samples (%d samples/hold); first hold %.2f to %.2f s",
    branch_name,
    nrow(lookup_df),
    nrow(sample_df),
    hold_counts$n_samples[[1]],
    min(sample_df$time_s[sample_df$hold_id == 1]),
    max(sample_df$time_s[sample_df$hold_id == 1])
  ))
}

main <- function() {
  script_dir <- find_entry_script_dir()
  repo_dir <- dirname(script_dir)
  output_dir <- file.path(repo_dir, "analysis_outputs")

  up_lookup_path <- file.path(output_dir, "lqr_steady_state_lookup_up.csv")
  down_lookup_path <- file.path(output_dir, "lqr_steady_state_lookup_down.csv")
  up_output_path <- file.path(output_dir, "lqr_steady_state_lookup_up_samples.csv")
  down_output_path <- file.path(output_dir, "lqr_steady_state_lookup_down_samples.csv")

  up_lookup_df <- load_lookup_table(up_lookup_path, expected_branch = "up") %>%
    orient_lookup_table(expected_branch = "up")
  down_lookup_df <- load_lookup_table(down_lookup_path, expected_branch = "down") %>%
    orient_lookup_table(expected_branch = "down")

  up_sample_df <- expand_lookup_to_samples(up_lookup_df)
  down_sample_df <- expand_lookup_to_samples(down_lookup_df)

  write_csv(up_sample_df, up_output_path)
  write_csv(down_sample_df, down_output_path)

  print_generation_summary("up", up_lookup_df, up_sample_df)
  print_generation_summary("down", down_lookup_df, down_sample_df)
  message(sprintf("Wrote: %s", up_output_path))
  message(sprintf("Wrote: %s", down_output_path))
}

main()
