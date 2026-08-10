#!/usr/bin/env Rscript

# Generate one common reference-position stream for real-system validation of
# controllers using the active raw open-loop steady-state lookup tables.

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

allocate_counts <- function(total_n, target_props) {
  raw_counts <- total_n * target_props
  counts <- floor(raw_counts)
  remainder <- total_n - sum(counts)
  if (remainder > 0) {
    order_idx <- order(raw_counts - counts, decreasing = TRUE)
    counts[order_idx[seq_len(remainder)]] <- counts[order_idx[seq_len(remainder)]] + 1L
  }
  stats::setNames(as.integer(counts), names(target_props))
}

allocate_direction_bin_counts <- function(bin_counts_target) {
  quota_df <- tibble(
    step_bin = names(bin_counts_target),
    total_count = as.integer(bin_counts_target)
  ) %>%
    mutate(up_count = total_count %/% 2L, down_count = total_count %/% 2L)

  for (idx in which(quota_df$total_count %% 2L == 1L)) {
    if (sum(quota_df$up_count) <= sum(quota_df$down_count)) {
      quota_df$up_count[idx] <- quota_df$up_count[idx] + 1L
    } else {
      quota_df$down_count[idx] <- quota_df$down_count[idx] + 1L
    }
  }

  quota_df %>%
    dplyr::select(step_bin, up_count, down_count) %>%
    pivot_longer(c(up_count, down_count), names_to = "direction", values_to = "count") %>%
    mutate(direction = recode(direction, up_count = "up", down_count = "down"))
}

validate_lookup_table <- function(path, controller, expected_branch) {
  required_columns <- c("branch", "s_ref_mm", "u_ff_v", "p_ff_bar")
  if (!file.exists(path)) stop(sprintf("Missing active lookup table: %s", path))

  lookup_df <- read_csv(path, show_col_types = FALSE)
  if (!identical(names(lookup_df), required_columns)) {
    stop(sprintf(
      "Lookup schema mismatch in %s. Expected exactly: %s",
      path,
      paste(required_columns, collapse = ", ")
    ))
  }
  if (nrow(lookup_df) < 2) stop(sprintf("Lookup table has fewer than two rows: %s", path))
  if (!all(lookup_df$branch == expected_branch)) {
    stop(sprintf("Lookup branch values do not match %s: %s", expected_branch, path))
  }
  if (!all(vapply(lookup_df[c("s_ref_mm", "u_ff_v", "p_ff_bar")], function(x) all(is.finite(x)), logical(1)))) {
    stop(sprintf("Lookup table contains non-finite values: %s", path))
  }
  if (anyDuplicated(lookup_df$s_ref_mm) || any(diff(lookup_df$s_ref_mm) <= 0)) {
    stop(sprintf("Lookup references must be unique and strictly ascending: %s", path))
  }

  tibble(
    controller = controller,
    branch = expected_branch,
    source_file = basename(path),
    n_references = nrow(lookup_df),
    s_ref_min_mm = min(lookup_df$s_ref_mm),
    s_ref_max_mm = max(lookup_df$s_ref_mm)
  )
}

plan_transition_count <- function(hold_s, target_duration_min, hard_max_duration_min, min_transitions, max_transitions) {
  target_transitions <- max(1L, floor(target_duration_min * 60 / hold_s) - 1L)
  max_feasible <- max(1L, floor(hard_max_duration_min * 60 / hold_s) - 1L)
  n_transitions <- min(max_transitions, target_transitions)
  used_hard_cap <- FALSE
  planning_warning <- NA_character_

  if (n_transitions < min_transitions) {
    used_hard_cap <- TRUE
    n_transitions <- min(max_transitions, max_feasible)
    if (n_transitions < min_transitions) {
      planning_warning <- "Minimum transition target could not be reached within the hard duration cap."
    }
  }

  list(
    n_transitions = as.integer(n_transitions),
    n_holds = as.integer(n_transitions + 1L),
    total_duration_s = (n_transitions + 1L) * hold_s,
    used_hard_cap = used_hard_cap,
    planning_warning = planning_warning
  )
}

generate_transition_sequence <- function(n_transitions, s_low_mm, s_high_mm, seed, quota_df, step_bin_defs) {
  span_mm <- s_high_mm - s_low_mm

  for (attempt in seq_len(200L)) {
    set.seed(seed + attempt - 1L)
    remaining <- quota_df %>%
      left_join(step_bin_defs, by = "step_bin") %>%
      mutate(
        total_count = count,
        remaining_strata = purrr::map(count, function(n) if (n > 0) seq_len(n) else integer(0))
      )
    current_s <- mean(c(s_low_mm, s_high_mm))
    rows <- vector("list", n_transitions)
    success <- TRUE

    for (i in seq_len(n_transitions)) {
      available_mm <- c(up = s_high_mm - current_s, down = current_s - s_low_mm)
      candidates <- remaining %>%
        mutate(candidate_id = row_number()) %>%
        filter(count > 0) %>%
        mutate(
          available_fraction = available_mm[direction] / span_mm,
          feasible_strata = pmap(
            list(lower_frac, upper_frac, total_count, remaining_strata, available_fraction),
            function(lower, upper, total, strata, available) {
              if (length(strata) == 0 || !is.finite(available)) return(integer(0))
              edges <- seq(lower, upper, length.out = total + 1L)
              strata[edges[strata] <= available + 1e-9]
            }
          ),
          feasible_count = lengths(feasible_strata)
        ) %>%
        filter(feasible_count > 0) %>%
        mutate(weight = feasible_count)

      if (nrow(candidates) == 0) {
        success <- FALSE
        break
      }

      pick <- candidates[sample(seq_len(nrow(candidates)), 1L, prob = candidates$weight), ]
      stratum <- sample(pick$feasible_strata[[1]], 1L)
      edges <- seq(pick$lower_frac, pick$upper_frac, length.out = pick$total_count + 1L)
      lower <- edges[stratum]
      upper <- min(edges[stratum + 1L], pick$available_fraction)
      if (!is.finite(upper) || upper < lower) {
        success <- FALSE
        break
      }

      step_fraction <- runif(1L, lower, upper)
      next_s <- current_s + ifelse(pick$direction == "up", 1, -1) * step_fraction * span_mm
      next_s <- min(max(next_s, s_low_mm), s_high_mm)
      delta_s <- next_s - current_s
      rows[[i]] <- tibble(
        trans_id = i,
        s_ref_mm = next_s,
        delta_s_mm = delta_s,
        delta_s_pct_span = 100 * abs(delta_s) / span_mm,
        direction = pick$direction,
        step_bin = pick$step_bin
      )
      remaining$count[pick$candidate_id] <- remaining$count[pick$candidate_id] - 1L
      remaining$remaining_strata[[pick$candidate_id]] <- setdiff(
        remaining$remaining_strata[[pick$candidate_id]], stratum
      )
      current_s <- next_s
    }

    if (success && all(remaining$count == 0L)) return(bind_rows(rows))
  }

  stop("Could not generate a bounded reference schedule that satisfies the requested quotas.")
}

build_schedule <- function(transitions, hold_s, start_s_mm) {
  bind_rows(
    tibble(
      hold_id = 1L,
      t_start_s = 0,
      t_end_s = hold_s,
      s_ref_mm = start_s_mm,
      delta_s_mm = 0,
      delta_s_pct_span = 0,
      direction = "start",
      step_bin = "start"
    ),
    transitions %>%
      transmute(
        hold_id = trans_id + 1L,
        t_start_s = trans_id * hold_s,
        t_end_s = (trans_id + 1L) * hold_s,
        s_ref_mm,
        delta_s_mm,
        delta_s_pct_span,
        direction,
        step_bin
      )
  )
}

expand_samples <- function(schedule_df, fs_hz) {
  hold_s <- schedule_df$t_end_s[[1]] - schedule_df$t_start_s[[1]]
  samples_per_hold <- as.integer(round(hold_s * fs_hz))
  if (samples_per_hold < 1L || abs(samples_per_hold / fs_hz - hold_s) > 1e-9) {
    stop("hold_s * fs_hz must be an exact positive integer to preserve fixed timing.")
  }

  map_dfr(seq_len(nrow(schedule_df)), function(i) {
    row <- schedule_df[i, ]
    tibble(
      time_s = row$t_start_s + (seq_len(samples_per_hold) - 1L) / fs_hz,
      hold_id = row$hold_id,
      s_ref_mm = row$s_ref_mm
    )
  }) %>%
    mutate(sample_id = row_number()) %>%
    relocate(sample_id)
}

regular_breaks <- function(step, anchor = 0) {
  force(step)
  force(anchor)
  function(x) {
    finite_x <- x[is.finite(x)]
    if (length(finite_x) == 0) return(numeric(0))
    seq(
      anchor + floor((min(finite_x) - anchor) / step) * step,
      anchor + ceiling((max(finite_x) - anchor) / step) * step,
      by = step
    )
  }
}

create_plots <- function(schedule_df, sample_df, step_bin_defs, plot_theme) {
  transitions <- schedule_df %>% filter(hold_id > 1)
  step_bin_actual <- transitions %>%
    count(step_bin, name = "count") %>%
    right_join(step_bin_defs, by = "step_bin") %>%
    mutate(
      count = replace_na(count, 0L),
      prop = count / max(1, sum(count)),
      step_bin = factor(step_bin, levels = step_bin_defs$step_bin),
      step_bin_label = factor(step_bin_label, levels = step_bin_defs$step_bin_label)
    )
  direction_df <- transitions %>%
    mutate(trans_id = row_number(), cum_up = cumsum(direction == "up"), cum_down = cumsum(direction == "down")) %>%
    dplyr::select(trans_id, cum_up, cum_down) %>%
    pivot_longer(c(cum_up, cum_down), names_to = "series", values_to = "count") %>%
    mutate(series = recode(series, cum_up = "Upward transitions", cum_down = "Downward transitions"))

  list(
    time = ggplot(schedule_df, aes(t_start_s / 60, s_ref_mm)) +
      geom_step(color = "#1F1F1F", linewidth = 0.35) +
      scale_x_continuous(breaks = regular_breaks(1)) +
      scale_y_continuous(breaks = regular_breaks(4)) +
      labs(x = "Time [min]", y = "Reference [mm]") +
      plot_theme + theme(legend.position = "none"),
    histogram = ggplot(sample_df, aes(s_ref_mm)) +
      geom_histogram(bins = 40, fill = "#3C5488", color = "white", linewidth = 0.1) +
      scale_x_continuous(breaks = regular_breaks(2)) +
      labs(x = "Reference displacement [mm]", y = "Count") +
      plot_theme + theme(legend.position = "none"),
    signed_steps = ggplot(transitions, aes(delta_s_mm, fill = direction)) +
      geom_histogram(bins = 40, alpha = 0.85, color = "white", linewidth = 0.1) +
      scale_fill_nejm(labels = c(down = "Downward", up = "Upward")) +
      scale_x_continuous(breaks = regular_breaks(2)) +
      labs(x = expression(Delta * "s_ref [mm]"), y = "Count", fill = "Direction") +
      plot_theme,
    step_bins = ggplot(step_bin_actual, aes(step_bin_label, prop)) +
      geom_linerange(aes(ymin = target_low_prop, ymax = target_high_prop), color = "grey40", linewidth = 1.5) +
      geom_point(size = 2.2, color = "#B24745") +
      scale_y_continuous(breaks = regular_breaks(0.1), labels = function(x) sprintf("%.0f%%", 100 * x)) +
      labs(x = expression("Step-size bin and " * "|" * Delta * "s_ref| / " * Delta * "s[ref],span"), y = "Proportion of transitions") +
      plot_theme + theme(legend.position = "none"),
    direction = ggplot(direction_df, aes(trans_id, count, color = series)) +
      geom_line(linewidth = 0.4) +
      scale_color_nejm(name = "Series") +
      scale_x_continuous(breaks = regular_breaks(50)) +
      scale_y_continuous(breaks = regular_breaks(25)) +
      labs(x = "Transition index", y = "Cumulative count", color = "Series") +
      plot_theme
  )
}

save_plots <- function(plots, output_dir) {
  dimensions <- list(time = c(15, 3.75), histogram = c(15, 3.75), signed_steps = c(15, 3.75), step_bins = c(15, 3.75), direction = c(15, 3.75))
  filenames <- c(
    time = "validation_reference_time.png",
    histogram = "validation_reference_histogram.png",
    signed_steps = "validation_reference_signed_step_histogram.png",
    step_bins = "validation_reference_step_bin_proportions.png",
    direction = "validation_reference_direction_balance.png"
  )
  for (name in names(plots)) {
    ggsave(file.path(output_dir, filenames[[name]]), plots[[name]], width = dimensions[[name]][1], height = dimensions[[name]][2], units = "cm", dpi = 600)
  }
}

script_dir <- find_entry_script_dir()
source(file.path(script_dir, "lib", "setup.R"), local = FALSE)
load_analysis_packages()
context <- initialize_analysis_context()
plot_theme <- create_analysis_plot_theme()

defaults <- list(
  hold_s = 2.0,
  fs_hz = 100L,
  target_duration_min = 15.0,
  hard_max_duration_min = 25.0,
  min_transitions = 300L,
  max_transitions = 1000L,
  seed = 12345L
)
settings <- parse_settings(defaults)

if (!is.finite(settings$hold_s) || settings$hold_s <= 0 || !is.finite(settings$fs_hz) || settings$fs_hz <= 0) {
  stop("hold_s and fs_hz must be positive.")
}

lookup_specs <- tribble(
  ~controller, ~branch, ~filename,
  "open_loop_lookup", "up", "lqr_steady_state_lookup_up.csv",
  "open_loop_lookup", "down", "lqr_steady_state_lookup_down.csv"
)
lookup_bounds <- pmap_dfr(lookup_specs, function(controller, branch, filename) {
  validate_lookup_table(file.path(context$output_dir, filename), controller, branch)
})

s_low_mm <- max(lookup_bounds$s_ref_min_mm)
s_high_mm <- min(lookup_bounds$s_ref_max_mm)
if (!is.finite(s_low_mm) || !is.finite(s_high_mm) || s_low_mm >= s_high_mm) {
  stop("The corrected lookup tables do not have a non-empty common reference range.")
}

planning <- plan_transition_count(
  settings$hold_s,
  settings$target_duration_min,
  settings$hard_max_duration_min,
  settings$min_transitions,
  settings$max_transitions
)
step_bin_defs <- tibble(
  step_bin = c("xsmall", "small", "medium", "large", "xlarge"),
  lower_frac = c(0.05, 0.10, 0.25, 0.40, 0.60),
  upper_frac = c(0.10, 0.25, 0.40, 0.60, 0.85),
  target_prop = c(0.10, 0.50, 0.25, 0.12, 0.03),
  target_low_prop = c(0.05, 0.45, 0.20, 0.10, 0.00),
  target_high_prop = c(0.15, 0.55, 0.30, 0.15, 0.05)
) %>%
  mutate(step_bin_label = sprintf("%s\n%.0f-%.0f%% of span", recode(step_bin, xsmall = "Very small", .default = str_to_title(step_bin)), 100 * lower_frac, 100 * upper_frac))

bin_counts <- allocate_counts(planning$n_transitions, setNames(step_bin_defs$target_prop, step_bin_defs$step_bin))
quota_df <- allocate_direction_bin_counts(bin_counts)
transitions <- generate_transition_sequence(
  planning$n_transitions,
  s_low_mm,
  s_high_mm,
  settings$seed,
  quota_df,
  step_bin_defs
)
schedule_df <- build_schedule(transitions, settings$hold_s, mean(c(s_low_mm, s_high_mm)))
sample_df <- expand_samples(schedule_df, settings$fs_hz)

transition_summary <- transitions %>%
  count(step_bin, name = "count") %>%
  right_join(step_bin_defs, by = "step_bin") %>%
  mutate(count = replace_na(count, 0L), prop = count / nrow(transitions), prop_ok = prop >= target_low_prop & prop <= target_high_prop)
up_count <- sum(transitions$direction == "up")
down_count <- sum(transitions$direction == "down")
imbalance_pct <- 100 * abs(up_count - down_count) / nrow(transitions)
samples_per_hold <- as.integer(round(settings$hold_s * settings$fs_hz))
bounds_ok <- all(sample_df$s_ref_mm >= s_low_mm - 1e-9 & sample_df$s_ref_mm <= s_high_mm + 1e-9)
balance_ok <- imbalance_pct <= 10
near_zero_ok <- all(abs(transitions$delta_s_mm) >= 0.05 * (s_high_mm - s_low_mm) - 1e-9)
timing_ok <- nrow(sample_df) == planning$n_holds * samples_per_hold &&
  abs(max(schedule_df$t_end_s) - planning$total_duration_s) <= 1e-9

summary_df <- tibble(
  seed = settings$seed,
  fs_hz = settings$fs_hz,
  hold_s = settings$hold_s,
  samples_per_hold = samples_per_hold,
  target_duration_min = settings$target_duration_min,
  hard_max_duration_min = settings$hard_max_duration_min,
  min_transitions = settings$min_transitions,
  max_transitions = settings$max_transitions,
  planned_transitions = nrow(transitions),
  planned_holds = nrow(schedule_df),
  planned_samples = nrow(sample_df),
  total_duration_min = max(schedule_df$t_end_s) / 60,
  used_hard_cap = planning$used_hard_cap,
  planning_warning = planning$planning_warning,
  reference_range_policy = "intersection_raw_open_loop_steady_state_lookup_up_down",
  lookup_source_count = nrow(lookup_bounds),
  lookup_source_files = paste(lookup_bounds$source_file, collapse = ";"),
  common_reference_low_mm = s_low_mm,
  common_reference_high_mm = s_high_mm,
  common_reference_span_mm = s_high_mm - s_low_mm,
  up_count = up_count,
  down_count = down_count,
  up_down_imbalance_pct = imbalance_pct,
  bounds_ok = bounds_ok,
  balance_ok = balance_ok,
  near_zero_ok = near_zero_ok,
  timing_ok = timing_ok
) %>%
  bind_cols(
    transition_summary %>%
      dplyr::select(step_bin, count, prop, prop_ok) %>%
      pivot_wider(names_from = step_bin, values_from = c(count, prop, prop_ok), names_sep = "_")
  ) %>%
  mutate(all_constraints_ok = bounds_ok && balance_ok && near_zero_ok && timing_ok && all(transition_summary$prop_ok))

write_csv(schedule_df, file.path(context$output_dir, "validation_reference_schedule.csv"))
write_csv(sample_df, file.path(context$output_dir, "validation_reference_samples.csv"))
write_csv(summary_df, file.path(context$output_dir, "validation_reference_signal_summary.csv"))
save_plots(create_plots(schedule_df, sample_df, step_bin_defs, plot_theme), context$output_dir)

cat("\nReference-position validation signal summary\n")
print(summary_df)
cat(sprintf("\nSaved validation-reference CSVs and plots to: %s\n", context$output_dir))
