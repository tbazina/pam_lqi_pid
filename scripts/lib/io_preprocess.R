# Extract experiment step size and settle time from the CSV file name.
parse_run_parameters <- function(path) {
  file <- basename(path)
  m <- str_match(
    file,
    "step_(\\d+)p(\\d+)v_settle_time_(\\d+)ms(?:_p_voltage_bar)?\\.csv$"
  )
  if (any(is.na(m))) {
    stop(sprintf("Could not parse parameters from file name: %s", file))
  }

  tibble(
    file = file,
    step_v = as.numeric(sprintf("%s.%s", m[, 2], m[, 3])),
    settle_ms = as.integer(m[, 4])
  )
}

# Infer monotonic up/down segments and attach lagged variables for diagnostics.
infer_direction_and_cycles <- function(df) {
  first_nonzero_sign <- sign(df$u[2] - df$u[1])
  if (first_nonzero_sign == 0) {
    nz <- which(diff(df$u) != 0)
    first_nonzero_sign <- if (length(nz) > 0) sign(diff(df$u)[nz[1]]) else 1
  }
  first_direction <- if (first_nonzero_sign >= 0) "up" else "down"

  direction_raw <- case_when(
    dplyr::lag(df$u, default = df$u[1]) < df$u ~ "up",
    dplyr::lag(df$u, default = df$u[1]) > df$u ~ "down",
    TRUE ~ NA_character_
  )
  direction_raw[1] <- first_direction

  direction <- tidyr::fill(
    tibble(direction = direction_raw),
    direction,
    .direction = "downup"
  )$direction
  segment_start <- c(TRUE, direction[-1] != direction[-nrow(df)])
  segment_id <- cumsum(segment_start)

  df %>%
    mutate(
      direction = factor(direction, levels = c("up", "down")),
      segment_id = segment_id,
      u_prev = lag(u),
      p_prev = lag(p),
      s_prev = lag(s)
    )
}

# Remove short partial edge segments so only complete cycles remain.
prune_incomplete_segments <- function(df, keep_fraction = 0.5) {
  seg_info <- df %>%
    count(segment_id, direction, name = "segment_n") %>%
    arrange(segment_id)

  ref_n <- median(seg_info$segment_n)
  min_keep_n <- max(5, floor(keep_fraction * ref_n))

  seg_map <- seg_info %>%
    filter(segment_n >= min_keep_n) %>%
    mutate(
      segment_id_clean = row_number(),
      cycle_id = ceiling(segment_id_clean / 2)
    )

  df %>%
    left_join(seg_map, by = c("segment_id", "direction")) %>%
    filter(!is.na(segment_id_clean)) %>%
    mutate(segment_id = segment_id_clean) %>%
    dplyr::select(-segment_n, -segment_id_clean)
}

# Load all experiment CSVs and combine them into one annotated tibble.
load_raw_data <- function(input_dir, input_pattern) {
  csv_files <- list.files(input_dir, pattern = input_pattern, full.names = TRUE)
  if (length(csv_files) == 0) {
    stop(sprintf("No matching CSV files found in: %s", input_dir))
  }

  map_dfr(csv_files, function(path) {
    params <- parse_run_parameters(path)
    df <- read_csv(path, show_col_types = FALSE, trim_ws = TRUE) %>%
      rename(u = `u[V]`, p = `p[bar]`, s = `s[mm]`) %>%
      mutate(row_id = row_number())

    bind_cols(df, params[rep(1, nrow(df)), ]) %>%
      infer_direction_and_cycles() %>%
      prune_incomplete_segments()
  }) %>%
    mutate(
      run_id = sprintf("step_%0.2fV_settle_%dms", step_v, settle_ms),
      run_id = factor(run_id, levels = unique(run_id))
    )
}
