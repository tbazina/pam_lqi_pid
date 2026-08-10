# Load the packages used across all analysis modules.
load_analysis_packages <- function() {
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
    library(ggsci)
    library(segmented)
    library(tibble)
  })
}

# Locate the repository root robustly for scripted and interactive sessions.
find_repo_dir <- function() {
  repo_markers <- c("experiment", "scripts", "analysis_outputs")

  has_repo_markers <- function(path) {
    if (is.na(path) || !nzchar(path) || !dir.exists(path)) {
      return(FALSE)
    }
    all(file.exists(file.path(path, repo_markers)))
  }

  find_upwards <- function(start_path) {
    if (is.na(start_path) || !nzchar(start_path)) {
      return(NA_character_)
    }

    current <- normalizePath(start_path, winslash = "/", mustWork = FALSE)
    if (file.exists(current) && !dir.exists(current)) {
      current <- dirname(current)
    }

    repeat {
      if (has_repo_markers(current)) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        return(NA_character_)
      }
      current <- parent
    }
  }

  repo_override <- Sys.getenv("LQR_PAM_REPO_DIR", unset = "")
  if (nzchar(repo_override)) {
    repo_dir <- find_upwards(repo_override)
    if (!is.na(repo_dir)) {
      return(repo_dir)
    }
  }

  script_args <- commandArgs(trailingOnly = FALSE)
  script_match <- script_args[grepl("^--file=", script_args)]
  if (length(script_match) > 0) {
    script_path <- sub("^--file=", "", script_match[1])
    repo_dir <- find_upwards(script_path)
    if (!is.na(repo_dir)) {
      return(repo_dir)
    }
  }

  repo_dir <- find_upwards(getwd())
  if (!is.na(repo_dir)) {
    return(repo_dir)
  }

  stop(
    paste(
      "Could not locate the repository root.",
      "Run the script from inside the repository or set LQR_PAM_REPO_DIR."
    )
  )
}

# Build the shared path context used by the analysis entrypoint.
initialize_analysis_context <- function() {
  repo_dir <- find_repo_dir()
  output_dir <- file.path(repo_dir, "analysis_outputs")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  list(
    repo_dir = repo_dir,
    input_dir = file.path(repo_dir, "experiment"),
    input_pattern = "^static_threshold_hysteresis_step_.*\\.csv$",
    output_dir = output_dir
  )
}

# Centralize run-time analysis settings in one place.
get_analysis_settings <- function() {
  list(
    bootstrap_B = 0L,
    bootstrap_seed = 123L
  )
}

# Convert internal run identifiers to compact figure labels.
humanize_run_id <- function(x) {
  x <- as.character(x)
  x <- sub("^step_", "Step ", x)
  x <- sub("V_settle_", " V, hold ", x, fixed = TRUE)
  sub("ms$", " ms", x)
}

# Define the common ggplot theme used for all exported figures.
create_analysis_plot_theme <- function() {
  theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(
        colour = "black",
        size = 9,
        margin = margin(0, 2, 0, 0, "mm")
      ),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.box.spacing = grid::unit(0.5, "mm"),
      legend.spacing.y = grid::unit(0, "mm"),
      legend.spacing.x = grid::unit(0, "mm"),
      legend.margin = margin(0, 0, 0, 0, "mm"),
      legend.key.spacing = grid::unit(0, "mm"),
      legend.key.size = grid::unit(3, "mm"),
      legend.text = element_text(
        colour = "black",
        size = 9,
        margin = margin(0, 0.3, 0, 0, "mm")
      ),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
      panel.background = element_blank(),
      panel.spacing.y = grid::unit(0.5, "mm"),
      panel.spacing.x = grid::unit(0.5, "mm"),
      axis.title = element_text(face = "plain", size = 10),
      axis.text = element_text(
        color = "black",
        size = 9,
        margin = margin(0, 0, 0, 0, "mm")
      ),
      axis.line = element_line(linewidth = 0.2, colour = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid = element_line(colour = "grey80", linewidth = 0.2),
      panel.border = element_rect(linewidth = 0.2),
      strip.background = element_rect(linewidth = 0.01),
      strip.text = element_text(
        colour = "black",
        face = "plain",
        size = 10,
        margin = margin(b = 0.3, t = 0.3, unit = "mm")
      ),
      axis.ticks = element_line(linewidth = 0.2),
      axis.ticks.length = grid::unit(0.1, "lines")
    )
}

# Use slightly smaller text for figures exported at the manuscript's wide size.
create_wide_analysis_plot_theme <- function() {
  create_analysis_plot_theme() +
    theme(
      axis.title = element_text(face = "plain", size = 9),
      axis.text = element_text(colour = "black", size = 8),
      legend.title = element_text(colour = "black", size = 8),
      legend.text = element_text(colour = "black", size = 8),
      strip.text = element_text(face = "plain", colour = "black", size = 9)
    )
}
