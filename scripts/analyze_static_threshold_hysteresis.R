#!/usr/bin/env Rscript

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

# Source one helper file from the local analysis library.
source_analysis_helper <- function(script_dir, helper_file) {
  source(file.path(script_dir, "lib", helper_file), local = FALSE)
}

# Load helper modules in dependency order, then initialize shared context.
script_dir <- find_entry_script_dir()

source_analysis_helper(script_dir, "setup.R")
load_analysis_packages()
source_analysis_helper(script_dir, "io_preprocess.R")
source_analysis_helper(script_dir, "cutoff_model.R")
source_analysis_helper(script_dir, "summaries_validation.R")
source_analysis_helper(script_dir, "plots_outputs.R")

analysis_context <- initialize_analysis_context()
analysis_settings <- get_analysis_settings()
analysis_plot_theme <- create_analysis_plot_theme()

repo_dir <- analysis_context$repo_dir
input_dir <- analysis_context$input_dir
input_pattern <- "^static_threshold_hysteresis_step_.*_p_voltage_bar\\.csv$"
output_dir <- analysis_context$output_dir

bootstrap_B <- analysis_settings$bootstrap_B
bootstrap_seed <- analysis_settings$bootstrap_seed

# Run the analysis pipeline from raw data load through fitted cutoff outputs.
raw_data <- load_raw_data(input_dir, input_pattern)
run_summary <- create_run_summary(raw_data)
cycle_summary <- create_cycle_summary(raw_data)

direction_outputs <- estimate_direction_results(
  raw_data = raw_data,
  bootstrap_B = bootstrap_B,
  bootstrap_seed = bootstrap_seed
)

cutoff_summary <- direction_outputs$cutoff_summary
cutoff_curves <- direction_outputs$cutoff_curves
cutoff_mean_curves <- direction_outputs$cutoff_mean_curves

run_overlap_summary <- create_run_overlap_summary(cutoff_summary)
recommended_control_window <- create_recommended_control_window(run_overlap_summary)
settling_summary <- create_settling_summary(raw_data, recommended_control_window)
cutoff_validation <- create_cutoff_validation(cutoff_summary)
displacement_validation_flags <- create_displacement_validation_flags(cutoff_summary)

# Persist tabular outputs before creating the full plot set.
write_analysis_outputs(
  output_dir = output_dir,
  run_summary = run_summary,
  cycle_summary = cycle_summary,
  cutoff_summary = cutoff_summary,
  settling_summary = settling_summary,
  cutoff_validation = cutoff_validation,
  recommended_control_window = recommended_control_window
)

# Build the raw-data, fitted-response, and diagnostic plot objects.
pressure_plot <- create_pressure_plot(
  raw_data = raw_data,
  cutoff_curves = cutoff_curves,
  analysis_plot_theme = analysis_plot_theme
)
raw_pressure_voltage_scatter_plot <- create_raw_pressure_voltage_scatter_plot(
  raw_data = raw_data,
  analysis_plot_theme = analysis_plot_theme
)
raw_pressure_stroke_scatter_plot <- create_raw_pressure_stroke_scatter_plot(
  raw_data = raw_data,
  analysis_plot_theme = analysis_plot_theme
)
raw_voltage_stroke_scatter_plot <- create_raw_voltage_stroke_scatter_plot(
  raw_data = raw_data,
  analysis_plot_theme = analysis_plot_theme
)
stroke_plot <- create_stroke_plot(
  raw_data = raw_data,
  cutoff_mean_curves = cutoff_mean_curves,
  analysis_plot_theme = analysis_plot_theme
)
usable_region_plot <- create_usable_region_plot(
  raw_data = raw_data,
  cutoff_curves = cutoff_curves,
  cutoff_summary = cutoff_summary,
  analysis_plot_theme = analysis_plot_theme
)
step_distribution_plot <- create_step_distribution_plot(
  raw_data = raw_data,
  analysis_plot_theme = analysis_plot_theme
)

# Save the full plot set and print a compact console summary.
save_analysis_plots(
  output_dir = output_dir,
  raw_pressure_voltage_scatter_plot = raw_pressure_voltage_scatter_plot,
  raw_pressure_stroke_scatter_plot = raw_pressure_stroke_scatter_plot,
  raw_voltage_stroke_scatter_plot = raw_voltage_stroke_scatter_plot,
  pressure_plot = pressure_plot,
  stroke_plot = stroke_plot,
  usable_region_plot = usable_region_plot,
  step_distribution_plot = step_distribution_plot
)

print_analysis_summary(
  output_dir = output_dir,
  run_summary = run_summary,
  cycle_summary = cycle_summary,
  cutoff_summary = cutoff_summary,
  run_overlap_summary = run_overlap_summary,
  recommended_control_window = recommended_control_window,
  displacement_validation_flags = displacement_validation_flags,
  cutoff_validation = cutoff_validation,
  settling_summary = settling_summary
)
