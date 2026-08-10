# Write the tabular outputs used for inspection and later reporting.
write_analysis_outputs <- function(
  output_dir,
  run_summary,
  cycle_summary,
  cutoff_summary,
  settling_summary,
  cutoff_validation,
  recommended_control_window
) {
  write_csv(run_summary, file.path(output_dir, "run_summary.csv"))
  write_csv(cycle_summary, file.path(output_dir, "cycle_summary.csv"))
  write_csv(cutoff_summary, file.path(output_dir, "usable_region_cutoffs.csv"))
  write_csv(settling_summary, file.path(output_dir, "settling_assessment.csv"))
  write_csv(cutoff_validation, file.path(output_dir, "cutoff_validation.csv"))
  write_csv(
    recommended_control_window,
    file.path(output_dir, "recommended_control_window.csv")
  )
}

# Use a small horizontal jitter so repeated voltages remain visible in scatter plots.
get_point_jitter_width <- function(raw_data, fraction = 0.2, fallback = 0.002) {
  u_step <- raw_data %>%
    distinct(u) %>%
    arrange(u) %>%
    mutate(delta_u = u - lag(u)) %>%
    summarise(min_delta_u = suppressWarnings(min(delta_u[delta_u > 0], na.rm = TRUE))) %>%
    pull(min_delta_u)

  if (!is.finite(u_step) || length(u_step) == 0) {
    return(fallback)
  }

  max(fallback, fraction * u_step)
}

# Plot raw pressure data with segmented piecewise curves and CI bands.
create_pressure_plot <- function(raw_data, cutoff_curves, analysis_plot_theme) {
  point_jitter_width <- get_point_jitter_width(raw_data)

  ggplot(raw_data, aes(u, p, color = direction)) +
    geom_ribbon(
      data = cutoff_curves %>% filter(!is.na(p_fit_ci_low)),
      aes(
        x = u,
        ymin = p_fit_ci_low,
        ymax = p_fit_ci_high,
        fill = direction
      ),
      inherit.aes = FALSE,
      alpha = 0.28,
      color = NA
    ) +
    geom_point(
      aes(shape = "Raw samples"),
      alpha = 0.22,
      size = 1.4,
      position = position_jitter(width = point_jitter_width, height = 0, seed = 123)
    ) +
    geom_line(
      data = cutoff_curves %>% filter(!is.na(p_fit)),
      aes(y = p_fit, linetype = direction),
      color = "#1F1F1F",
      linewidth = 0.32
    ) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Branch", labels = c(up = "Upward", down = "Downward")) +
    scale_fill_nejm(guide = "none") +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    scale_linetype_manual(
      name = "Segmented model",
      values = c(up = "solid", down = "longdash"),
      labels = c(up = "Upward model", down = "Downward model")
    ) +
    labs(
      x = "Input voltage [V]",
      y = "Measured pressure [bar]"
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 2, byrow = TRUE, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5)),
      linetype = guide_legend(order = 3, nrow = 2, byrow = TRUE, override.aes = list(color = "#1F1F1F", linewidth = 0.5))
    ) +
    analysis_plot_theme +
    theme(legend.box.just = "center")
}

# Plot raw pressure-voltage scatter before averaging and fitting.
create_raw_pressure_voltage_scatter_plot <- function(raw_data, analysis_plot_theme) {
  point_jitter_width <- get_point_jitter_width(raw_data)

  ggplot(raw_data, aes(u, p, color = direction)) +
    geom_point(
      aes(shape = "Raw samples"),
      alpha = 0.24,
      size = 1.4,
      position = position_jitter(width = point_jitter_width, height = 0, seed = 123)
    ) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Branch", labels = c(up = "Upward", down = "Downward")) +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    labs(
      x = "Input voltage [V]",
      y = "Measured pressure [bar]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5))
    ) +
    analysis_plot_theme
}

# Plot raw pressure-stroke scatter before averaging and fitting.
create_raw_pressure_stroke_scatter_plot <- function(raw_data, analysis_plot_theme) {
  ggplot(raw_data, aes(s, p, color = direction)) +
    geom_point(aes(shape = "Raw samples"), alpha = 0.24, size = 0.85) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Branch", labels = c(up = "Upward", down = "Downward")) +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    labs(
      x = "Stroke [mm]",
      y = "Measured pressure [bar]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5))
    ) +
    analysis_plot_theme
}

# Plot raw voltage-stroke scatter before averaging and fitting.
create_raw_voltage_stroke_scatter_plot <- function(raw_data, analysis_plot_theme) {
  ggplot(raw_data, aes(u, s, color = direction)) +
    geom_point(aes(shape = "Raw samples"), alpha = 0.38, size = 0.38) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Branch", labels = c(up = "Upward", down = "Downward")) +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    labs(
      x = "Input voltage [V]",
      y = "Stroke [mm]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5))
    ) +
    analysis_plot_theme
}

# Plot raw stroke data with mean curves overlaid.
create_stroke_plot <- function(raw_data, cutoff_mean_curves, analysis_plot_theme) {
  mean_curve_plot_df <- cutoff_mean_curves %>%
    arrange(run_id, direction, u)

  ggplot(raw_data, aes(u, s, color = direction)) +
    geom_point(aes(shape = "Raw samples"), alpha = 0.32, size = 0.32) +
    geom_line(
      data = mean_curve_plot_df,
      aes(y = s_mean, group = direction, linetype = direction),
      inherit.aes = TRUE,
      color = "#1F1F1F",
      linewidth = 0.4
    ) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Raw branch", labels = c(up = "Up", down = "Down")) +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    scale_linetype_manual(
      name = "Mean stroke",
      values = c(up = "solid", down = "longdash"),
      labels = c(up = "Upward mean", down = "Downward mean")
    ) +
    labs(
      x = "Input voltage [V]",
      y = "Stroke [mm]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5)),
      linetype = guide_legend(order = 3, override.aes = list(color = "#1F1F1F", linewidth = 0.6))
    ) +
    analysis_plot_theme
}

# Plot raw pressure, segmented model, CI ribbon, and accepted cutoff lines.
create_usable_region_plot <- function(
  raw_data,
  cutoff_curves,
  cutoff_summary,
  analysis_plot_theme
) {
  point_jitter_width <- get_point_jitter_width(raw_data)

  ggplot(raw_data, aes(u, p, color = direction)) +
    geom_ribbon(
      data = cutoff_curves %>% filter(!is.na(p_fit_ci_low)),
      aes(
        x = u,
        ymin = p_fit_ci_low,
        ymax = p_fit_ci_high,
        fill = direction
      ),
      inherit.aes = FALSE,
      alpha = 0.28,
      color = NA
    ) +
    geom_point(
      aes(shape = "Raw samples"),
      alpha = 0.38,
      size = 0.38,
      position = position_jitter(width = point_jitter_width, height = 0, seed = 123)
    ) +
    geom_line(
      data = cutoff_curves %>% filter(!is.na(p_fit)),
      aes(u, p_fit, linetype = direction),
      color = "#1F1F1F",
      linewidth = 0.5
    ) +
    geom_vline(
      data = cutoff_summary %>% filter(fit_valid),
      aes(xintercept = cutoff_low_v, color = direction),
      linetype = "dashed",
      linewidth = 0.2
    ) +
    geom_vline(
      data = cutoff_summary %>% filter(fit_valid),
      aes(xintercept = cutoff_high_v, color = direction),
      linetype = "dotted",
      linewidth = 0.3
    ) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_color_nejm(name = "Branch", labels = c(up = "Upward", down = "Downward")) +
    scale_fill_nejm(guide = "none") +
    scale_shape_manual(name = "Raw data", values = c("Raw samples" = 16)) +
    scale_linetype_manual(
      name = "Segmented model",
      values = c(up = "solid", down = "longdash"),
      labels = c(up = "Up model", down = "Down model")
    ) +
    labs(
      x = "Input voltage [V]",
      y = "Measured pressure [bar]"
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 0.8, size = 1.5)),
      shape = guide_legend(order = 2, override.aes = list(color = "grey20", alpha = 0.8, size = 1.5)),
      linetype = guide_legend(order = 3, override.aes = list(color = "#1F1F1F", linewidth = 0.6))
    ) +
    analysis_plot_theme
}

# Plot the realized signed input-step distribution in each run.
create_step_distribution_plot <- function(raw_data, analysis_plot_theme) {
  raw_data %>%
    group_by(run_id) %>%
    mutate(delta_u = u - lag(u)) %>%
    ungroup() %>%
    filter(!is.na(delta_u)) %>%
    ggplot(aes(delta_u, fill = run_id)) +
    geom_histogram(bins = 40, alpha = 0.85, linewidth = 0.1) +
    facet_wrap(~run_id, scales = "fixed", ncol = 3, labeller = as_labeller(humanize_run_id)) +
    scale_fill_nejm(guide = "none") +
    labs(
      x = expression(Delta * "u [V]"),
      y = "Count"
    ) +
    analysis_plot_theme
}

# Save all generated figures to the analysis output directory.
save_analysis_plots <- function(
  output_dir,
  raw_pressure_voltage_scatter_plot,
  raw_pressure_stroke_scatter_plot,
  raw_voltage_stroke_scatter_plot,
  pressure_plot,
  stroke_plot,
  usable_region_plot,
  step_distribution_plot
) {
  ggsave(
    filename = file.path(output_dir, "pressure_vs_voltage_by_run.png"),
    plot = pressure_plot,
    width = 15,
    height = 5.1724,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "raw_pressure_vs_voltage_scatter.png"),
    plot = raw_pressure_voltage_scatter_plot,
    width = 22,
    height = 8,
    units = "cm",
    dpi = 360
  )

  ggsave(
    filename = file.path(output_dir, "raw_pressure_vs_stroke_scatter.png"),
    plot = raw_pressure_stroke_scatter_plot,
    width = 15,
    height = 5,
    units = "cm",
    dpi = 600
  )

  ggsave(
    filename = file.path(output_dir, "raw_voltage_vs_stroke_scatter.png"),
    plot = raw_voltage_stroke_scatter_plot,
    width = 22,
    height = 8,
    units = "cm",
    dpi = 360
  )

  ggsave(
    filename = file.path(output_dir, "stroke_vs_voltage_by_run.png"),
    plot = stroke_plot,
    width = 22,
    height = 8,
    units = "cm",
    dpi = 360
  )

  ggsave(
    filename = file.path(output_dir, "usable_voltage_cutoffs.png"),
    plot = usable_region_plot,
    width = 22,
    height = 8,
    units = "cm",
    dpi = 360
  )

  ggsave(
    filename = file.path(output_dir, "step_distribution.png"),
    plot = step_distribution_plot,
    width = 20,
    height = 7,
    units = "cm",
    dpi = 320
  )
}

# Print the key summaries in the same order as the analysis pipeline.
print_analysis_summary <- function(
  output_dir,
  run_summary,
  cycle_summary,
  cutoff_summary,
  run_overlap_summary,
  recommended_control_window,
  displacement_validation_flags,
  cutoff_validation,
  settling_summary
) {
  cat("\nRun summary\n")
  print(run_summary)
  cat("\nCycle summary (first 16 rows)\n")
  print(head(cycle_summary, 16))
  cat("\nDirection-level piecewise cutoff summary\n")
  print(cutoff_summary)
  cat("\nRun-level overlap summary\n")
  print(run_overlap_summary)
  cat("\nGlobal recommended control window\n")
  print(filter(recommended_control_window, level == "global"))
  cat("\nDisplacement validation flags\n")
  print(displacement_validation_flags)
  cat("\nCutoff validation across runs\n")
  print(cutoff_validation)
  cat("\nSettling assessment\n")
  print(settling_summary)
  cat(sprintf("\nSaved summary tables and plots to: %s\n", output_dir))
}
