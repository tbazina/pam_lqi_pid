# pam_lqi_pid

Repository for data-driven position control of a FESTO DMSP-10-250 RM-CM
McKibben pneumatic artificial muscle (PAM). The project combines open-loop
actuator characterization, branchwise feedforward lookup construction,
control-oriented model identification, simulated controller tuning, and
LabVIEW/PAM validation.

## Active workflow

The current controller-design workflow uses the single active experiment
`experiment/excitation_experiment_100hz_15min_4s_hold.csv`:

1. Static threshold, hysteresis, sensor calibration, sampling, settling, and
   repeatability analyses establish the usable operating region and timing.
2. The conditioned 4 s voltage excitation data are used to build raw branchwise
   displacement-to-voltage and displacement-to-pressure lookups.
3. Branchwise discrete plant models are identified from the same open-loop data.
4. Open-loop creep-voltage and dynamic-pressure profiles are identified. The
   transient layers are applied during simulated evaluation before controller
   ranking.
5. Python drivers tune and compare base LQI, velocity-state LQI, and branchwise
   PID controllers on synthetic staircase references. A lookup-plus-creep
   feedforward-only baseline is evaluated separately.
6. Controller gains, limits, lookups, observables, matrices, and metrics
   are exported to `analysis_outputs/`.

Acquisition conditioning is performed in LabVIEW before offline processing:
700 displacement oversamples and 90 pressure oversamples are averaged, rounded
to `0.02 mm` and `0.02 bar`, respectively, and passed through a three-sample median
filter. 

The design is single-dataset and in-sample. Earlier v2/v3/v4 train-test runs,
controller-specific lookup post-correction, and real final-validation logs are
retained as archival.

## Repository layout

```text
experiment/          Raw excitation, calibration, post-correction, and validation CSV files.
analysis_outputs/    Generated lookups, summaries, controller exports, runtimes, and figures.
scripts/             Python controller-design and R analysis/generation entry points.
scripts/lib/         Shared Python and R helpers.
```

## Important scripts

Run commands from the repository root.

Python controller design:

- `scripts/design_lifted_lqr_controller.py` searches base and velocity-state
  LQI configurations.
- `scripts/design_branchwise_pid_controller.py` searches the branchwise PID
  baseline.
- `scripts/evaluate_feedforward_only_controller.py` evaluates the separate
  feedforward-only comparison baseline.
- `scripts/lib/lifted_lqr_design.py` contains shared lookup, model,
  simulation, metrics, selection, and export logic.
- `scripts/lib/pid_baseline_design.py` contains PID search and export logic.

R analysis and generation:

- `scripts/analyze_static_threshold_hysteresis.R`
- `scripts/analyze_pressure_sensor_calibration.R`
- `scripts/analyze_effective_settling_time.R`
- `scripts/analyze_max_sampling_time.R`
- `scripts/generate_excitation_signal.R`
- `scripts/analyze_local_static_map_repeatability.R`
- `scripts/analyze_open_loop_creep_compensation.R`
- `scripts/analyze_open_loop_dynamic_pressure_reference.R`
- `scripts/generate_reference_validation_signal.R`
- `scripts/analyze_final_controller_validation.R`
- `scripts/plot_open_loop_transient_phenomena_example.R`
- `scripts/plot_steady_state_feedforward_lookup.R`

Shared R plotting and preprocessing helpers are in `scripts/lib/`.

The R environment uses `readr`, `dplyr`, `stringr`, `purrr`, `tidyr`,
`ggplot2`, `ggsci`, `segmented`, and `tibble`. Python controller scripts use
`numpy`, `pandas`, and `scipy`.
