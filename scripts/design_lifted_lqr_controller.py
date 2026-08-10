#!/home/tomislav/bin/micromamba/prefix/envs/pykmd_mkl/bin/python

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def resolve_script_dir() -> Path:
    if '__file__' in globals():
        return Path(__file__).resolve().parent

    cwd = Path.cwd().resolve()
    if (cwd / 'scripts' / 'lib' / 'lifted_lqr_design.py').exists():
        return cwd / 'scripts'
    if (cwd / 'lib' / 'lifted_lqr_design.py').exists():
        return cwd
    return cwd


SCRIPT_DIR = resolve_script_dir()
REPO_DIR = SCRIPT_DIR.parent
LIB_DIR = SCRIPT_DIR / 'lib'
EXPECTED_PYTHON = Path(
    '/home/tomislav/bin/micromamba/prefix/envs/pykmd_mkl/bin/python'
).resolve()
os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from lifted_lqr_design import (  # noqa: E402
    ControllerDesignConfig,
    build_control_limits,
    build_controller_sweep_metrics,
    build_design_summary,
    build_feedforward_alignment_outputs,
    build_lookup_maps,
    build_lookup_tables,
    compare_to_excitation_schedule,
    evaluate_controllers,
    evaluate_pressure_filter_sensitivity,
    extract_steady_state_points,
    identify_candidate_models,
    infer_sampling_time_s,
    load_excitation_experiment,
    load_transient_profiles,
    pick_final_variant,
    segment_holds,
    write_controller_limits,
    write_feedforward_alignment_outputs,
    write_lookup_tables,
    write_matrix_exports,
    write_nonselected_feedback_gains,
    write_observable_definition,
)

# Focused refinement from the latest 31,500-configuration run:
# - both family winners are at the high q_s and low R boundaries
# - both winners use the largest tested branch deadband and smallest tested position deadband
# - retain the previous winning values inside each shifted grid while probing beyond them
Q_S_GRID = [6.5, 6.75, 7.0, 7.25, 7.5]
Q_P_GRID = [5.0, 5.5, 6.0, 6.5, 7.0]
Q_I_GRID = [0.0005, 0.0007, 0.0009, 0.0011, 0.0013]
R_GRID = [10.0, 10.5, 11.0, 11.5]
BRANCH_DEADBAND_GRID_MM = [0.14, 0.15, 0.16]
POSITION_DEADBAND_GRID_MM = [0.04, 0.05, 0.06]
VEL_POSITION_PREFILTER_TAU_GRID_S = [0.08, 0.10, 0.12]
VEL_TAU_GRID_S = [0.33, 0.36, 0.39]
VEL_SLOPE_WINDOW_GRID_SAMPLES = [27, 29, 31]
Q_V_GRID = [0.02, 0.0225, 0.025, 0.0275, 0.03]

TOP_BASE_CONFIGS_FOR_VEL = 200
RISE_TIME_BOUNDS = (0.10, 0.90)
SETTLING_BAND_FRAC = 0.02
SETTLING_BAND_FLOOR_MM = 0.15
STEP_METRIC_MIN_STEP_MM = 0.50
SETTLED_WINDOW_FRAC = 0.30
LQI_SIMULATION_INTEGRAL_STATE_LIMIT = 20.0


def read_single_global_row(path: Path, level_value: str = 'global') -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f'Missing required analysis file: {path}')
    df = pd.read_csv(path)
    row = df.loc[df['level'] == level_value].head(1)
    if row.empty:
        raise ValueError(f'No level={level_value!r} row found in {path}')
    return row.iloc[0]


def warn_if_nonstandard_python() -> None:
    current_python = Path(sys.executable).resolve()
    if current_python == EXPECTED_PYTHON:
        return
    print(
        'Warning: running with a non-default interpreter.\n'
        f'  preferred: {EXPECTED_PYTHON}\n'
        f'  current:   {current_python}\n'
        'Continuing with the current interpreter because environment layouts may differ across machines.'
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Design minimal branchwise LQI controllers from one excitation experiment.'
    )
    parser.add_argument(
        '--input_csv',
        default='experiment/excitation_experiment_100hz_15min_4s_hold.csv',
        help='Single excitation experiment CSV relative to repo root.',
    )
    parser.add_argument(
        '--steady_tail_fraction',
        type=float,
        default=0.05,
        help='Tail fraction used to extract steady-state points.',
    )
    args, unknown_args = parser.parse_known_args()
    if unknown_args:
        print(f'Ignoring unrecognized IDE/launcher arguments: {" ".join(unknown_args)}')
    return args


def build_design_config() -> ControllerDesignConfig:
    return ControllerDesignConfig(
        q_s_grid=tuple(Q_S_GRID),
        q_p_grid=tuple(Q_P_GRID),
        q_i_grid=tuple(Q_I_GRID),
        r_grid=tuple(R_GRID),
        branch_deadband_grid_mm=tuple(BRANCH_DEADBAND_GRID_MM),
        position_deadband_grid_mm=tuple(POSITION_DEADBAND_GRID_MM),
        vel_tau_grid_s=tuple(VEL_TAU_GRID_S),
        vel_slope_window_grid_samples=tuple(VEL_SLOPE_WINDOW_GRID_SAMPLES),
        vel_position_prefilter_tau_grid_s=tuple(VEL_POSITION_PREFILTER_TAU_GRID_S),
        q_v_grid=tuple(Q_V_GRID),
        top_base_configs_for_vel=TOP_BASE_CONFIGS_FOR_VEL,
        rise_time_bounds=tuple(RISE_TIME_BOUNDS),
        settling_band_frac=SETTLING_BAND_FRAC,
        settling_band_floor_mm=SETTLING_BAND_FLOOR_MM,
        step_metric_min_step_mm=STEP_METRIC_MIN_STEP_MM,
        settled_window_frac=SETTLED_WINDOW_FRAC,
    )


def _format_grid(values: list[float] | tuple[float, ...]) -> str:
    return ', '.join(f'{value:g}' for value in values)


def print_search_plan(
    args: argparse.Namespace, design_config: ControllerDesignConfig
) -> None:
    base_config_count = (
        len(design_config.q_s_grid)
        * len(design_config.q_p_grid)
        * len(design_config.q_i_grid)
        * len(design_config.r_grid)
        * len(design_config.branch_deadband_grid_mm)
        * len(design_config.position_deadband_grid_mm)
    )
    vel_config_count = (
        design_config.top_base_configs_for_vel
        * len(design_config.vel_position_prefilter_tau_grid_s)
        * len(design_config.vel_tau_grid_s)
        * len(design_config.vel_slope_window_grid_samples)
        * len(design_config.q_v_grid)
    )
    print('Starting minimal branchwise LQI controller design')
    print(f'Input dataset: {args.input_csv}')
    print(f'Steady-state tail fraction: {args.steady_tail_fraction:g}')
    print(
        'Base sweep grids: '
        f'q_s=[{_format_grid(design_config.q_s_grid)}], '
        f'q_p=[{_format_grid(design_config.q_p_grid)}], '
        f'q_i=[{_format_grid(design_config.q_i_grid)}], '
        f'r=[{_format_grid(design_config.r_grid)}], '
        f'branch_deadband_mm=[{_format_grid(design_config.branch_deadband_grid_mm)}], '
        f'position_deadband_mm=[{_format_grid(design_config.position_deadband_grid_mm)}]'
    )
    print(
        'Velocity-state staged sweep: '
        f'top_base_configs={design_config.top_base_configs_for_vel}, '
        f'position_prefilter_tau_s=[{_format_grid(design_config.vel_position_prefilter_tau_grid_s)}], '
        f'tau_v_s=[{_format_grid(design_config.vel_tau_grid_s)}], '
        f'slope_window_samples=[{_format_grid(design_config.vel_slope_window_grid_samples)}], '
        f'q_v=[{_format_grid(design_config.q_v_grid)}]'
    )
    print(f'Planned base configurations: {base_config_count}')
    print(f'Maximum staged velocity-state configurations: {vel_config_count}')


def main() -> None:
    warn_if_nonstandard_python()
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(line_buffering=True)
    if hasattr(sys.stderr, 'reconfigure'):
        sys.stderr.reconfigure(line_buffering=True)
    args = parse_args()
    design_config = build_design_config()
    print_search_plan(args, design_config)
    input_path = (REPO_DIR / args.input_csv).resolve()
    output_dir = REPO_DIR / 'analysis_outputs'
    output_dir.mkdir(parents=True, exist_ok=True)

    control_window_row = read_single_global_row(
        output_dir / 'recommended_control_window.csv'
    )
    excitation_summary_row = pd.read_csv(
        output_dir / 'excitation_signal_summary_4s_hold.csv'
    ).iloc[0]

    data_df = segment_holds(load_excitation_experiment(input_path))
    ts_s = infer_sampling_time_s(data_df)
    print(
        'Loaded single excitation dataset: '
        f'{len(data_df)} samples, {data_df["hold_id"].nunique()} holds, inferred Ts={ts_s:.6f} s'
    )

    schedule_check = compare_to_excitation_schedule(
        data_df, output_dir / 'excitation_schedule_4s_hold.csv'
    )

    steady_df = extract_steady_state_points(
        data_df, ts_s, tail_fraction=args.steady_tail_fraction
    )
    steady_df.assign(dataset_role='single_design_dataset').to_csv(
        output_dir / 'lqr_steady_state_points.csv', index=False
    )
    print(
        'Steady-state extraction complete: '
        f'{int(steady_df["steady_enough"].sum())}/{len(steady_df)} holds accepted'
    )

    lookup_maps = build_lookup_maps(steady_df)
    transient_profiles = load_transient_profiles(
        output_dir,
        use_creep_compensation=True,
        use_dynamic_pressure_reference=True,
    )
    print(
        "Loaded open-loop transient profiles: creep voltage + dynamic pressure reference "
        f"through {transient_profiles.creep_by_branch['up'].horizon_s:.2f} s"
    )
    alignment_points_df, alignment_summary_df = build_feedforward_alignment_outputs(
        steady_df, lookup_maps
    )
    alignment_points_df = alignment_points_df.assign(dataset_role='single_design_dataset')
    alignment_summary_df = alignment_summary_df.assign(dataset_role='single_design_dataset')
    lookup_tables = build_lookup_tables(lookup_maps)

    write_feedforward_alignment_outputs(
        alignment_points_df=alignment_points_df,
        alignment_summary_df=alignment_summary_df,
        output_dir=output_dir,
    )
    write_lookup_tables(lookup_tables, output_dir)

    full_hold_ids = (
        data_df.groupby('hold_id', sort=True)['sample_id']
        .count()
        .loc[lambda s: s > 1]
        .index.to_list()
    )

    models = identify_candidate_models(
        data_df, steady_df, full_hold_ids, design_config
    )
    print(
        'Identified candidate branchwise models: '
        f'full_holds={len(full_hold_ids)} (single-dataset design and simulation)'
    )

    provisional_limits = {
        'u_min_v': float(control_window_row['window_low_v']),
        'u_max_v': float(control_window_row['window_high_v']),
        'excitation_low_v': float(excitation_summary_row['excitation_low_v']),
        'excitation_high_v': float(excitation_summary_row['excitation_high_v']),
        'integral_state_limit': LQI_SIMULATION_INTEGRAL_STATE_LIMIT,
        'branch_selection_rule': 'hold_previous_inside_deadband_then_sign_of_position_error',
        'branch_deadband_mm': float(BRANCH_DEADBAND_GRID_MM[0]),
        'position_deadband_mm': float(POSITION_DEADBAND_GRID_MM[0]),
        'anti_windup_mode': 'conditional_freeze',
    }

    controller_results, all_controller_results = evaluate_controllers(
        models=models,
        design_df=data_df,
        design_hold_ids=full_hold_ids,
        evaluation_df=data_df,
        evaluation_hold_ids=full_hold_ids,
        lookup_maps=lookup_maps,
        design_ts_s=ts_s,
        evaluation_ts_s=ts_s,
        control_limits=provisional_limits,
        lookup_steady_df=steady_df,
        config=design_config,
        transient_profiles=transient_profiles,
    )
    print(
        f'Controller evaluation complete: {len(all_controller_results)} aggregate configurations tested'
    )
    chosen_variant = pick_final_variant(controller_results)
    chosen_result = controller_results[chosen_variant]
    pressure_filter_ok, pressure_filter_tau = evaluate_pressure_filter_sensitivity(
        result=chosen_result,
        lookup_maps=lookup_maps,
        ts_s=ts_s,
        control_limits=provisional_limits,
        steady_df=steady_df,
        config=design_config,
        transient_profiles=transient_profiles,
    )
    chosen_result.uses_pressure_filter = pressure_filter_ok
    chosen_result.pressure_filter_tau_s = pressure_filter_tau

    final_limits = build_control_limits(
        control_window_row=control_window_row,
        excitation_summary_row=excitation_summary_row,
        chosen_result=chosen_result,
    )

    summary_df = build_design_summary(
        controller_results=controller_results,
        chosen_variant=chosen_variant,
        ts_s=ts_s,
        steady_df=steady_df,
        data_df=data_df,
        schedule_check=schedule_check,
        tail_fraction=args.steady_tail_fraction,
        input_dataset=str(Path(args.input_csv)),
    )
    summary_df.to_csv(output_dir / 'lqr_controller_design_summary.csv', index=False)
    build_controller_sweep_metrics(
        all_results=all_controller_results,
        family_best_results=controller_results,
        chosen_variant=chosen_variant,
        input_dataset=str(Path(args.input_csv)),
    ).to_csv(output_dir / 'lqr_controller_sweep_metrics.csv', index=False)

    write_observable_definition(output_dir, chosen_result)
    write_matrix_exports(chosen_result, output_dir)
    write_nonselected_feedback_gains(output_dir, controller_results, chosen_variant)
    write_controller_limits(output_dir, final_limits)

    print('\nMinimal branchwise LQI controller design summary')
    print(summary_df)
    print(f'\nChosen variant: {chosen_variant}')
    print(f'Input dataset path: {input_path}')
    chosen_summary_row = summary_df.loc[summary_df['variant'] == chosen_variant].iloc[0]
    print(
        f'Chosen selection_score_primary: {chosen_summary_row["selection_score_primary"]:.6f}'
    )
    non_selected_rows = summary_df.loc[summary_df['variant'] != chosen_variant]
    if not non_selected_rows.empty:
        for _, row in non_selected_rows.iterrows():
            print(
                f'Non-selected {row["variant"]} selection_score_primary: {row["selection_score_primary"]:.6f}'
            )
    tracking_gate_failures = summary_df.loc[
        ~summary_df['selection_tracking_gate_passed'], 'variant'
    ].tolist()
    if tracking_gate_failures:
        print(f'Tracking gate failed for: {", ".join(tracking_gate_failures)}')
    if chosen_result.search_grid_boundary_hit:
        print(
            f'Chosen config hits search-grid boundary fields: {chosen_result.search_grid_boundary_fields}'
        )
    print(f'Outputs written to: {output_dir}')


if __name__ == '__main__':
    main()
