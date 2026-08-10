#!/home/tomislav/bin/micromamba/prefix/envs/pykmd_mkl/bin/python

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def resolve_script_dir() -> Path:
    if "__file__" in globals():
        return Path(__file__).resolve().parent

    cwd = Path.cwd().resolve()
    if (cwd / "scripts" / "lib" / "lifted_lqr_design.py").exists():
        return cwd / "scripts"
    if (cwd / "lib" / "lifted_lqr_design.py").exists():
        return cwd
    return cwd


SCRIPT_DIR = resolve_script_dir()
REPO_DIR = SCRIPT_DIR.parent
LIB_DIR = SCRIPT_DIR / "lib"
EXPECTED_PYTHON = Path("/home/tomislav/bin/micromamba/prefix/envs/pykmd_mkl/bin/python").resolve()
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from lifted_lqr_design import (  # noqa: E402
    ControllerDesignConfig,
    build_lookup_maps,
    compare_to_excitation_schedule,
    extract_steady_state_points,
    identify_candidate_models,
    infer_sampling_time_s,
    load_excitation_experiment,
    load_transient_profiles,
    segment_holds,
)
from pid_baseline_design import (  # noqa: E402
    PIDDesignConfig,
    build_pid_control_limits,
    build_pid_design_summary,
    build_pid_sweep_metrics,
    evaluate_pid_candidates,
    write_pid_controller_limits,
    write_pid_feedback_gains,
)

# Passive low-gain PID refinement:
# - the latest run uses the lower-gain search region while retaining the existing damping,
#   smoothing, and deadband grids
KP_GRID = [0.55, 0.575, 0.60, 0.625, 0.65]
KI_GRID = [0.00070, 0.00090, 0.00110, 0.00130, 0.00150]
KD_GRID = [0.030, 0.032, 0.034, 0.036, 0.038]
SHARED_SEED_TAU_D_S = 0.65
SHARED_SEED_SLOPE_WINDOW_SAMPLES = 15
SHARED_SEED_POSITION_PREFILTER_TAU_S = 0.004
D_TAU_GRID_S = [0.80, 0.85, 0.90, 0.95]
D_SLOPE_WINDOW_GRID_SAMPLES = [17, 19, 21, 23]
D_POSITION_PREFILTER_TAU_GRID_S = [0.004, 0.006, 0.008, 0.010]
BRANCH_DEADBAND_GRID_MM = [0.475, 0.50, 0.525, 0.55, 0.575]
POSITION_DEADBAND_GRID_MM = [0.050, 0.055, 0.060]
TOP_SHARED_PID_CONFIGS = 10
BRANCH_GAIN_FACTORS = [0.85, 1.00, 1.15]

RISE_TIME_BOUNDS = (0.10, 0.90)
SETTLING_BAND_FRAC = 0.02
SETTLING_BAND_FLOOR_MM = 0.15
STEP_METRIC_MIN_STEP_MM = 0.50
SETTLED_WINDOW_FRAC = 0.30


def read_single_global_row(path: Path, level_value: str = "global") -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Missing required analysis file: {path}")
    df = pd.read_csv(path)
    row = df.loc[df["level"] == level_value].head(1)
    if row.empty:
        raise ValueError(f"No level={level_value!r} row found in {path}")
    return row.iloc[0]


def warn_if_nonstandard_python() -> None:
    current_python = Path(sys.executable).resolve()
    if current_python == EXPECTED_PYTHON:
        return
    print(
        "Warning: running with a non-default interpreter.\n"
        f"  preferred: {EXPECTED_PYTHON}\n"
        f"  current:   {current_python}\n"
        "Continuing with the current interpreter because environment layouts may differ across machines."
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Design a branchwise PID baseline from one excitation experiment.")
    parser.add_argument(
        "--input_csv",
        default="experiment/excitation_experiment_100hz_15min_4s_hold.csv",
        help="Single excitation experiment CSV relative to repo root.",
    )
    parser.add_argument(
        "--steady_tail_fraction",
        type=float,
        default=0.05,
        help="Tail fraction used to extract steady-state points.",
    )
    args, unknown_args = parser.parse_known_args()
    if unknown_args:
        print(f"Ignoring unrecognized IDE/launcher arguments: {' '.join(unknown_args)}")
    return args


def build_metric_config() -> ControllerDesignConfig:
    return ControllerDesignConfig(
        q_s_grid=(0.0,),
        q_p_grid=(0.0,),
        q_i_grid=(0.0,),
        r_grid=(0.0,),
        branch_deadband_grid_mm=tuple(BRANCH_DEADBAND_GRID_MM),
        position_deadband_grid_mm=tuple(POSITION_DEADBAND_GRID_MM),
        vel_tau_grid_s=(0.0,),
        vel_slope_window_grid_samples=(0,),
        vel_position_prefilter_tau_grid_s=(0.0,),
        q_v_grid=(0.0,),
        top_base_configs_for_vel=1,
        rise_time_bounds=tuple(RISE_TIME_BOUNDS),
        settling_band_frac=SETTLING_BAND_FRAC,
        settling_band_floor_mm=SETTLING_BAND_FLOOR_MM,
        step_metric_min_step_mm=STEP_METRIC_MIN_STEP_MM,
        settled_window_frac=SETTLED_WINDOW_FRAC,
    )


def build_pid_config() -> PIDDesignConfig:
    return PIDDesignConfig(
        kp_grid=tuple(KP_GRID),
        ki_grid=tuple(KI_GRID),
        kd_grid=tuple(KD_GRID),
        shared_seed_tau_d_s=float(SHARED_SEED_TAU_D_S),
        shared_seed_slope_window_samples=int(SHARED_SEED_SLOPE_WINDOW_SAMPLES),
        shared_seed_position_prefilter_tau_s=float(SHARED_SEED_POSITION_PREFILTER_TAU_S),
        d_tau_grid_s=tuple(D_TAU_GRID_S),
        d_slope_window_grid_samples=tuple(D_SLOPE_WINDOW_GRID_SAMPLES),
        d_position_prefilter_tau_grid_s=tuple(D_POSITION_PREFILTER_TAU_GRID_S),
        branch_deadband_grid_mm=tuple(BRANCH_DEADBAND_GRID_MM),
        position_deadband_grid_mm=tuple(POSITION_DEADBAND_GRID_MM),
        top_shared_pid_configs=TOP_SHARED_PID_CONFIGS,
        branch_gain_factors=tuple(BRANCH_GAIN_FACTORS),
    )


def _format_grid(values: list[float] | tuple[float, ...]) -> str:
    return ", ".join(f"{value:g}" for value in values)


def print_search_plan(args: argparse.Namespace, pid_config: PIDDesignConfig) -> None:
    shared_seed_count = (
        len(pid_config.kp_grid)
        * len(pid_config.ki_grid)
        * len(pid_config.kd_grid)
        * len(pid_config.branch_deadband_grid_mm)
        * len(pid_config.position_deadband_grid_mm)
    )
    smoothing_count = (
        pid_config.top_shared_pid_configs
        * len(pid_config.d_tau_grid_s)
        * len(pid_config.d_slope_window_grid_samples)
        * len(pid_config.d_position_prefilter_tau_grid_s)
    )
    refinement_count = (
        pid_config.top_shared_pid_configs
        * len(pid_config.branch_gain_factors) ** 6
    )
    print("Starting branchwise PID baseline design")
    print(f"Input dataset: {args.input_csv}")
    print(f"Steady-state tail fraction: {args.steady_tail_fraction:g}")
    print(
        "Shared-seed sweep: "
        f"kp=[{_format_grid(pid_config.kp_grid)}], "
        f"ki=[{_format_grid(pid_config.ki_grid)}], "
        f"kd=[{_format_grid(pid_config.kd_grid)}], "
        f"tau_d_s_fixed={pid_config.shared_seed_tau_d_s:g}, "
        f"slope_window_samples_fixed={pid_config.shared_seed_slope_window_samples:g}, "
        f"position_prefilter_tau_s_fixed={pid_config.shared_seed_position_prefilter_tau_s:g}, "
        f"branch_deadband_mm=[{_format_grid(pid_config.branch_deadband_grid_mm)}], "
        f"position_deadband_mm=[{_format_grid(pid_config.position_deadband_grid_mm)}]"
    )
    print(
        "Smoothing refinement: "
        f"top_shared_pid_configs={pid_config.top_shared_pid_configs}, "
        f"tau_d_s=[{_format_grid(pid_config.d_tau_grid_s)}], "
        f"slope_window_samples=[{_format_grid(pid_config.d_slope_window_grid_samples)}], "
        f"position_prefilter_tau_s=[{_format_grid(pid_config.d_position_prefilter_tau_grid_s)}]"
    )
    print(
        "Branchwise gain refinement: "
        f"top_shared_pid_configs={pid_config.top_shared_pid_configs}, "
        f"branch_gain_factors=[{_format_grid(pid_config.branch_gain_factors)}]"
    )
    print(f"Planned shared-seed configurations: {shared_seed_count}")
    print(f"Maximum smoothing-refinement configurations: {smoothing_count}")
    print(f"Maximum branchwise refinement configurations: {refinement_count}")


def main() -> None:
    warn_if_nonstandard_python()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(line_buffering=True)
    args = parse_args()
    metric_config = build_metric_config()
    pid_config = build_pid_config()
    print_search_plan(args, pid_config)

    input_path = (REPO_DIR / args.input_csv).resolve()
    output_dir = REPO_DIR / "analysis_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    control_window_row = read_single_global_row(output_dir / "recommended_control_window.csv")
    excitation_summary_row = pd.read_csv(output_dir / "excitation_signal_summary_4s_hold.csv").iloc[0]

    data_df = segment_holds(load_excitation_experiment(input_path))
    ts_s = infer_sampling_time_s(data_df)
    print(
        "Loaded single excitation dataset: "
        f"{len(data_df)} samples, {data_df['hold_id'].nunique()} holds, inferred Ts={ts_s:.6f} s"
    )
    compare_to_excitation_schedule(data_df, output_dir / "excitation_schedule_4s_hold.csv")

    steady_df = extract_steady_state_points(data_df, ts_s, tail_fraction=args.steady_tail_fraction)
    print(
        "Steady-state extraction complete: "
        f"{int(steady_df['steady_enough'].sum())}/{len(steady_df)} holds accepted"
    )
    lookup_maps = build_lookup_maps(steady_df)
    transient_profiles = load_transient_profiles(
        output_dir,
        use_creep_compensation=True,
        use_dynamic_pressure_reference=False,
    )
    print(
        "Loaded open-loop transient profiles: creep voltage "
        f"through {transient_profiles.creep_by_branch['up'].horizon_s:.2f} s"
    )

    full_hold_ids = (
        data_df.groupby("hold_id", sort=True)["sample_id"]
        .count()
        .loc[lambda s: s > 1]
        .index.to_list()
    )
    base_model = identify_candidate_models(data_df, steady_df, full_hold_ids, metric_config)[
        "lookup_plus_branchwise_minimal_LQI"
    ]
    print(
        "Identified branchwise PID tuning plant: "
        f"full_holds={len(full_hold_ids)} with same-dataset lookup centering"
    )

    provisional_limits = {
        "u_min_v": float(control_window_row["window_low_v"]),
        "u_max_v": float(control_window_row["window_high_v"]),
        "excitation_low_v": float(excitation_summary_row["excitation_low_v"]),
        "excitation_high_v": float(excitation_summary_row["excitation_high_v"]),
        "integral_state_limit": 10.0,
        "branch_selection_rule": "hold_previous_inside_deadband_then_sign_of_position_error",
        "branch_deadband_mm": float(BRANCH_DEADBAND_GRID_MM[0]),
        "position_deadband_mm": float(POSITION_DEADBAND_GRID_MM[0]),
        "anti_windup_mode": "conditional_freeze",
    }

    chosen_result, all_results, sweep_rows = evaluate_pid_candidates(
        model=base_model,
        lookup_maps=lookup_maps,
        ts_s=ts_s,
        control_limits=provisional_limits,
        steady_df=steady_df,
        metric_config=metric_config,
        pid_config=pid_config,
        transient_profiles=transient_profiles,
    )
    print(f"PID evaluation complete: {len(all_results)} aggregate configurations tested")

    summary_df = build_pid_design_summary(
        chosen_result=chosen_result,
        all_results=all_results,
        ts_s=ts_s,
        steady_df=steady_df,
        data_df=data_df,
        tail_fraction=args.steady_tail_fraction,
        input_dataset=str(Path(args.input_csv)),
    )
    summary_df.to_csv(output_dir / "pid_controller_design_summary.csv", index=False)
    build_pid_sweep_metrics(
        sweep_rows,
        all_results,
        chosen_result,
        input_dataset=str(Path(args.input_csv)),
    ).to_csv(
        output_dir / "pid_controller_sweep_metrics.csv",
        index=False,
    )
    write_pid_feedback_gains(output_dir, chosen_result)
    write_pid_controller_limits(
        output_dir,
        build_pid_control_limits(
            control_window_row=control_window_row,
            excitation_summary_row=excitation_summary_row,
            chosen_result=chosen_result,
        ),
    )

    print("\nBranchwise PID baseline design summary")
    print(summary_df)
    print(f"\nChosen config: {chosen_result.config_id}")
    print(f"Chosen selection_score_primary: {chosen_result.selection_score_primary:.6f}")
    print(f"Input dataset path: {input_path}")
    if chosen_result.search_grid_boundary_hit:
        print(f"Chosen config hits search-grid boundary fields: {chosen_result.search_grid_boundary_fields}")
    print(f"Outputs written to: {output_dir}")


if __name__ == "__main__":
    main()
