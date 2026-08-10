#!/usr/bin/env python3

"""Evaluate the raw branchwise lookup without residual feedback.

The feedforward-only baseline uses the same identified branchwise plant,
staircase references, saturation, creep layer, and metric implementation as
the active LQI/PID simulations. It is reported for comparison only and is not
included in the feedback-controller selection score.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
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
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
if str(LIB_DIR) not in sys.path:
    sys.path.insert(0, str(LIB_DIR))

from lifted_lqr_design import (  # noqa: E402
    ControllerDesignConfig,
    IdentifiedModel,
    LookupMap,
    _simulate_controller,
    build_reference_scenarios,
    extract_steady_state_points,
    evaluate_prediction_metrics,
    identify_candidate_models,
    infer_sampling_time_s,
    load_excitation_experiment,
    load_transient_profiles,
    segment_holds,
)


DEFAULT_INPUT_CSV = "experiment/excitation_experiment_100hz_15min_4s_hold.csv"
LOOKUP_UP_NAME = "lqr_steady_state_lookup_up.csv"
LOOKUP_DOWN_NAME = "lqr_steady_state_lookup_down.csv"
SUMMARY_NAME = "feedforward_only_controller_summary.csv"
SWEEP_NAME = "feedforward_only_controller_sweep_metrics.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate the branchwise feedforward-only lookup on the identified plant."
    )
    parser.add_argument(
        "--input_csv",
        default=DEFAULT_INPUT_CSV,
        help="Single excitation experiment CSV relative to the repository root.",
    )
    parser.add_argument(
        "--steady_tail_fraction",
        type=float,
        default=0.05,
        help="Tail fraction used to identify the simulation plant.",
    )
    parser.add_argument(
        "--branch_deadband_mm",
        type=float,
        default=None,
        help="Override branch deadband; default is read from lqr_controller_limits.csv.",
    )
    args, unknown_args = parser.parse_known_args()
    if unknown_args:
        print(f"Ignoring unrecognized IDE/launcher arguments: {' '.join(unknown_args)}")
    return args


def _single_value_limits(path: Path) -> pd.Series:
    if not path.exists():
        raise FileNotFoundError(f"Missing controller limits export: {path}")
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"Controller limits export is empty: {path}")
    return df.iloc[0]


def _load_lookup_map(path: Path, expected_branch: str) -> tuple[LookupMap, LookupMap]:
    if not path.exists():
        raise FileNotFoundError(f"Missing active lookup table: {path}")
    df = pd.read_csv(path)
    required = {"branch", "s_ref_mm", "u_ff_v", "p_ff_bar"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Lookup {path} is missing columns: {sorted(missing)}")
    if not np.all(df["branch"].astype(str).eq(expected_branch)):
        raise ValueError(f"Lookup {path} contains an unexpected branch label.")
    values = df[["s_ref_mm", "u_ff_v", "p_ff_bar"]].apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        raise ValueError(f"Lookup {path} contains non-numeric values.")
    if len(values) < 2 or values["s_ref_mm"].duplicated().any():
        raise ValueError(f"Lookup {path} requires at least two unique references.")
    if not np.all(np.diff(values["s_ref_mm"].to_numpy(dtype=float)) > 0.0):
        raise ValueError(f"Lookup {path} must be strictly increasing in s_ref_mm.")
    x = values["s_ref_mm"].to_numpy(dtype=float)
    u = values["u_ff_v"].to_numpy(dtype=float)
    p = values["p_ff_bar"].to_numpy(dtype=float)
    return (
        LookupMap(expected_branch, "s_ref_mm", "u_ff_v", x, u),
        LookupMap(expected_branch, "s_ref_mm", "p_ff_bar", x, p),
    )


def load_active_lookup_maps(output_dir: Path) -> dict[str, LookupMap]:
    up_u, up_p = _load_lookup_map(output_dir / LOOKUP_UP_NAME, "up")
    down_u, down_p = _load_lookup_map(output_dir / LOOKUP_DOWN_NAME, "down")
    return {
        "up_inverse_u": up_u,
        "up_inverse_p": up_p,
        "down_inverse_u": down_u,
        "down_inverse_p": down_p,
    }


def build_identification_config() -> ControllerDesignConfig:
    singleton = (0.0,)
    return ControllerDesignConfig(
        q_s_grid=singleton,
        q_p_grid=singleton,
        q_i_grid=singleton,
        r_grid=singleton,
        branch_deadband_grid_mm=singleton,
        position_deadband_grid_mm=singleton,
        vel_tau_grid_s=singleton,
        vel_slope_window_grid_samples=(2,),
        vel_position_prefilter_tau_grid_s=singleton,
        q_v_grid=singleton,
        top_base_configs_for_vel=1,
        rise_time_bounds=(0.10, 0.90),
        settling_band_frac=0.02,
        settling_band_floor_mm=0.15,
        step_metric_min_step_mm=0.50,
        settled_window_frac=0.30,
    )


def _base_zero_feedback_gains(model: IdentifiedModel) -> dict[str, np.ndarray]:
    return {
        branch: np.zeros((model.A.shape[0] + 1,), dtype=float)
        for branch in ("up", "down")
    }


def _build_control_limits(
    output_dir: Path,
    limits_row: pd.Series,
    branch_deadband_override: float | None,
) -> dict[str, float | str | bool]:
    branch_deadband = (
        float(branch_deadband_override)
        if branch_deadband_override is not None
        else float(limits_row["branch_deadband_mm"])
    )
    return {
        "u_min_v": float(limits_row["u_min_v"]),
        "u_max_v": float(limits_row["u_max_v"]),
        "integral_state_limit": 0.0,
        "branch_deadband_mm": branch_deadband,
        "position_deadband_mm": float(limits_row["position_deadband_mm"]),
    }


def _metric_columns() -> list[str]:
    return [
        "tracking_rms_mm",
        "tracking_overshoot_mm",
        "control_effort_mean_abs_v",
        "settled_error_rms_mm",
        "settled_error_p2p_mm",
        "settled_command_p2p_v",
        "settled_velocity_rms_mm_s",
        "rise_time_s_10_90",
        "settling_time_s_2pct",
        "peak_time_s",
        "overshoot_mm",
        "overshoot_pct",
        "steady_state_bias_mm",
        "IAE_mm_s",
        "ISE_mm2_s",
        "ITAE_mm_s2",
        "control_residual_rms_v",
        "command_total_variation_v",
        "saturation_fraction",
    ]


def _simulation_valid(metrics: dict[str, float]) -> bool:
    """Allow undefined step metrics while requiring core metrics to be finite."""
    optional_step_metrics = {
        "rise_time_s_10_90",
        "settling_time_s_2pct",
        "peak_time_s",
        "overshoot_mm",
        "overshoot_pct",
    }
    return all(
        np.isfinite(float(value))
        for name, value in metrics.items()
        if name not in optional_step_metrics
    )


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    args = parse_args()
    input_path = (REPO_DIR / args.input_csv).resolve()
    output_dir = REPO_DIR / "analysis_outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Starting feedforward-only lookup evaluation")
    print(f"Input dataset: {args.input_csv}")
    print(f"Lookup tables: {LOOKUP_UP_NAME}, {LOOKUP_DOWN_NAME}")

    data_df = segment_holds(load_excitation_experiment(input_path))
    ts_s = infer_sampling_time_s(data_df)
    steady_df = extract_steady_state_points(
        data_df, ts_s, tail_fraction=args.steady_tail_fraction
    )
    full_hold_ids = (
        data_df.groupby("hold_id", sort=True)["sample_id"]
        .count()
        .loc[lambda counts: counts > 1]
        .index.to_list()
    )
    lookup_maps = load_active_lookup_maps(output_dir)
    transient_profiles = load_transient_profiles(
        output_dir,
        use_creep_compensation=True,
        use_dynamic_pressure_reference=False,
    )
    limits_row = _single_value_limits(output_dir / "lqr_controller_limits.csv")
    control_limits = _build_control_limits(
        output_dir, limits_row, args.branch_deadband_mm
    )

    models = identify_candidate_models(
        data_df, steady_df, full_hold_ids, build_identification_config()
    )
    model = models["lookup_plus_branchwise_minimal_LQI"]
    prediction_metrics = evaluate_prediction_metrics(
        model,
        data_df,
        full_hold_ids,
        steady_df,
        ts_s,
    )
    branch_gains = _base_zero_feedback_gains(model)
    lookup_tables = {
        branch: pd.DataFrame(
            {
                "branch": branch,
                "s_ref_mm": lookup_maps[f"{branch}_inverse_u"].grid_x,
                "u_ff_v": lookup_maps[f"{branch}_inverse_u"].grid_y,
                "p_ff_bar": lookup_maps[f"{branch}_inverse_p"].grid_y,
            }
        )
        for branch in ("up", "down")
    }
    hold_steps = max(60, int(round(1.1 / ts_s)))
    scenarios = build_reference_scenarios(lookup_tables, hold_steps)
    config = build_identification_config()

    scenario_metrics: dict[str, dict[str, float]] = {}
    sweep_rows: list[dict[str, object]] = []
    for scenario_name, ref_values in scenarios.items():
        print(f"Evaluating scenario: {scenario_name}")
        metrics = _simulate_controller(
            model=model,
            lookup_maps=lookup_maps,
            ref_mm=ref_values,
            ts_s=ts_s,
            control_limits=control_limits,
            steady_df=steady_df,
            config=config,
            branch_K_aug=branch_gains,
            transient_profiles=transient_profiles,
            use_creep_compensation=True,
            use_dynamic_pressure_reference=False,
        )
        scenario_metrics[scenario_name] = metrics
        sweep_rows.append(
            {
                "variant": "lookup_only_feedforward",
                "scenario": scenario_name,
                "input_dataset": str(Path(args.input_csv)),
                "data_split_policy": "single_dataset_no_hold_split",
                "selection_dataset_role": "single_dataset_simulation",
                "controller_family": "feedforward_only_baseline",
                "uses_velocity_state": False,
                "uses_creep_compensation": True,
                "uses_dynamic_pressure_reference": False,
                "transient_profile_horizon_s": 3.99,
                "transient_lookup_source": "open_loop_steady_state_lookup",
                "n_runtime_observables": 0,
                "branch_deadband_mm": control_limits["branch_deadband_mm"],
                "position_deadband_mm": control_limits["position_deadband_mm"],
                "validation_ok": _simulation_valid(metrics),
                "selection_score_primary": np.nan,
                "selection_eligible": False,
                "selection_note": "not_ranked_no_feedback_gain_or_integrator",
                **prediction_metrics,
                **metrics,
            }
        )

    from lifted_lqr_design import aggregate_scenario_metrics

    aggregate = aggregate_scenario_metrics(scenario_metrics)
    sweep_rows.append(
        {
            "variant": "lookup_only_feedforward",
            "scenario": "all",
            "input_dataset": str(Path(args.input_csv)),
            "data_split_policy": "single_dataset_no_hold_split",
            "selection_dataset_role": "single_dataset_simulation",
            "controller_family": "feedforward_only_baseline",
            "uses_velocity_state": False,
            "uses_creep_compensation": True,
            "uses_dynamic_pressure_reference": False,
            "transient_profile_horizon_s": 3.99,
            "transient_lookup_source": "open_loop_steady_state_lookup",
            "n_runtime_observables": 0,
            "branch_deadband_mm": control_limits["branch_deadband_mm"],
            "position_deadband_mm": control_limits["position_deadband_mm"],
            "validation_ok": _simulation_valid(aggregate),
            "selection_score_primary": np.nan,
            "selection_eligible": False,
            "selection_note": "not_ranked_no_feedback_gain_or_integrator",
            **prediction_metrics,
            **aggregate,
        }
    )

    sweep_df = pd.DataFrame(sweep_rows)
    sweep_df.to_csv(output_dir / SWEEP_NAME, index=False)

    summary = {
        "variant": "lookup_only_feedforward",
        "chosen_variant": False,
        "global_selected": False,
        "evaluated_variant": True,
        "controller_family": "feedforward_only_baseline",
        "input_dataset": str(Path(args.input_csv)),
        "data_split_policy": "single_dataset_no_hold_split",
        "prediction_dataset_role": "single_dataset_in_sample",
        "selection_dataset_role": "single_dataset_simulation",
        "sampling_time_s": ts_s,
        "n_samples": len(data_df),
        "n_holds_total": int(data_df["hold_id"].nunique()),
        "n_holds_full": len(full_hold_ids),
        "n_steady_points_accepted": int(steady_df["steady_enough"].sum()),
        "uses_creep_compensation": True,
        "uses_dynamic_pressure_reference": False,
        "transient_profile_horizon_s": 3.99,
        "transient_lookup_source": "open_loop_steady_state_lookup",
        "n_runtime_observables": 0,
        "branch_deadband_mm": control_limits["branch_deadband_mm"],
        "position_deadband_mm": control_limits["position_deadband_mm"],
        "integral_state_limit": 0.0,
        "stable_model": model.stable_model,
        "spectral_radius": model.spectral_radius,
        **prediction_metrics,
        "gain_l1_norm": 0.0,
        "max_abs_gain_s": 0.0,
        "max_abs_gain_p": 0.0,
        "min_abs_integral_gain": np.nan,
        "selection_score_primary": np.nan,
        "selection_eligible": False,
        "selection_note": "comparison_baseline_not_ranked_against_feedback_controllers",
        **aggregate,
    }
    pd.DataFrame([summary]).to_csv(output_dir / SUMMARY_NAME, index=False)
    print(f"Aggregate evaluation complete: tracking_rms_mm={aggregate['tracking_rms_mm']:.6f}")
    print(f"Wrote: {output_dir / SUMMARY_NAME}")
    print(f"Wrote: {output_dir / SWEEP_NAME}")


if __name__ == "__main__":
    main()
