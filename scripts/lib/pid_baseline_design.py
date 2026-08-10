from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from lifted_lqr_design import (
    ControllerDesignConfig,
    IdentifiedModel,
    LookupMap,
    TransientProfiles,
    SELECTION_GAIN_PENALTY_EXPONENT,
    SELECTION_MIN_ABS_INTEGRAL_GAIN,
    SELECTION_PRIMARY_METRIC_EXPORT_NAMES,
    SELECTION_SCORE_FORMULA,
    SELECTION_TIE_BREAK_EXPORT_NAMES,
    SELECTION_TRACKING_GATE_FACTOR,
    aggregate_scenario_metrics,
    alpha_from_tau,
    build_global_means,
    build_lookup_tables,
    build_reference_scenarios,
    compute_first_order_lowpass_signal,
    compute_latest_linear_fit_slope,
    compute_segment_metrics,
    rank_results_for_selection,
    segment_bounds,
    update_first_order_lowpass,
)


@dataclass(frozen=True)
class PIDDesignConfig:
    kp_grid: Tuple[float, ...]
    ki_grid: Tuple[float, ...]
    kd_grid: Tuple[float, ...]
    shared_seed_tau_d_s: float
    shared_seed_slope_window_samples: int
    shared_seed_position_prefilter_tau_s: float
    d_tau_grid_s: Tuple[float, ...]
    d_slope_window_grid_samples: Tuple[int, ...]
    d_position_prefilter_tau_grid_s: Tuple[float, ...]
    branch_deadband_grid_mm: Tuple[float, ...]
    position_deadband_grid_mm: Tuple[float, ...]
    top_shared_pid_configs: int
    branch_gain_factors: Tuple[float, ...]


@dataclass
class PIDControllerResult:
    config_id: str
    stage: str
    variant: str
    controller_family: str
    kp_up: float
    ki_up: float
    kd_up: float
    kp_down: float
    ki_down: float
    kd_down: float
    tau_d_s: float
    slope_window_samples: int
    position_prefilter_tau_s: float
    derivative_filter_alpha: float
    derivative_position_prefilter_alpha: float
    branch_deadband_mm: float
    position_deadband_mm: float
    kp_shared: float
    ki_shared: float
    kd_shared: float
    kp_up_factor: float | None
    ki_up_factor: float | None
    kd_up_factor: float | None
    kp_down_factor: float | None
    ki_down_factor: float | None
    kd_down_factor: float | None
    tracking_rms_mm: float
    tracking_overshoot_mm: float
    control_effort_mean_abs_v: float
    settled_error_rms_mm: float
    settled_error_p2p_mm: float
    settled_command_p2p_v: float
    settled_velocity_rms_mm_s: float
    rise_time_s_10_90: float
    settling_time_s_2pct: float
    peak_time_s: float
    overshoot_mm: float
    overshoot_pct: float
    steady_state_bias_mm: float
    abs_steady_state_bias_mm: float
    iae_mm_s: float
    ise_mm2_s: float
    itae_mm_s2: float
    control_residual_rms_v: float
    command_total_variation_v: float
    saturation_fraction: float
    validation_ok: bool
    invalid_reason: str
    gain_l1_norm: float
    max_abs_gain_s: float
    min_abs_integral_gain: float
    search_grid_boundary_hit: bool
    search_grid_boundary_fields: str
    selection_score_primary: float = float("inf")
    selection_tracking_gate_passed: bool = False
    selection_eligibility_passed: bool = False
    uses_creep_compensation: bool = True
    uses_dynamic_pressure_reference: bool = False
    transient_profile_horizon_s: float | None = 3.99
    transient_lookup_source: str = "open_loop_steady_state_lookup"


_PID_WORKER_MODEL: IdentifiedModel | None = None
_PID_WORKER_REFERENCE_SCENARIOS: Dict[str, np.ndarray] | None = None
_PID_WORKER_BASE_LIMITS: Dict[str, float | str | bool] | None = None
_PID_WORKER_METRIC_CONFIG: ControllerDesignConfig | None = None
_PID_WORKER_TS_S: float | None = None
_PID_WORKER_X_BAR_2: np.ndarray | None = None
_PID_WORKER_U_BAR: float | None = None
_PID_WORKER_SIM_CONTEXT: Dict[str, np.ndarray | float] | None = None
_PID_WORKER_PID_CONFIG: PIDDesignConfig | None = None
_PID_WORKER_TRANSIENT_PROFILES: TransientProfiles | None = None


def _pid_gain_l1_norm(
    kp_up: float,
    ki_up: float,
    kd_up: float,
    kp_down: float,
    ki_down: float,
    kd_down: float,
) -> float:
    return float(
        abs(kp_up) + abs(ki_up) + abs(kd_up) + abs(kp_down) + abs(ki_down) + abs(kd_down)
    )


def _pid_max_abs_gain_s(
    kp_up: float,
    kp_down: float,
) -> float:
    return float(max(abs(kp_up), abs(kp_down)))


def _pid_min_abs_integral_gain(
    ki_up: float,
    ki_down: float,
) -> float:
    return float(min(abs(ki_up), abs(ki_down)))


def _init_pid_worker(
    model: IdentifiedModel,
    reference_scenarios: Dict[str, np.ndarray],
    base_limits: Dict[str, float | str | bool],
    metric_config: ControllerDesignConfig,
    ts_s: float,
    x_bar_2: np.ndarray,
    u_bar: float,
    sim_context: Dict[str, np.ndarray | float],
    pid_config: PIDDesignConfig,
    transient_profiles: TransientProfiles,
) -> None:
    global _PID_WORKER_MODEL
    global _PID_WORKER_REFERENCE_SCENARIOS
    global _PID_WORKER_BASE_LIMITS
    global _PID_WORKER_METRIC_CONFIG
    global _PID_WORKER_TS_S
    global _PID_WORKER_X_BAR_2
    global _PID_WORKER_U_BAR
    global _PID_WORKER_SIM_CONTEXT
    global _PID_WORKER_PID_CONFIG
    global _PID_WORKER_TRANSIENT_PROFILES

    _PID_WORKER_MODEL = model
    _PID_WORKER_REFERENCE_SCENARIOS = reference_scenarios
    _PID_WORKER_BASE_LIMITS = base_limits
    _PID_WORKER_METRIC_CONFIG = metric_config
    _PID_WORKER_TS_S = ts_s
    _PID_WORKER_X_BAR_2 = x_bar_2
    _PID_WORKER_U_BAR = u_bar
    _PID_WORKER_SIM_CONTEXT = sim_context
    _PID_WORKER_PID_CONFIG = pid_config
    _PID_WORKER_TRANSIENT_PROFILES = transient_profiles


def _pid_runtime_limits(
    control_limits: Dict[str, float | str | bool],
    ki_up: float,
    ki_down: float,
) -> Dict[str, float | str | bool]:
    limits = dict(control_limits)
    ki_max = max(abs(float(ki_up)), abs(float(ki_down)))
    if ki_max < 1e-6:
        integral_limit = 10.0
    else:
        excitation_span = float(limits["excitation_high_v"]) - float(limits["excitation_low_v"])
        integral_limit = 0.40 * excitation_span / ki_max
    limits["integral_state_limit"] = float(integral_limit)
    return limits


def _practical_displacement_bounds(lookup_maps: Dict[str, LookupMap]) -> Tuple[float, float]:
    s_values = np.concatenate(
        [
            lookup_maps["up_inverse_u"].grid_x,
            lookup_maps["down_inverse_u"].grid_x,
        ]
    )
    s_min = float(np.min(s_values))
    s_max = float(np.max(s_values))
    span = float(s_max - s_min)
    return float(s_min - 3.0 * span), float(s_max + 3.0 * span)


def _build_pid_simulation_context(
    model: IdentifiedModel,
    lookup_maps: Dict[str, LookupMap],
) -> Dict[str, np.ndarray | float]:
    s_lower_bound, s_upper_bound = _practical_displacement_bounds(lookup_maps)
    return {
        "u_lookup_up_x": lookup_maps["up_inverse_u"].grid_x.astype(float, copy=False),
        "u_lookup_up_y": lookup_maps["up_inverse_u"].grid_y.astype(float, copy=False),
        "u_lookup_down_x": lookup_maps["down_inverse_u"].grid_x.astype(float, copy=False),
        "u_lookup_down_y": lookup_maps["down_inverse_u"].grid_y.astype(float, copy=False),
        "A_up": model.branch_models["up"][0],
        "B_up": model.branch_models["up"][1][:, 0],
        "A_down": model.branch_models["down"][0],
        "B_down": model.branch_models["down"][1][:, 0],
        "s_lower_bound": float(s_lower_bound),
        "s_upper_bound": float(s_upper_bound),
    }


def _simulate_pid_controller(
    model: IdentifiedModel,
    ref_mm: np.ndarray,
    ts_s: float,
    control_limits: Dict[str, float | str | bool],
    metric_config: ControllerDesignConfig,
    branch_gains: Dict[str, Tuple[float, float, float]],
    tau_d_s: float,
    slope_window_samples: int,
    position_prefilter_tau_s: float,
    x_bar_2: np.ndarray,
    u_bar: float,
    sim_context: Dict[str, np.ndarray | float],
    transient_profiles: TransientProfiles,
) -> Tuple[Dict[str, float], bool, str]:
    x = np.zeros(model.A.shape[0], dtype=float)
    e_int = 0.0
    s_meas_hist: List[float] = [float(x_bar_2[0])]
    s_prefilt_hist: List[float] = [float(x_bar_2[0])]
    v_hat = 0.0
    s_prefilt_state = float(x_bar_2[0])
    active_branch = "up"
    branch_deadband = float(control_limits["branch_deadband_mm"])
    position_deadband = float(control_limits["position_deadband_mm"])
    alpha_d = alpha_from_tau(ts_s, tau_d_s)
    alpha_s = alpha_from_tau(ts_s, position_prefilter_tau_s)
    u_lookup_up_x = sim_context["u_lookup_up_x"]
    u_lookup_up_y = sim_context["u_lookup_up_y"]
    u_lookup_down_x = sim_context["u_lookup_down_x"]
    u_lookup_down_y = sim_context["u_lookup_down_y"]
    A_up = sim_context["A_up"]
    B_up = sim_context["B_up"]
    A_down = sim_context["A_down"]
    B_down = sim_context["B_down"]
    s_lower_bound = float(sim_context["s_lower_bound"])
    s_upper_bound = float(sim_context["s_upper_bound"])

    s_hist_mm: List[float] = []
    err_hist_mm: List[float] = []
    u_cmd_hist_v: List[float] = []
    u_ff_hist_v: List[float] = []
    v_slope_raw_hist_mm_s: List[float] = []
    sat_hist: List[bool] = []
    initial_output_mm = float(x_bar_2[0])

    invalid_reason = ""
    valid = True
    previous_ref_mm: float | None = None
    step_delta_s_mm = 0.0
    step_branch: str | None = None
    time_since_step_s = 0.0

    for s_ref_value in ref_mm:
        s_ref = float(s_ref_value)
        s_meas_now = float(x[0] + x_bar_2[0])
        s_tilde = s_meas_now - s_ref
        s_tilde_fb = 0.0 if abs(s_tilde) < position_deadband else s_tilde

        if s_ref > s_meas_now + branch_deadband:
            active_branch = "up"
        elif s_ref < s_meas_now - branch_deadband:
            active_branch = "down"

        if previous_ref_mm is None:
            time_since_step_s = 0.0
            step_delta_s_mm = 0.0
            step_branch = None
        elif abs(s_ref - previous_ref_mm) > 1e-12:
            step_delta_s_mm = s_ref - previous_ref_mm
            step_branch = "up" if step_delta_s_mm > 0.0 else "down"
            time_since_step_s = 0.0

        if active_branch == "up":
            u_ff_static = float(np.interp(s_ref, u_lookup_up_x, u_lookup_up_y))
            A_eff = A_up
            B_eff = B_up
        else:
            u_ff_static = float(np.interp(s_ref, u_lookup_down_x, u_lookup_down_y))
            A_eff = A_down
            B_eff = B_down
        u_ff = u_ff_static
        direction_matches = step_branch == active_branch and abs(step_delta_s_mm) > 1e-12
        creep_profile = transient_profiles.creep_by_branch.get(active_branch)
        if creep_profile is None:
            raise ValueError("Creep compensation is enabled but no branch profile is loaded.")
        if direction_matches and time_since_step_s <= creep_profile.horizon_s + 1e-12:
            creep_fraction = creep_profile.evaluate(time_since_step_s)
            if active_branch == "up":
                lookup_x, lookup_y = u_lookup_up_x, u_lookup_up_y
            else:
                lookup_x, lookup_y = u_lookup_down_x, u_lookup_down_y
            s_virtual = np.clip(s_ref + step_delta_s_mm * creep_fraction, lookup_x[0], lookup_x[-1])
            u_ff = float(np.interp(s_virtual, lookup_x, lookup_y))
        s_prefilt_state = update_first_order_lowpass(s_prefilt_state, s_meas_now, alpha_s)
        v_slope_raw = compute_latest_linear_fit_slope(
            s_prefilt_hist + [s_prefilt_state],
            ts_s,
            slope_window_samples,
        )
        v_hat = alpha_d * v_hat + (1.0 - alpha_d) * v_slope_raw

        kp_branch, ki_branch, kd_branch = branch_gains[active_branch]
        u_cmd_raw = float(u_ff - kp_branch * s_tilde_fb - ki_branch * e_int - kd_branch * v_hat)
        u_cmd = float(np.clip(u_cmd_raw, float(control_limits["u_min_v"]), float(control_limits["u_max_v"])))
        saturated = bool(
            u_cmd >= float(control_limits["u_max_v"]) - 1e-12
            or u_cmd <= float(control_limits["u_min_v"]) + 1e-12
        )

        u_tilde = u_cmd - u_bar
        x_next = A_eff @ x + B_eff * u_tilde
        if not (np.all(np.isfinite(x_next)) and np.isfinite(u_cmd) and np.isfinite(v_hat)):
            valid = False
            invalid_reason = "nonfinite_state"
            break

        s_meas_next = float(x_next[0] + x_bar_2[0])
        if s_meas_next < s_lower_bound or s_meas_next > s_upper_bound:
            valid = False
            invalid_reason = "displacement_out_of_range"
            break

        error_now = s_meas_now - s_ref
        if u_cmd >= float(control_limits["u_max_v"]) - 1e-12 and error_now < 0:
            pass
        elif u_cmd <= float(control_limits["u_min_v"]) + 1e-12 and error_now > 0:
            pass
        else:
            e_int = float(e_int + (s_meas_next - s_ref) * ts_s)
        e_int = float(
            np.clip(
                e_int,
                -float(control_limits["integral_state_limit"]),
                float(control_limits["integral_state_limit"]),
            )
        )

        s_hist_mm.append(s_meas_next)
        err_hist_mm.append(s_meas_next - s_ref)
        u_cmd_hist_v.append(u_cmd)
        u_ff_hist_v.append(u_ff)
        v_slope_raw_hist_mm_s.append(v_slope_raw)
        sat_hist.append(saturated)

        s_meas_hist.append(s_meas_now)
        s_prefilt_hist.append(s_prefilt_state)
        x = x_next
        previous_ref_mm = s_ref
        time_since_step_s += ts_s

    if not valid or not s_hist_mm:
        nan_metrics = {
            "tracking_rms_mm": float("nan"),
            "tracking_overshoot_mm": float("nan"),
            "control_effort_mean_abs_v": float("nan"),
            "settled_error_rms_mm": float("nan"),
            "settled_error_p2p_mm": float("nan"),
            "settled_command_p2p_v": float("nan"),
            "settled_velocity_rms_mm_s": float("nan"),
            "rise_time_s_10_90": float("nan"),
            "settling_time_s_2pct": float("nan"),
            "peak_time_s": float("nan"),
            "overshoot_mm": float("nan"),
            "overshoot_pct": float("nan"),
            "steady_state_bias_mm": float("nan"),
            "IAE_mm_s": float("nan"),
            "ISE_mm2_s": float("nan"),
            "ITAE_mm_s2": float("nan"),
            "control_residual_rms_v": float("nan"),
            "command_total_variation_v": float("nan"),
            "saturation_fraction": float("nan"),
        }
        return nan_metrics, False, invalid_reason or "simulation_invalid"

    s_hist_arr = np.asarray(s_hist_mm, dtype=float)
    err_hist_arr = np.asarray(err_hist_mm, dtype=float)
    u_cmd_hist_arr = np.asarray(u_cmd_hist_v, dtype=float)
    u_ff_hist_arr = np.asarray(u_ff_hist_v, dtype=float)
    v_raw_hist_arr = np.asarray(v_slope_raw_hist_mm_s, dtype=float)
    sat_hist_arr = np.asarray(sat_hist, dtype=bool)

    segment_metric_rows: List[Dict[str, float]] = []
    for start, end in segment_bounds(ref_mm):
        y0_mm = initial_output_mm if start == 0 else float(s_hist_arr[start - 1])
        segment_metric_rows.append(
            compute_segment_metrics(
                y0_mm=y0_mm,
                ref_level_mm=float(ref_mm[start]),
                y_seg_mm=s_hist_arr[start:end],
                err_seg_mm=err_hist_arr[start:end],
                u_cmd_seg_v=u_cmd_hist_arr[start:end],
                u_ff_seg_v=u_ff_hist_arr[start:end],
                v_seg_mm_s=v_raw_hist_arr[start:end],
                sat_seg=sat_hist_arr[start:end],
                ts_s=ts_s,
                config=metric_config,
            )
        )

    scenario_itae = float(np.sum([row["ITAE_mm_s2"] for row in segment_metric_rows]))
    scenario_metrics = {
        "tracking_rms_mm": float(np.sqrt(np.mean(np.square(err_hist_arr)))),
        "tracking_overshoot_mm": float(np.max(np.abs(err_hist_arr))),
        "control_effort_mean_abs_v": float(np.mean(np.abs(u_cmd_hist_arr - u_ff_hist_arr))),
        "settled_error_rms_mm": float(np.nanmean([row["settled_error_rms_mm"] for row in segment_metric_rows])),
        "settled_error_p2p_mm": float(np.nanmean([row["settled_error_p2p_mm"] for row in segment_metric_rows])),
        "settled_command_p2p_v": float(np.nanmean([row["settled_command_p2p_v"] for row in segment_metric_rows])),
        "settled_velocity_rms_mm_s": float(
            np.nanmean([row["settled_velocity_rms_mm_s"] for row in segment_metric_rows])
        ),
        "rise_time_s_10_90": float(np.nanmean([row["rise_time_s_10_90"] for row in segment_metric_rows])),
        "settling_time_s_2pct": float(np.nanmean([row["settling_time_s_2pct"] for row in segment_metric_rows])),
        "peak_time_s": float(np.nanmean([row["peak_time_s"] for row in segment_metric_rows])),
        "overshoot_mm": float(np.nanmean([row["overshoot_mm"] for row in segment_metric_rows])),
        "overshoot_pct": float(np.nanmean([row["overshoot_pct"] for row in segment_metric_rows])),
        "steady_state_bias_mm": float(np.nanmean([row["steady_state_bias_mm"] for row in segment_metric_rows])),
        "IAE_mm_s": float(np.sum(np.abs(err_hist_arr)) * ts_s),
        "ISE_mm2_s": float(np.sum(np.square(err_hist_arr)) * ts_s),
        "ITAE_mm_s2": scenario_itae,
        "control_residual_rms_v": float(np.sqrt(np.mean(np.square(u_cmd_hist_arr - u_ff_hist_arr)))),
        "command_total_variation_v": float(np.sum(np.abs(np.diff(u_cmd_hist_arr)))) if len(u_cmd_hist_arr) > 1 else 0.0,
        "saturation_fraction": float(np.mean(sat_hist_arr)),
    }
    return scenario_metrics, True, ""


def _pid_validation_ok(
    scenario_metrics: Dict[str, Dict[str, float]],
    scenario_validity: Dict[str, Tuple[bool, str]],
    model: IdentifiedModel,
) -> Tuple[bool, str]:
    invalid_reasons = [reason for valid, reason in scenario_validity.values() if not valid and reason]
    if invalid_reasons:
        return False, "|".join(sorted(set(invalid_reasons)))
    if not model.stable_model:
        return False, "unstable_model"

    primary_metric_names = [
        "settled_error_p2p_mm",
        "ITAE_mm_s2",
        "overshoot_pct",
        "settling_time_s_2pct",
        "tracking_rms_mm",
    ]
    aggregate_metrics = scenario_metrics["all"]
    if any(not np.isfinite(float(aggregate_metrics[name])) for name in primary_metric_names):
        return False, "nonfinite_primary_metric"

    scenario_saturation = [
        float(metrics["saturation_fraction"])
        for name, metrics in scenario_metrics.items()
        if name != "all"
    ]
    if scenario_saturation and all(value >= 0.99 for value in scenario_saturation):
        return False, "saturation_lockup"

    return True, ""


def _pid_boundary_fields_for_result(result: PIDControllerResult, config: PIDDesignConfig) -> str:
    fields: List[str] = []
    if result.stage == "shared_seed":
        if result.kp_shared in (min(config.kp_grid), max(config.kp_grid)):
            fields.append("kp_shared")
        if result.ki_shared in (min(config.ki_grid), max(config.ki_grid)):
            fields.append("ki_shared")
        if result.kd_shared in (min(config.kd_grid), max(config.kd_grid)):
            fields.append("kd_shared")
    else:
        if result.kp_shared in (min(config.kp_grid), max(config.kp_grid)):
            fields.append("kp_seed")
        if result.ki_shared in (min(config.ki_grid), max(config.ki_grid)):
            fields.append("ki_seed")
        if result.kd_shared in (min(config.kd_grid), max(config.kd_grid)):
            fields.append("kd_seed")
        factor_fields = {
            "kp_up_factor": result.kp_up_factor,
            "ki_up_factor": result.ki_up_factor,
            "kd_up_factor": result.kd_up_factor,
            "kp_down_factor": result.kp_down_factor,
            "ki_down_factor": result.ki_down_factor,
            "kd_down_factor": result.kd_down_factor,
        }
        for name, value in factor_fields.items():
            if value is not None and value in (min(config.branch_gain_factors), max(config.branch_gain_factors)):
                fields.append(name)

    tau_grid = config.d_tau_grid_s + (config.shared_seed_tau_d_s,)
    if result.tau_d_s in (min(tau_grid), max(tau_grid)):
        fields.append("tau_d_s")
    slope_grid = config.d_slope_window_grid_samples + (config.shared_seed_slope_window_samples,)
    if result.slope_window_samples in (min(slope_grid), max(slope_grid)):
        fields.append("slope_window_samples")
    prefilter_grid = config.d_position_prefilter_tau_grid_s + (config.shared_seed_position_prefilter_tau_s,)
    if result.position_prefilter_tau_s in (min(prefilter_grid), max(prefilter_grid)):
        fields.append("position_prefilter_tau_s")
    if result.branch_deadband_mm in (
        min(config.branch_deadband_grid_mm),
        max(config.branch_deadband_grid_mm),
    ):
        fields.append("branch_deadband_mm")
    if result.position_deadband_mm in (
        min(config.position_deadband_grid_mm),
        max(config.position_deadband_grid_mm),
    ):
        fields.append("position_deadband_mm")
    return "|".join(fields)


def _build_pid_result(
    config_id: str,
    stage: str,
    kp_up: float,
    ki_up: float,
    kd_up: float,
    kp_down: float,
    ki_down: float,
    kd_down: float,
    tau_d_s: float,
    slope_window_samples: int,
    branch_deadband_mm: float,
    position_deadband_mm: float,
    kp_shared: float,
    ki_shared: float,
    kd_shared: float,
    kp_up_factor: float | None,
    ki_up_factor: float | None,
    kd_up_factor: float | None,
    kp_down_factor: float | None,
    ki_down_factor: float | None,
    kd_down_factor: float | None,
    aggregate_metrics: Dict[str, float],
    validation_ok: bool,
    invalid_reason: str,
    pid_config: PIDDesignConfig,
    position_prefilter_tau_s: float,
    derivative_filter_alpha: float,
    derivative_position_prefilter_alpha: float,
) -> PIDControllerResult:
    gain_l1_norm = _pid_gain_l1_norm(kp_up, ki_up, kd_up, kp_down, ki_down, kd_down)
    max_abs_gain_s = _pid_max_abs_gain_s(kp_up, kp_down)
    min_abs_integral_gain = _pid_min_abs_integral_gain(ki_up, ki_down)
    result = PIDControllerResult(
        config_id=config_id,
        stage=stage,
        variant="lookup_plus_branchwise_PID",
        controller_family="branchwise_pid_baseline",
        kp_up=float(kp_up),
        ki_up=float(ki_up),
        kd_up=float(kd_up),
        kp_down=float(kp_down),
        ki_down=float(ki_down),
        kd_down=float(kd_down),
        tau_d_s=float(tau_d_s),
        slope_window_samples=int(slope_window_samples),
        position_prefilter_tau_s=float(position_prefilter_tau_s),
        derivative_filter_alpha=float(derivative_filter_alpha),
        derivative_position_prefilter_alpha=float(derivative_position_prefilter_alpha),
        branch_deadband_mm=float(branch_deadband_mm),
        position_deadband_mm=float(position_deadband_mm),
        kp_shared=float(kp_shared),
        ki_shared=float(ki_shared),
        kd_shared=float(kd_shared),
        kp_up_factor=kp_up_factor,
        ki_up_factor=ki_up_factor,
        kd_up_factor=kd_up_factor,
        kp_down_factor=kp_down_factor,
        ki_down_factor=ki_down_factor,
        kd_down_factor=kd_down_factor,
        tracking_rms_mm=float(aggregate_metrics["tracking_rms_mm"]),
        tracking_overshoot_mm=float(aggregate_metrics["tracking_overshoot_mm"]),
        control_effort_mean_abs_v=float(aggregate_metrics["control_effort_mean_abs_v"]),
        settled_error_rms_mm=float(aggregate_metrics["settled_error_rms_mm"]),
        settled_error_p2p_mm=float(aggregate_metrics["settled_error_p2p_mm"]),
        settled_command_p2p_v=float(aggregate_metrics["settled_command_p2p_v"]),
        settled_velocity_rms_mm_s=float(aggregate_metrics["settled_velocity_rms_mm_s"]),
        rise_time_s_10_90=float(aggregate_metrics["rise_time_s_10_90"]),
        settling_time_s_2pct=float(aggregate_metrics["settling_time_s_2pct"]),
        peak_time_s=float(aggregate_metrics["peak_time_s"]),
        overshoot_mm=float(aggregate_metrics["overshoot_mm"]),
        overshoot_pct=float(aggregate_metrics["overshoot_pct"]),
        steady_state_bias_mm=float(aggregate_metrics["steady_state_bias_mm"]),
        abs_steady_state_bias_mm=abs(float(aggregate_metrics["steady_state_bias_mm"])),
        iae_mm_s=float(aggregate_metrics["IAE_mm_s"]),
        ise_mm2_s=float(aggregate_metrics["ISE_mm2_s"]),
        itae_mm_s2=float(aggregate_metrics["ITAE_mm_s2"]),
        control_residual_rms_v=float(aggregate_metrics["control_residual_rms_v"]),
        command_total_variation_v=float(aggregate_metrics["command_total_variation_v"]),
        saturation_fraction=float(aggregate_metrics["saturation_fraction"]),
        validation_ok=validation_ok,
        invalid_reason=invalid_reason,
        gain_l1_norm=gain_l1_norm,
        max_abs_gain_s=max_abs_gain_s,
        min_abs_integral_gain=min_abs_integral_gain,
        search_grid_boundary_hit=False,
        search_grid_boundary_fields="",
    )
    result.search_grid_boundary_fields = _pid_boundary_fields_for_result(result, pid_config)
    result.search_grid_boundary_hit = result.search_grid_boundary_fields != ""
    return result


def _evaluate_pid_candidate(
    *,
    stage: str,
    kp_up: float,
    ki_up: float,
    kd_up: float,
    kp_down: float,
    ki_down: float,
    kd_down: float,
    tau_d_s: float,
    slope_window_samples: int,
    branch_deadband_mm: float,
    position_deadband_mm: float,
    kp_shared: float,
    ki_shared: float,
    kd_shared: float,
    kp_up_factor: float | None,
    ki_up_factor: float | None,
    kd_up_factor: float | None,
    kp_down_factor: float | None,
    ki_down_factor: float | None,
    kd_down_factor: float | None,
    config_id: str,
    position_prefilter_tau_s: float,
) -> Tuple[PIDControllerResult, List[Dict[str, float | str | bool]]]:
    if (
        _PID_WORKER_MODEL is None
        or _PID_WORKER_REFERENCE_SCENARIOS is None
        or _PID_WORKER_BASE_LIMITS is None
        or _PID_WORKER_METRIC_CONFIG is None
        or _PID_WORKER_TS_S is None
        or _PID_WORKER_X_BAR_2 is None
        or _PID_WORKER_U_BAR is None
        or _PID_WORKER_SIM_CONTEXT is None
        or _PID_WORKER_PID_CONFIG is None
        or _PID_WORKER_TRANSIENT_PROFILES is None
    ):
        raise RuntimeError("PID worker context is not initialized.")

    runtime_limits = _pid_runtime_limits(_PID_WORKER_BASE_LIMITS, ki_up, ki_down)
    runtime_limits["branch_deadband_mm"] = float(branch_deadband_mm)
    runtime_limits["position_deadband_mm"] = float(position_deadband_mm)
    branch_gains = {
        "up": (float(kp_up), float(ki_up), float(kd_up)),
        "down": (float(kp_down), float(ki_down), float(kd_down)),
    }
    scenario_metrics: Dict[str, Dict[str, float]] = {}
    scenario_validity: Dict[str, Tuple[bool, str]] = {}
    for scenario_name, ref_values in _PID_WORKER_REFERENCE_SCENARIOS.items():
        metrics, valid, reason = _simulate_pid_controller(
            model=_PID_WORKER_MODEL,
            ref_mm=ref_values,
            ts_s=float(_PID_WORKER_TS_S),
            control_limits=runtime_limits,
            metric_config=_PID_WORKER_METRIC_CONFIG,
            branch_gains=branch_gains,
            tau_d_s=float(tau_d_s),
            slope_window_samples=int(slope_window_samples),
            position_prefilter_tau_s=float(position_prefilter_tau_s),
            x_bar_2=_PID_WORKER_X_BAR_2,
            u_bar=float(_PID_WORKER_U_BAR),
            sim_context=_PID_WORKER_SIM_CONTEXT,
            transient_profiles=_PID_WORKER_TRANSIENT_PROFILES,
        )
        scenario_metrics[scenario_name] = metrics
        scenario_validity[scenario_name] = (valid, reason)
    scenario_metrics["all"] = aggregate_scenario_metrics(scenario_metrics)
    validation_ok, invalid_reason = _pid_validation_ok(
        scenario_metrics=scenario_metrics,
        scenario_validity=scenario_validity,
        model=_PID_WORKER_MODEL,
    )
    result = _build_pid_result(
        config_id=config_id,
        stage=stage,
        kp_up=float(kp_up),
        ki_up=float(ki_up),
        kd_up=float(kd_up),
        kp_down=float(kp_down),
        ki_down=float(ki_down),
        kd_down=float(kd_down),
        tau_d_s=float(tau_d_s),
        slope_window_samples=int(slope_window_samples),
        branch_deadband_mm=float(branch_deadband_mm),
        position_deadband_mm=float(position_deadband_mm),
        kp_shared=float(kp_shared),
        ki_shared=float(ki_shared),
        kd_shared=float(kd_shared),
        kp_up_factor=kp_up_factor,
        ki_up_factor=ki_up_factor,
        kd_up_factor=kd_up_factor,
        kp_down_factor=kp_down_factor,
        ki_down_factor=ki_down_factor,
        kd_down_factor=kd_down_factor,
        aggregate_metrics=scenario_metrics["all"],
        validation_ok=validation_ok,
        invalid_reason=invalid_reason,
        pid_config=_PID_WORKER_PID_CONFIG,
        position_prefilter_tau_s=float(position_prefilter_tau_s),
        derivative_filter_alpha=alpha_from_tau(float(_PID_WORKER_TS_S), float(tau_d_s)),
        derivative_position_prefilter_alpha=alpha_from_tau(float(_PID_WORKER_TS_S), float(position_prefilter_tau_s)),
    )
    return result, _pid_sweep_rows_for_result(result=result, scenario_metrics=scenario_metrics)


def _evaluate_pid_seed_task(
    task: Tuple[float, float, float, float, float]
) -> Tuple[PIDControllerResult, List[Dict[str, float | str | bool]]]:
    kp_shared, ki_shared, kd_shared, branch_deadband_mm, position_deadband_mm = task
    if _PID_WORKER_PID_CONFIG is None:
        raise RuntimeError("PID worker config is not initialized.")
    tau_d_s = _PID_WORKER_PID_CONFIG.shared_seed_tau_d_s
    slope_window_samples = _PID_WORKER_PID_CONFIG.shared_seed_slope_window_samples
    position_prefilter_tau_s = _PID_WORKER_PID_CONFIG.shared_seed_position_prefilter_tau_s
    config_id = (
        f"pid_seed_kp{kp_shared:g}_ki{ki_shared:g}_kd{kd_shared:g}"
        f"_tau{tau_d_s:g}_n{slope_window_samples}_sp{position_prefilter_tau_s:g}"
        f"_bd{branch_deadband_mm:g}_pd{position_deadband_mm:g}"
    )
    return _evaluate_pid_candidate(
        stage="shared_seed",
        kp_up=float(kp_shared),
        ki_up=float(ki_shared),
        kd_up=float(kd_shared),
        kp_down=float(kp_shared),
        ki_down=float(ki_shared),
        kd_down=float(kd_shared),
        tau_d_s=float(tau_d_s),
        slope_window_samples=int(slope_window_samples),
        branch_deadband_mm=float(branch_deadband_mm),
        position_deadband_mm=float(position_deadband_mm),
        kp_shared=float(kp_shared),
        ki_shared=float(ki_shared),
        kd_shared=float(kd_shared),
        kp_up_factor=None,
        ki_up_factor=None,
        kd_up_factor=None,
        kp_down_factor=None,
        ki_down_factor=None,
        kd_down_factor=None,
        config_id=config_id,
        position_prefilter_tau_s=float(position_prefilter_tau_s),
    )


def _evaluate_pid_refinement_task(
    task: Tuple[PIDControllerResult, float, float, float]
) -> Tuple[PIDControllerResult, List[Dict[str, float | str | bool]]]:
    seed_result, tau_d_s, slope_window_samples, position_prefilter_tau_s = task
    config_id = (
        f"pid_smooth_seed{seed_result.config_id}"
        f"_tau{tau_d_s:g}_n{slope_window_samples}_sp{position_prefilter_tau_s:g}"
    )
    return _evaluate_pid_candidate(
        stage="smoothing_refinement",
        kp_up=float(seed_result.kp_shared),
        ki_up=float(seed_result.ki_shared),
        kd_up=float(seed_result.kd_shared),
        kp_down=float(seed_result.kp_shared),
        ki_down=float(seed_result.ki_shared),
        kd_down=float(seed_result.kd_shared),
        tau_d_s=float(tau_d_s),
        slope_window_samples=int(slope_window_samples),
        branch_deadband_mm=float(seed_result.branch_deadband_mm),
        position_deadband_mm=float(seed_result.position_deadband_mm),
        kp_shared=float(seed_result.kp_shared),
        ki_shared=float(seed_result.ki_shared),
        kd_shared=float(seed_result.kd_shared),
        kp_up_factor=None,
        ki_up_factor=None,
        kd_up_factor=None,
        kp_down_factor=None,
        ki_down_factor=None,
        kd_down_factor=None,
        config_id=config_id,
        position_prefilter_tau_s=float(position_prefilter_tau_s),
    )


def _evaluate_pid_branchwise_refinement_task(
    task: Tuple[PIDControllerResult, float, float, float, float, float, float]
) -> Tuple[PIDControllerResult, List[Dict[str, float | str | bool]]]:
    seed_result, kp_up_factor, ki_up_factor, kd_up_factor, kp_down_factor, ki_down_factor, kd_down_factor = task
    kp_up = seed_result.kp_shared * kp_up_factor
    ki_up = seed_result.ki_shared * ki_up_factor
    kd_up = seed_result.kd_shared * kd_up_factor
    kp_down = seed_result.kp_shared * kp_down_factor
    ki_down = seed_result.ki_shared * ki_down_factor
    kd_down = seed_result.kd_shared * kd_down_factor
    config_id = (
        f"pid_refine_seed{seed_result.config_id}"
        f"_ku{kp_up_factor:g}_iu{ki_up_factor:g}_du{kd_up_factor:g}"
        f"_kd{kp_down_factor:g}_id{ki_down_factor:g}_dd{kd_down_factor:g}"
    )
    return _evaluate_pid_candidate(
        stage="branchwise_refinement",
        kp_up=float(kp_up),
        ki_up=float(ki_up),
        kd_up=float(kd_up),
        kp_down=float(kp_down),
        ki_down=float(ki_down),
        kd_down=float(kd_down),
        tau_d_s=float(seed_result.tau_d_s),
        slope_window_samples=int(seed_result.slope_window_samples),
        branch_deadband_mm=float(seed_result.branch_deadband_mm),
        position_deadband_mm=float(seed_result.position_deadband_mm),
        kp_shared=float(seed_result.kp_shared),
        ki_shared=float(seed_result.ki_shared),
        kd_shared=float(seed_result.kd_shared),
        kp_up_factor=float(kp_up_factor),
        ki_up_factor=float(ki_up_factor),
        kd_up_factor=float(kd_up_factor),
        kp_down_factor=float(kp_down_factor),
        ki_down_factor=float(ki_down_factor),
        kd_down_factor=float(kd_down_factor),
        config_id=config_id,
        position_prefilter_tau_s=float(seed_result.position_prefilter_tau_s),
    )


def _should_print_pid_progress(completed: int, total: int, progress_step: int) -> bool:
    early_marks = {1, 2, 5, 10, 25, 50}
    return completed in early_marks or completed % progress_step == 0 or completed == total


def _collect_pid_results_with_progress(
    futures: List,
    label: str,
    total: int,
    progress_step: int,
    heartbeat_s: float = 10.0,
) -> Tuple[List[PIDControllerResult], List[List[Dict[str, float | str | bool]]]]:
    results: List[PIDControllerResult] = []
    sweep_row_groups: List[List[Dict[str, float | str | bool]]] = []
    pending = set(futures)
    completed = 0
    print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
    while pending:
        done, pending = wait(pending, timeout=heartbeat_s, return_when=FIRST_COMPLETED)
        if not done:
            print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
            continue
        for future in done:
            result, result_sweep_rows = future.result()
            results.append(result)
            sweep_row_groups.append(result_sweep_rows)
            completed += 1
            if _should_print_pid_progress(completed, total, progress_step):
                print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
    return results, sweep_row_groups


def evaluate_pid_candidates(
    model: IdentifiedModel,
    lookup_maps: Dict[str, LookupMap],
    ts_s: float,
    control_limits: Dict[str, float | str | bool],
    steady_df: pd.DataFrame,
    metric_config: ControllerDesignConfig,
    pid_config: PIDDesignConfig,
    transient_profiles: TransientProfiles,
) -> Tuple[PIDControllerResult, List[PIDControllerResult], List[Dict[str, float | str | bool]]]:
    lookup_tables = build_lookup_tables(lookup_maps)
    hold_steps = max(60, int(round(1.1 / ts_s)))
    reference_scenarios = build_reference_scenarios(lookup_tables, hold_steps)
    x_bar_2, u_bar = build_global_means(steady_df)
    sim_context = _build_pid_simulation_context(model, lookup_maps)

    all_results: List[PIDControllerResult] = []
    sweep_rows: List[Dict[str, float | str | bool]] = []
    shared_results: List[PIDControllerResult] = []

    ctx = mp.get_context("fork")
    seed_tasks = [
        (kp_shared, ki_shared, kd_shared, branch_deadband_mm, position_deadband_mm)
        for kp_shared in pid_config.kp_grid
        for ki_shared in pid_config.ki_grid
        for kd_shared in pid_config.kd_grid
        for branch_deadband_mm in pid_config.branch_deadband_grid_mm
        for position_deadband_mm in pid_config.position_deadband_grid_mm
    ]
    max_workers = max(1, (os.cpu_count() or 2) - 1)
    shared_progress_step = max(1, min(250, len(seed_tasks) // 10 if len(seed_tasks) > 0 else 1))
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_pid_worker,
        initargs=(
            model,
            reference_scenarios,
            control_limits,
            metric_config,
            ts_s,
            x_bar_2,
            u_bar,
            sim_context,
            pid_config,
            transient_profiles,
        ),
    ) as executor:
        futures = [executor.submit(_evaluate_pid_seed_task, task) for task in seed_tasks]
        shared_result_list, shared_sweep_groups = _collect_pid_results_with_progress(
            futures=futures,
            label="PID shared-seed progress",
            total=len(seed_tasks),
            progress_step=shared_progress_step,
        )
        for result, result_sweep_rows in zip(shared_result_list, shared_sweep_groups):
            shared_results.append(result)
            all_results.append(result)
            sweep_rows.extend(result_sweep_rows)

    ranked_shared, _, _ = rank_results_for_selection(shared_results)
    top_shared_results = ranked_shared[: pid_config.top_shared_pid_configs]

    smoothed_results: List[PIDControllerResult] = []
    smoothing_tasks = [
        (
            seed_result,
            tau_d_s,
            slope_window_samples,
            position_prefilter_tau_s,
        )
        for seed_result in top_shared_results
        for tau_d_s in pid_config.d_tau_grid_s
        for slope_window_samples in pid_config.d_slope_window_grid_samples
        for position_prefilter_tau_s in pid_config.d_position_prefilter_tau_grid_s
    ]
    smoothing_progress_step = max(1, min(250, len(smoothing_tasks) // 10 if len(smoothing_tasks) > 0 else 1))
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_pid_worker,
        initargs=(
            model,
            reference_scenarios,
            control_limits,
            metric_config,
            ts_s,
            x_bar_2,
            u_bar,
            sim_context,
            pid_config,
            transient_profiles,
        ),
    ) as executor:
        futures = [executor.submit(_evaluate_pid_refinement_task, task) for task in smoothing_tasks]
        smoothed_result_list, smoothed_sweep_groups = _collect_pid_results_with_progress(
            futures=futures,
            label="PID smoothing-refinement progress",
            total=len(smoothing_tasks),
            progress_step=smoothing_progress_step,
        )
        for result, result_sweep_rows in zip(smoothed_result_list, smoothed_sweep_groups):
            smoothed_results.append(result)
            all_results.append(result)
            sweep_rows.extend(result_sweep_rows)

    if not smoothed_results:
        raise ValueError("PID smoothing refinement did not produce any candidate results.")

    ranked_smoothed, _, _ = rank_results_for_selection(smoothed_results)
    top_smoothed_results = ranked_smoothed[: pid_config.top_shared_pid_configs]

    refined_results: List[PIDControllerResult] = []
    refinement_tasks = [
        (
            seed_result,
            kp_up_factor,
            ki_up_factor,
            kd_up_factor,
            kp_down_factor,
            ki_down_factor,
            kd_down_factor,
        )
        for seed_result in top_smoothed_results
        for kp_up_factor in pid_config.branch_gain_factors
        for ki_up_factor in pid_config.branch_gain_factors
        for kd_up_factor in pid_config.branch_gain_factors
        for kp_down_factor in pid_config.branch_gain_factors
        for ki_down_factor in pid_config.branch_gain_factors
        for kd_down_factor in pid_config.branch_gain_factors
    ]
    refinement_progress_step = max(1, min(250, len(refinement_tasks) // 10 if len(refinement_tasks) > 0 else 1))
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_pid_worker,
        initargs=(
            model,
            reference_scenarios,
            control_limits,
            metric_config,
            ts_s,
            x_bar_2,
            u_bar,
            sim_context,
            pid_config,
            transient_profiles,
        ),
    ) as executor:
        futures = [executor.submit(_evaluate_pid_branchwise_refinement_task, task) for task in refinement_tasks]
        refined_result_list, refined_sweep_groups = _collect_pid_results_with_progress(
            futures=futures,
            label="PID branchwise-refinement progress",
            total=len(refinement_tasks),
            progress_step=refinement_progress_step,
        )
        for result, result_sweep_rows in zip(refined_result_list, refined_sweep_groups):
            refined_results.append(result)
            all_results.append(result)
            sweep_rows.extend(result_sweep_rows)

    if not refined_results:
        raise ValueError("PID branchwise refinement did not produce any candidate results.")

    ranked_all, gate_used, best_tracking_rms = rank_results_for_selection(all_results)
    chosen_result = ranked_all[0]

    for result in all_results:
        result.selection_score_primary = float(result.selection_score_primary)
        result.selection_tracking_gate_passed = bool(result.selection_tracking_gate_passed)
        result.selection_eligibility_passed = bool(result.selection_eligibility_passed)

    chosen_result.selection_tracking_gate_passed = bool(chosen_result.selection_tracking_gate_passed)
    chosen_result.selection_eligibility_passed = bool(chosen_result.selection_eligibility_passed)
    if not np.isfinite(best_tracking_rms):
        raise ValueError("PID selection did not produce a finite best tracking RMS.")
    if gate_used not in (True, False):
        raise ValueError("PID selection gate state is invalid.")

    return chosen_result, all_results, sweep_rows


def build_pid_control_limits(
    control_window_row: pd.Series,
    excitation_summary_row: pd.Series,
    chosen_result: PIDControllerResult,
) -> Dict[str, float | str | bool]:
    excitation_low = float(excitation_summary_row["excitation_low_v"])
    excitation_high = float(excitation_summary_row["excitation_high_v"])
    ki_max = max(abs(chosen_result.ki_up), abs(chosen_result.ki_down))
    if ki_max < 1e-6:
        integral_limit = 10.0
    else:
        integral_limit = 0.40 * (excitation_high - excitation_low) / ki_max
    return {
        "u_min_v": float(control_window_row["window_low_v"]),
        "u_max_v": float(control_window_row["window_high_v"]),
        "excitation_low_v": excitation_low,
        "excitation_high_v": excitation_high,
        "integral_state_limit": float(integral_limit),
        "branch_selection_rule": "hold_previous_inside_deadband_then_sign_of_position_error",
        "branch_deadband_mm": float(chosen_result.branch_deadband_mm),
        "anti_windup_mode": "conditional_freeze",
        "position_deadband_mm": float(chosen_result.position_deadband_mm),
        "derivative_filter_tau_s": float(chosen_result.tau_d_s),
        "slope_window_samples": int(chosen_result.slope_window_samples),
        "derivative_filter_alpha": float(chosen_result.derivative_filter_alpha),
        "derivative_position_prefilter_tau_s": float(chosen_result.position_prefilter_tau_s),
        "derivative_position_prefilter_alpha": float(chosen_result.derivative_position_prefilter_alpha),
        "min_abs_integral_gain": float(chosen_result.min_abs_integral_gain),
        "uses_creep_compensation": chosen_result.uses_creep_compensation,
        "uses_dynamic_pressure_reference": chosen_result.uses_dynamic_pressure_reference,
        "transient_profile_horizon_s": chosen_result.transient_profile_horizon_s,
        "transient_lookup_source": chosen_result.transient_lookup_source,
        "uses_lookup_prefilter": True,
    }


def write_pid_feedback_gains(output_dir: Path, chosen_result: PIDControllerResult) -> None:
    rows = [
        {"branch": "up", "term": "kp_v_per_mm", "value": float(chosen_result.kp_up)},
        {"branch": "up", "term": "ki_v_per_mm_s", "value": float(chosen_result.ki_up)},
        {"branch": "up", "term": "kd_v_s_per_mm", "value": float(chosen_result.kd_up)},
        {"branch": "down", "term": "kp_v_per_mm", "value": float(chosen_result.kp_down)},
        {"branch": "down", "term": "ki_v_per_mm_s", "value": float(chosen_result.ki_down)},
        {"branch": "down", "term": "kd_v_s_per_mm", "value": float(chosen_result.kd_down)},
    ]
    pd.DataFrame(rows).to_csv(output_dir / "pid_feedback_gains.csv", index=False)


def write_pid_controller_limits(
    output_dir: Path,
    control_limits: Dict[str, float | str | bool],
) -> None:
    pd.DataFrame([control_limits]).to_csv(output_dir / "pid_controller_limits.csv", index=False)


def build_pid_design_summary(
    chosen_result: PIDControllerResult,
    all_results: List[PIDControllerResult],
    ts_s: float,
    steady_df: pd.DataFrame,
    data_df: pd.DataFrame,
    tail_fraction: float,
    input_dataset: str,
) -> pd.DataFrame:
    ranked_results, selection_gate_used, selection_best_tracking_rms_mm = rank_results_for_selection(all_results)
    ranked_ids = {result.config_id for result in ranked_results}
    row = {
        "variant": chosen_result.variant,
        "chosen_variant": True,
        "global_selected": True,
        "input_dataset": input_dataset,
        "data_split_policy": "single_dataset_no_hold_split",
        "controller_family": chosen_result.controller_family,
        "stage": chosen_result.stage,
        "sampling_time_s": float(ts_s),
        "n_samples": int(len(data_df)),
        "n_holds_total": int(data_df["hold_id"].nunique()),
        "n_holds_full": int((data_df.groupby("hold_id")["sample_id"].count() > 1).sum()),
        "steady_tail_fraction": float(tail_fraction),
        "n_steady_points_total": int(len(steady_df)),
        "n_steady_points_accepted": int(steady_df["steady_enough"].sum()),
        "selection_dataset_role": "single_dataset_simulation",
        "uses_lookup_prefilter": True,
        "uses_filtered_derivative": True,
        "uses_creep_compensation": chosen_result.uses_creep_compensation,
        "uses_dynamic_pressure_reference": chosen_result.uses_dynamic_pressure_reference,
        "transient_profile_horizon_s": chosen_result.transient_profile_horizon_s,
        "transient_lookup_source": chosen_result.transient_lookup_source,
        "kp_up_v_per_mm": float(chosen_result.kp_up),
        "ki_up_v_per_mm_s": float(chosen_result.ki_up),
        "kd_up_v_s_per_mm": float(chosen_result.kd_up),
        "kp_down_v_per_mm": float(chosen_result.kp_down),
        "ki_down_v_per_mm_s": float(chosen_result.ki_down),
        "kd_down_v_s_per_mm": float(chosen_result.kd_down),
        "tau_d_s": float(chosen_result.tau_d_s),
        "slope_window_samples": int(chosen_result.slope_window_samples),
        "position_prefilter_tau_s": float(chosen_result.position_prefilter_tau_s),
        "derivative_filter_alpha": float(chosen_result.derivative_filter_alpha),
        "derivative_position_prefilter_alpha": float(chosen_result.derivative_position_prefilter_alpha),
        "branch_deadband_mm": float(chosen_result.branch_deadband_mm),
        "position_deadband_mm": float(chosen_result.position_deadband_mm),
        "gain_l1_norm": float(chosen_result.gain_l1_norm),
        "max_abs_gain_s": float(chosen_result.max_abs_gain_s),
        "min_abs_integral_gain": float(chosen_result.min_abs_integral_gain),
        "tracking_rms_mm": float(chosen_result.tracking_rms_mm),
        "tracking_overshoot_mm": float(chosen_result.tracking_overshoot_mm),
        "control_effort_mean_abs_v": float(chosen_result.control_effort_mean_abs_v),
        "settled_error_rms_mm": float(chosen_result.settled_error_rms_mm),
        "settled_error_p2p_mm": float(chosen_result.settled_error_p2p_mm),
        "settled_command_p2p_v": float(chosen_result.settled_command_p2p_v),
        "settled_velocity_rms_mm_s": float(chosen_result.settled_velocity_rms_mm_s),
        "rise_time_s_10_90": float(chosen_result.rise_time_s_10_90),
        "settling_time_s_2pct": float(chosen_result.settling_time_s_2pct),
        "peak_time_s": float(chosen_result.peak_time_s),
        "overshoot_mm": float(chosen_result.overshoot_mm),
        "overshoot_pct": float(chosen_result.overshoot_pct),
        "steady_state_bias_mm": float(chosen_result.steady_state_bias_mm),
        "abs_steady_state_bias_mm": float(chosen_result.abs_steady_state_bias_mm),
        "IAE_mm_s": float(chosen_result.iae_mm_s),
        "ISE_mm2_s": float(chosen_result.ise_mm2_s),
        "ITAE_mm_s2": float(chosen_result.itae_mm_s2),
        "control_residual_rms_v": float(chosen_result.control_residual_rms_v),
        "command_total_variation_v": float(chosen_result.command_total_variation_v),
        "saturation_fraction": float(chosen_result.saturation_fraction),
        "validation_ok": bool(chosen_result.validation_ok),
        "invalid_reason": chosen_result.invalid_reason,
        "search_grid_boundary_hit": bool(chosen_result.search_grid_boundary_hit),
        "search_grid_boundary_fields": chosen_result.search_grid_boundary_fields,
        "selection_gate_tracking_factor": SELECTION_TRACKING_GATE_FACTOR,
        "selection_gain_penalty_exponent": SELECTION_GAIN_PENALTY_EXPONENT,
        "selection_score_formula": SELECTION_SCORE_FORMULA,
        "selection_primary_metrics": "|".join(SELECTION_PRIMARY_METRIC_EXPORT_NAMES),
        "selection_tie_break_metrics": "|".join(SELECTION_TIE_BREAK_EXPORT_NAMES),
        "selection_best_tracking_rms_mm": float(selection_best_tracking_rms_mm),
        "selection_tracking_gate_passed": bool(chosen_result.selection_tracking_gate_passed),
        "selection_eligibility_passed": bool(chosen_result.selection_eligibility_passed),
        "selection_gate_used": bool(selection_gate_used),
        "selection_score_primary": float(chosen_result.selection_score_primary),
        "selection_ranked_pool_member": chosen_result.config_id in ranked_ids,
    }
    return pd.DataFrame([row])


def build_pid_sweep_metrics(
    sweep_rows: List[Dict[str, float | str | bool]],
    all_results: List[PIDControllerResult],
    chosen_result: PIDControllerResult,
    input_dataset: str,
) -> pd.DataFrame:
    result_by_config = {result.config_id: result for result in all_results}
    rows: List[Dict[str, float | str | bool]] = []
    for row in sweep_rows:
        result = result_by_config[str(row["config_id"])]
        row_copy = dict(row)
        row_copy["chosen_variant"] = True
        row_copy["global_selected"] = result.config_id == chosen_result.config_id
        row_copy["input_dataset"] = input_dataset
        row_copy["data_split_policy"] = "single_dataset_no_hold_split"
        row_copy["selection_dataset_role"] = "single_dataset_simulation"
        row_copy["selection_gain_penalty_exponent"] = SELECTION_GAIN_PENALTY_EXPONENT
        row_copy["selection_score_formula"] = SELECTION_SCORE_FORMULA
        row_copy["selection_primary_metrics"] = "|".join(SELECTION_PRIMARY_METRIC_EXPORT_NAMES)
        row_copy["selection_tie_break_metrics"] = "|".join(SELECTION_TIE_BREAK_EXPORT_NAMES)
        row_copy["selection_score_primary"] = float(result.selection_score_primary)
        row_copy["selection_tracking_gate_passed"] = bool(result.selection_tracking_gate_passed)
        row_copy["selection_eligibility_passed"] = bool(result.selection_eligibility_passed)
        rows.append(row_copy)
    return pd.DataFrame(rows)


def _pid_sweep_rows_for_result(
    result: PIDControllerResult,
    scenario_metrics: Dict[str, Dict[str, float]],
) -> List[Dict[str, float | str | bool]]:
    rows: List[Dict[str, float | str | bool]] = []
    for scenario_name, metrics in scenario_metrics.items():
        rows.append(
            {
                "variant": result.variant,
                "controller_family": result.controller_family,
                "stage": result.stage,
                "config_id": result.config_id,
                "validation_ok": bool(result.validation_ok),
                "uses_creep_compensation": result.uses_creep_compensation,
                "uses_dynamic_pressure_reference": result.uses_dynamic_pressure_reference,
                "transient_profile_horizon_s": result.transient_profile_horizon_s,
                "transient_lookup_source": result.transient_lookup_source,
                "invalid_reason": result.invalid_reason,
                "kp_up_v_per_mm": float(result.kp_up),
                "ki_up_v_per_mm_s": float(result.ki_up),
                "kd_up_v_s_per_mm": float(result.kd_up),
                "kp_down_v_per_mm": float(result.kp_down),
                "ki_down_v_per_mm_s": float(result.ki_down),
                "kd_down_v_s_per_mm": float(result.kd_down),
                "kp_shared": float(result.kp_shared),
                "ki_shared": float(result.ki_shared),
                "kd_shared": float(result.kd_shared),
                "kp_up_factor": result.kp_up_factor,
                "ki_up_factor": result.ki_up_factor,
                "kd_up_factor": result.kd_up_factor,
                "kp_down_factor": result.kp_down_factor,
                "ki_down_factor": result.ki_down_factor,
                "kd_down_factor": result.kd_down_factor,
                "tau_d_s": float(result.tau_d_s),
                "slope_window_samples": int(result.slope_window_samples),
                "position_prefilter_tau_s": float(result.position_prefilter_tau_s),
                "derivative_filter_alpha": float(result.derivative_filter_alpha),
                "derivative_position_prefilter_alpha": float(result.derivative_position_prefilter_alpha),
                "branch_deadband_mm": float(result.branch_deadband_mm),
                "position_deadband_mm": float(result.position_deadband_mm),
                "gain_l1_norm": float(result.gain_l1_norm),
                "max_abs_gain_s": float(result.max_abs_gain_s),
                "min_abs_integral_gain": float(result.min_abs_integral_gain),
                "search_grid_boundary_hit": bool(result.search_grid_boundary_hit),
                "search_grid_boundary_fields": result.search_grid_boundary_fields,
                "scenario": scenario_name,
                "abs_steady_state_bias_mm": float(result.abs_steady_state_bias_mm),
                **metrics,
            }
        )
    return rows
