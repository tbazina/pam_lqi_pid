from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.linalg import solve, solve_discrete_are

SELECTION_TRACKING_GATE_FACTOR = 1.12
SELECTION_MIN_ABS_INTEGRAL_GAIN = 0.001
SELECTION_PRIMARY_METRIC_ATTRS = (
    "abs_steady_state_bias_mm",
    "settled_error_p2p_mm",
    "itae_mm_s2",
    "tracking_rms_mm",
    "overshoot_pct",
    "settling_time_s_2pct",
    "settled_command_p2p_v",
    "settled_velocity_rms_mm_s",
    "max_abs_gain_s",
)
SELECTION_PRIMARY_METRIC_EXPORT_NAMES = (
    "abs_steady_state_bias_mm",
    "settled_error_p2p_mm",
    "ITAE_mm_s2",
    "tracking_rms_mm",
    "overshoot_pct",
    "settling_time_s_2pct",
    "settled_command_p2p_v",
    "settled_velocity_rms_mm_s",
    "max_abs_gain_s",
)
SELECTION_PRIMARY_WEIGHTS = {
    "abs_steady_state_bias_mm": 0.18,
    "settled_error_p2p_mm": 0.10,
    "itae_mm_s2": 0.10,
    "tracking_rms_mm": 0.04,
    "overshoot_pct": 0.08,
    "settling_time_s_2pct": 0.07,
    "settled_command_p2p_v": 0.17,
    "settled_velocity_rms_mm_s": 0.15,
    "max_abs_gain_s": 0.20,
}
SELECTION_GAIN_PENALTY_EXPONENT = 3.0
SELECTION_NORMALIZATION_FLOORS = {
    # Prevent near-zero bias minima from dominating the full score when candidates
    # differ only by numerically tiny steady-state offsets.
    "abs_steady_state_bias_mm": 0.01,
}
SELECTION_TIE_BREAK_ATTRS = (
    "saturation_fraction",
    "command_total_variation_v",
    "gain_l1_norm",
)
SELECTION_TIE_BREAK_EXPORT_NAMES = (
    "saturation_fraction",
    "command_total_variation_v",
    "gain_l1_norm",
)
SELECTION_SCORE_FORMULA = (
    "0.18*norm_abs_steady_state_bias_mm + 0.10*norm_settled_error_p2p_mm + "
    "0.10*norm_ITAE_mm_s2 + 0.04*norm_tracking_rms_mm + 0.08*norm_overshoot_pct + "
    "0.07*norm_settling_time_s_2pct + 0.17*norm_settled_command_p2p_v + "
    "0.15*norm_settled_velocity_rms_mm_s + 0.20*(norm_max_abs_gain_s^3)"
)


@dataclass(frozen=True)
class ControllerDesignConfig:
    q_s_grid: Tuple[float, ...]
    q_p_grid: Tuple[float, ...]
    q_i_grid: Tuple[float, ...]
    r_grid: Tuple[float, ...]
    branch_deadband_grid_mm: Tuple[float, ...]
    position_deadband_grid_mm: Tuple[float, ...]
    vel_tau_grid_s: Tuple[float, ...]
    vel_slope_window_grid_samples: Tuple[int, ...]
    vel_position_prefilter_tau_grid_s: Tuple[float, ...]
    q_v_grid: Tuple[float, ...]
    top_base_configs_for_vel: int
    rise_time_bounds: Tuple[float, float]
    settling_band_frac: float
    settling_band_floor_mm: float
    step_metric_min_step_mm: float
    settled_window_frac: float


@dataclass
class LookupMap:
    branch: str
    domain_name: str
    range_name: str
    grid_x: np.ndarray
    grid_y: np.ndarray

    def evaluate(self, x: np.ndarray | float) -> np.ndarray:
        x_arr = np.asarray(x, dtype=float)
        x_clip = np.clip(x_arr, self.grid_x[0], self.grid_x[-1])
        return np.interp(x_clip, self.grid_x, self.grid_y)


@dataclass(frozen=True)
class RuntimeTransientProfile:
    time_s: np.ndarray
    values: np.ndarray
    horizon_s: float

    def evaluate(self, elapsed_s: float) -> float:
        return float(np.interp(float(elapsed_s), self.time_s, self.values))


@dataclass(frozen=True)
class TransientProfiles:
    creep_by_branch: Dict[str, RuntimeTransientProfile]
    pressure_by_branch: Dict[str, RuntimeTransientProfile]
    source_lookup: str = "open_loop_steady_state_lookup"


def _read_runtime_profile(
    path: Path,
    value_column: str,
    profile_name: str,
) -> RuntimeTransientProfile:
    if not path.exists():
        raise FileNotFoundError(f"Missing {profile_name} profile: {path}")
    df = pd.read_csv(path)
    required = {"time_since_step_s", value_column, "runtime_horizon_s"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{profile_name} profile {path} is missing columns: {sorted(missing)}")
    time_s = pd.to_numeric(df["time_since_step_s"], errors="coerce").to_numpy(dtype=float)
    values = pd.to_numeric(df[value_column], errors="coerce").to_numpy(dtype=float)
    horizons = pd.to_numeric(df["runtime_horizon_s"], errors="coerce").to_numpy(dtype=float)
    if len(time_s) < 2 or np.any(~np.isfinite(time_s)) or np.any(~np.isfinite(values)):
        raise ValueError(f"{profile_name} profile {path} contains invalid values.")
    if np.any(np.diff(time_s) <= 0) or not np.allclose(np.diff(time_s), 0.01, atol=1e-9):
        raise ValueError(f"{profile_name} profile {path} must be strictly increasing at 10 ms spacing.")
    if not np.allclose(time_s[0], 0.0, atol=1e-9) or not np.allclose(time_s[-1], 3.99, atol=1e-9):
        raise ValueError(f"{profile_name} profile {path} must cover 0.00 to 3.99 s.")
    horizon_s = float(horizons[0])
    if not np.isfinite(horizon_s) or not np.allclose(horizons, horizon_s, atol=1e-9):
        raise ValueError(f"{profile_name} profile {path} has inconsistent runtime horizons.")
    return RuntimeTransientProfile(time_s=time_s, values=values, horizon_s=horizon_s)


def load_transient_profiles(
    output_dir: Path,
    use_creep_compensation: bool = True,
    use_dynamic_pressure_reference: bool = True,
) -> TransientProfiles:
    creep_by_branch: Dict[str, RuntimeTransientProfile] = {}
    pressure_by_branch: Dict[str, RuntimeTransientProfile] = {}
    if use_creep_compensation:
        for branch in ("up", "down"):
            creep_by_branch[branch] = _read_runtime_profile(
                output_dir / f"open_loop_creep_compensation_{branch}.csv",
                "creep_fraction",
                f"{branch} creep",
            )
    if use_dynamic_pressure_reference:
        for branch in ("up", "down"):
            pressure_by_branch[branch] = _read_runtime_profile(
                output_dir / f"open_loop_dynamic_pressure_reference_{branch}.csv",
                "pressure_progress_fraction",
                f"{branch} dynamic-pressure",
            )
    return TransientProfiles(
        creep_by_branch=creep_by_branch,
        pressure_by_branch=pressure_by_branch,
    )


@dataclass
class IdentifiedModel:
    variant: str
    controller_family: str
    state_names: List[str]
    runtime_observable_names: List[str]
    A: np.ndarray
    B: np.ndarray
    x_mean: np.ndarray
    u_mean: float
    uses_lookup_prefilter: bool
    stable_model: bool
    spectral_radius: float
    branch_models: Dict[str, Tuple[np.ndarray, np.ndarray]]
    uses_velocity_state: bool
    velocity_filter_tau_s: float | None
    velocity_slope_window_samples: int | None
    velocity_position_prefilter_tau_s: float | None


@dataclass
class ControllerResult:
    config_id: str
    variant: str
    model: IdentifiedModel
    evaluation_model: IdentifiedModel
    K_aug: np.ndarray
    A_aug: np.ndarray
    B_aug: np.ndarray
    C_track: np.ndarray
    q_diag: np.ndarray
    q_s: float
    q_p: float
    q_i: float
    q_v: float | None
    r_value: float
    branch_deadband_mm: float
    position_deadband_mm: float
    one_step_rmse_s_mm: float
    one_step_rmse_p_bar: float
    ten_step_rmse_s_mm: float
    ten_step_rmse_p_bar: float
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
    integral_state_limit: float
    validation_ok: bool
    branch_K_aug: Dict[str, np.ndarray]
    branch_A_aug: Dict[str, np.ndarray]
    branch_B_aug: Dict[str, np.ndarray]
    gain_l1_norm: float
    max_abs_gain_s: float
    max_abs_gain_p: float
    min_abs_integral_gain: float
    uses_velocity_state: bool
    velocity_filter_tau_s: float | None
    velocity_slope_window_samples: int | None
    velocity_position_prefilter_tau_s: float | None
    velocity_filter_alpha: float | None
    velocity_position_prefilter_alpha: float | None
    uses_creep_compensation: bool
    uses_dynamic_pressure_reference: bool
    transient_profile_horizon_s: float | None
    transient_lookup_source: str
    uses_pressure_filter: bool
    pressure_filter_tau_s: float | None
    scenario_metrics: Dict[str, Dict[str, float]]
    search_grid_boundary_hit: bool
    search_grid_boundary_fields: str
    selection_score_primary: float = float("inf")
    selection_tracking_gate_passed: bool = False
    selection_eligibility_passed: bool = False
    evaluation_model_stable: bool = True
    evaluation_spectral_radius: float = float("nan")


_LQI_WORKER_MODEL: IdentifiedModel | None = None
_LQI_WORKER_EVAL_MODEL: IdentifiedModel | None = None
_LQI_WORKER_LOOKUP_MAPS: Dict[str, LookupMap] | None = None
_LQI_WORKER_REFERENCE_SCENARIOS: Dict[str, np.ndarray] | None = None
_LQI_WORKER_TS_S: float | None = None
_LQI_WORKER_EVAL_TS_S: float | None = None
_LQI_WORKER_CONTROL_LIMITS: Dict[str, float | str | bool] | None = None
_LQI_WORKER_STEADY_DF: pd.DataFrame | None = None
_LQI_WORKER_CONFIG: ControllerDesignConfig | None = None
_LQI_WORKER_PREDICTION_METRICS: Dict[str, float] | None = None
_LQI_WORKER_TRANSIENT_PROFILES: TransientProfiles | None = None


def load_excitation_experiment(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [str(col).strip() for col in df.columns]

    rename_map = {
        "u[V]": "u_v",
        "s[mm]": "s_mm",
        "p[bar]": "p_bar",
        "t[ms]": "t_ms",
    }
    missing_cols = [src for src in rename_map if src not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df.rename(columns=rename_map)
    for col in rename_map.values():
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["u_v", "s_mm", "p_bar", "t_ms"]).reset_index(drop=True)
    df["sample_id"] = np.arange(1, len(df) + 1, dtype=int)
    return df


def infer_sampling_time_s(df: pd.DataFrame) -> float:
    valid_t_ms = df.loc[(df["sample_id"] > 1) & (df["t_ms"] < 1000), "t_ms"]
    if valid_t_ms.empty:
        raise ValueError("Could not infer sampling time from excitation experiment.")
    return float(np.median(valid_t_ms.to_numpy(dtype=float)) / 1000.0)


def segment_holds(df: pd.DataFrame) -> pd.DataFrame:
    u_values = df["u_v"].to_numpy(dtype=float)
    change = np.r_[True, np.abs(np.diff(u_values)) > 1e-12]
    hold_id = np.cumsum(change)
    df = df.copy()
    df["hold_id"] = hold_id

    hold_summary = (
        df.groupby("hold_id", sort=True)
        .agg(
            start_idx=("sample_id", "min"),
            end_idx=("sample_id", "max"),
            n_samples=("sample_id", "count"),
            u_cmd_v=("u_v", "first"),
        )
        .reset_index()
    )
    hold_summary["u_prev_v"] = hold_summary["u_cmd_v"].shift(1)
    hold_summary["u_prev_v"] = hold_summary["u_prev_v"].fillna(hold_summary["u_cmd_v"])
    hold_summary["du_v"] = hold_summary["u_cmd_v"] - hold_summary["u_prev_v"]
    hold_summary["branch"] = np.where(
        hold_summary["du_v"] > 1e-12,
        "up",
        np.where(hold_summary["du_v"] < -1e-12, "down", "start"),
    )

    return df.merge(
        hold_summary[["hold_id", "n_samples", "u_prev_v", "du_v", "branch"]],
        on="hold_id",
        how="left",
    )


def compare_to_excitation_schedule(
    df: pd.DataFrame, schedule_path: Path | None
) -> Dict[str, float | int | bool]:
    result = {
        "schedule_available": False,
        "schedule_match_holds": False,
        "schedule_match_commands": False,
        "schedule_holds": np.nan,
        "experiment_holds": int(df["hold_id"].nunique()),
    }
    if schedule_path is None or not schedule_path.exists():
        return result

    schedule_df = pd.read_csv(schedule_path)
    schedule_holds = schedule_df["hold_id"].nunique()
    hold_summary = (
        df.groupby("hold_id", sort=True)
        .agg(u_cmd_v=("u_v", "first"), n_samples=("sample_id", "count"))
        .reset_index()
    )
    hold_cmds = hold_summary.loc[hold_summary["n_samples"] > 1, "u_cmd_v"].round(3).to_numpy(dtype=float)
    schedule_cmds = schedule_df.sort_values("hold_id")["u_cmd_v"].round(3).to_numpy(dtype=float)

    same_count = len(hold_cmds) == len(schedule_cmds)
    same_cmds = same_count and np.allclose(hold_cmds, schedule_cmds, atol=1e-3)
    result.update(
        {
            "schedule_available": True,
            "schedule_match_holds": same_count,
            "schedule_match_commands": same_cmds,
            "schedule_holds": int(schedule_holds),
        }
    )
    return result


def extract_steady_state_points(
    df: pd.DataFrame, ts_s: float, tail_fraction: float = 0.25
) -> pd.DataFrame:
    if not (0 < tail_fraction <= 1):
        raise ValueError("tail_fraction must be in (0, 1].")

    rows: List[Dict[str, float | int | str | bool]] = []
    for hold_id, hold_df in df.groupby("hold_id", sort=True):
        hold_df = hold_df.sort_values("sample_id")
        n_samples = int(len(hold_df))
        if n_samples <= 1:
            continue

        tail_n = max(5, int(np.ceil(tail_fraction * n_samples)))
        tail_df = hold_df.iloc[-tail_n:].copy()
        dt_s = float(tail_df["t_ms"].sum() / 1000.0)
        if dt_s <= 0:
            continue

        s_vals = tail_df["s_mm"].to_numpy(dtype=float)
        p_vals = tail_df["p_bar"].to_numpy(dtype=float)
        s_slope = float((s_vals[-1] - s_vals[0]) / dt_s)
        p_slope = float((p_vals[-1] - p_vals[0]) / dt_s)
        s_std = float(np.std(s_vals, ddof=0))
        p_std = float(np.std(p_vals, ddof=0))
        steady_enough = (
            abs(s_slope) <= 1.2
            and abs(p_slope) <= 0.25
            and s_std <= 0.6
            and p_std <= 0.12
        )

        first_row = hold_df.iloc[0]
        rows.append(
            {
                "hold_id": int(hold_id),
                "u_ss_v": float(first_row["u_v"]),
                "u_prev_v": float(first_row["u_prev_v"]),
                "du_v": float(first_row["du_v"]),
                "branch": str(first_row["branch"]),
                "n_samples": n_samples,
                "tail_fraction_used": float(tail_fraction),
                "tail_n_samples": tail_n,
                "s_ss_mm": float(np.mean(s_vals)),
                "p_ss_bar": float(np.mean(p_vals)),
                "s_tail_std_mm": s_std,
                "p_tail_std_bar": p_std,
                "s_tail_slope_mm_s": s_slope,
                "p_tail_slope_bar_s": p_slope,
                "steady_enough": steady_enough,
            }
        )

    steady_df = pd.DataFrame(rows)
    if steady_df.empty:
        raise ValueError("No steady-state hold points were extracted.")
    return steady_df


def _moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1 or len(values) < 3:
        return values.copy()
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(values, kernel, mode="same")


def fit_monotone_lookup(
    data_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    branch: str,
    domain_name: str,
    range_name: str,
    n_bins: int = 60,
    grid_n: int = 200,
    smooth_window: int = 7,
) -> LookupMap:
    branch_df = data_df.loc[data_df["branch"] == branch, [x_col, y_col]].dropna()
    if len(branch_df) < 8:
        raise ValueError(f"Not enough steady points for branch {branch} and map {x_col}->{y_col}.")

    branch_df = branch_df.sort_values(x_col)
    x_vals = branch_df[x_col].to_numpy(dtype=float)
    y_vals = branch_df[y_col].to_numpy(dtype=float)

    edges = np.linspace(x_vals.min(), x_vals.max(), min(n_bins, len(branch_df)) + 1)
    bin_ids = np.digitize(x_vals, edges[1:-1], right=False)

    grouped_x: List[float] = []
    grouped_y: List[float] = []
    for bin_id in range(bin_ids.min(), bin_ids.max() + 1):
        mask = bin_ids == bin_id
        if np.any(mask):
            grouped_x.append(float(np.mean(x_vals[mask])))
            grouped_y.append(float(np.median(y_vals[mask])))

    x_unique = np.asarray(grouped_x, dtype=float)
    y_unique = np.asarray(grouped_y, dtype=float)
    order = np.argsort(x_unique)
    x_unique = x_unique[order]
    y_unique = y_unique[order]

    grid_x = np.linspace(x_unique[0], x_unique[-1], grid_n)
    grid_y = np.interp(grid_x, x_unique, y_unique)
    grid_y = _moving_average(grid_y, min(smooth_window, max(1, grid_n // 20 * 2 + 1)))
    grid_y = np.maximum.accumulate(grid_y)
    return LookupMap(
        branch=branch,
        domain_name=domain_name,
        range_name=range_name,
        grid_x=grid_x,
        grid_y=grid_y,
    )


def build_lookup_maps(steady_df: pd.DataFrame) -> Dict[str, LookupMap]:
    accepted_df = steady_df.loc[steady_df["steady_enough"] & steady_df["branch"].isin(["up", "down"])].copy()
    if accepted_df.empty:
        raise ValueError("No steady-state points satisfied the acceptance thresholds.")

    maps: Dict[str, LookupMap] = {}
    for branch in ("up", "down"):
        branch_df = accepted_df.loc[accepted_df["branch"] == branch].copy()
        if len(branch_df) < 8:
            raise ValueError(f"Branch {branch} has too few accepted steady points.")
        maps[f"{branch}_inverse_u"] = fit_monotone_lookup(
            branch_df, "s_ss_mm", "u_ss_v", branch, "s_ref_mm", "u_ff_v"
        )
        maps[f"{branch}_inverse_p"] = fit_monotone_lookup(
            branch_df, "s_ss_mm", "p_ss_bar", branch, "s_ref_mm", "p_ff_bar"
        )
    return maps


def build_global_means(steady_df: pd.DataFrame) -> Tuple[np.ndarray, float]:
    accepted_df = steady_df.loc[steady_df["steady_enough"] & steady_df["branch"].isin(["up", "down"])]
    s_bar = float(accepted_df["s_ss_mm"].mean())
    p_bar = float(accepted_df["p_ss_bar"].mean())
    u_bar = float(accepted_df["u_ss_v"].mean())
    return np.array([s_bar, p_bar], dtype=float), u_bar


def build_lookup_tables(lookup_maps: Dict[str, LookupMap]) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}
    for branch in ("up", "down"):
        tables[branch] = pd.DataFrame(
            {
                "branch": branch,
                "s_ref_mm": lookup_maps[f"{branch}_inverse_u"].grid_x,
                "u_ff_v": lookup_maps[f"{branch}_inverse_u"].grid_y,
                "p_ff_bar": lookup_maps[f"{branch}_inverse_p"].grid_y,
            }
        )
    return tables


def build_feedforward_alignment_outputs(
    steady_df: pd.DataFrame,
    lookup_maps: Dict[str, LookupMap],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    accepted_df = steady_df.loc[steady_df["steady_enough"] & steady_df["branch"].isin(["up", "down"])].copy()
    point_frames: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, float | int | str]] = []

    for branch in ("up", "down"):
        branch_df = accepted_df.loc[accepted_df["branch"] == branch].copy().reset_index(drop=True)
        u_ff = lookup_maps[f"{branch}_inverse_u"].evaluate(branch_df["s_ss_mm"].to_numpy(dtype=float))
        p_ff = lookup_maps[f"{branch}_inverse_p"].evaluate(branch_df["s_ss_mm"].to_numpy(dtype=float))
        merged = branch_df.assign(
            u_ff_v=u_ff,
            p_ff_bar=p_ff,
            u_abs_error_v=np.abs(branch_df["u_ss_v"].to_numpy(dtype=float) - u_ff),
            p_abs_error_bar=np.abs(branch_df["p_ss_bar"].to_numpy(dtype=float) - p_ff),
        )
        point_frames.append(merged)
        summary_rows.append(_alignment_summary_row("branch", branch, merged))

    points_df = pd.concat(point_frames, ignore_index=True)
    summary_rows.append(_alignment_summary_row("global", "global", points_df))
    return points_df, pd.DataFrame(summary_rows)


def _alignment_summary_row(level: str, branch: str, metric_df: pd.DataFrame) -> Dict[str, float | int | str]:
    def rmse(col: str) -> float:
        vals = metric_df[col].to_numpy(dtype=float)
        return float(np.sqrt(np.mean(np.square(vals))))

    def mae(col: str) -> float:
        return float(np.mean(metric_df[col].to_numpy(dtype=float)))

    def p95(col: str) -> float:
        return float(np.quantile(metric_df[col].to_numpy(dtype=float), 0.95))

    return {
        "level": level,
        "branch": branch,
        "n_points": int(len(metric_df)),
        "s_min_mm": float(metric_df["s_ss_mm"].min()),
        "s_max_mm": float(metric_df["s_ss_mm"].max()),
        "u_rmse_v": rmse("u_abs_error_v"),
        "u_mae_v": mae("u_abs_error_v"),
        "u_p95_abs_v": p95("u_abs_error_v"),
        "p_rmse_bar": rmse("p_abs_error_bar"),
        "p_mae_bar": mae("p_abs_error_bar"),
        "p_p95_abs_bar": p95("p_abs_error_bar"),
    }


def fit_linear_state_space_model(
    X: np.ndarray,
    U: np.ndarray,
    Y: np.ndarray,
    ridge_lambda: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray]:
    reg = np.hstack([X, U])
    gram = reg.T @ reg
    theta = np.linalg.solve(
        gram + ridge_lambda * np.eye(gram.shape[0], dtype=float),
        reg.T @ Y,
    )
    x_dim = X.shape[1]
    A = theta[:x_dim, :].T
    B = theta[x_dim:, :].T
    return A, B


def compute_sliding_linear_fit_slope(
    s_mm: np.ndarray,
    ts_s: float,
    window_samples: int,
) -> np.ndarray:
    if window_samples < 2:
        raise ValueError("window_samples must be at least 2.")
    s_arr = np.asarray(s_mm, dtype=float)
    slopes = np.zeros(len(s_arr), dtype=float)
    for idx in range(len(s_arr)):
        start = max(0, idx - window_samples + 1)
        segment = s_arr[start : idx + 1]
        if len(segment) < 2:
            slopes[idx] = 0.0
            continue
        t_seg = np.arange(len(segment), dtype=float) * ts_s
        t_centered = t_seg - float(np.mean(t_seg))
        denom = float(np.dot(t_centered, t_centered))
        if denom <= 0.0:
            slopes[idx] = 0.0
            continue
        s_centered = segment - float(np.mean(segment))
        slopes[idx] = float(np.dot(t_centered, s_centered) / denom)
    return slopes


def compute_first_order_lowpass_signal(
    x: np.ndarray,
    ts_s: float,
    tau_s: float,
) -> np.ndarray:
    x_arr = np.asarray(x, dtype=float)
    if x_arr.size == 0:
        return np.zeros(0, dtype=float)
    if tau_s <= 0.0:
        return x_arr.copy()
    alpha = _alpha_from_tau(ts_s, tau_s)
    y = np.empty_like(x_arr, dtype=float)
    y[0] = float(x_arr[0])
    for idx in range(1, len(x_arr)):
        y[idx] = alpha * y[idx - 1] + (1.0 - alpha) * float(x_arr[idx])
    return y


def update_first_order_lowpass(
    prev_value: float,
    current_value: float,
    alpha: float,
) -> float:
    return float(alpha * prev_value + (1.0 - alpha) * current_value)


def compute_latest_linear_fit_slope(
    s_history: Sequence[float] | np.ndarray,
    ts_s: float,
    window_samples: int,
) -> float:
    if window_samples < 2:
        raise ValueError("window_samples must be at least 2.")
    segment = np.asarray(s_history, dtype=float)
    if segment.size > window_samples:
        segment = segment[-window_samples:]
    if segment.size < 2:
        return 0.0
    t_seg = np.arange(segment.size, dtype=float) * ts_s
    t_centered = t_seg - float(np.mean(t_seg))
    denom = float(np.dot(t_centered, t_centered))
    if denom <= 0.0:
        return 0.0
    s_centered = segment - float(np.mean(segment))
    return float(np.dot(t_centered, s_centered) / denom)


def compute_filtered_velocity_signal(
    s_mm: np.ndarray,
    ts_s: float,
    tau_s: float,
    window_samples: int,
    position_prefilter_tau_s: float = 0.0,
) -> np.ndarray:
    alpha = _alpha_from_tau(ts_s, tau_s)
    s_for_slope = compute_first_order_lowpass_signal(s_mm, ts_s, position_prefilter_tau_s)
    v_slope_raw = compute_sliding_linear_fit_slope(s_for_slope, ts_s, window_samples)
    v_hat = np.zeros(len(s_mm), dtype=float)
    for idx in range(1, len(s_mm)):
        v_hat[idx] = alpha * v_hat[idx - 1] + (1.0 - alpha) * float(v_slope_raw[idx])
    return v_hat


def _build_state_arrays(
    df: pd.DataFrame,
    steady_df: pd.DataFrame,
    ts_s: float,
    state_names: Sequence[str],
    velocity_filter_tau_s: float | None = None,
    velocity_slope_window_samples: int | None = None,
    velocity_position_prefilter_tau_s: float | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    x_bar, u_bar = build_global_means(steady_df)
    state_columns: List[np.ndarray] = []
    for name in state_names:
        if name == "s_tilde":
            state_columns.append(df["s_mm"].to_numpy(dtype=float) - x_bar[0])
        elif name == "p_tilde":
            state_columns.append(df["p_bar"].to_numpy(dtype=float) - x_bar[1])
        elif name == "v_hat_mm_s":
            if velocity_filter_tau_s is None or velocity_slope_window_samples is None:
                raise ValueError("velocity filter tau and slope window are required for the velocity state.")
            state_columns.append(
                compute_filtered_velocity_signal(
                    df["s_mm"].to_numpy(dtype=float),
                    ts_s,
                    velocity_filter_tau_s,
                    velocity_slope_window_samples,
                    0.0 if velocity_position_prefilter_tau_s is None else velocity_position_prefilter_tau_s,
                )
            )
        else:
            raise ValueError(f"Unsupported state name: {name}")
    state_matrix = np.column_stack(state_columns)
    u_arr = df["u_v"].to_numpy(dtype=float) - u_bar
    return state_matrix, u_arr


def _identify_model_from_states(
    df: pd.DataFrame,
    state_matrix: np.ndarray,
    u_arr: np.ndarray,
    train_hold_ids: Iterable[int],
    variant: str,
    controller_family: str,
    state_names: List[str],
    runtime_observable_names: List[str],
    u_mean: float,
    uses_velocity_state: bool,
    velocity_filter_tau_s: float | None,
    velocity_slope_window_samples: int | None,
    velocity_position_prefilter_tau_s: float | None,
) -> IdentifiedModel:
    hold_ids_arr = np.asarray(list(train_hold_ids), dtype=int)
    hold_id_arr = df["hold_id"].to_numpy(dtype=int)
    branch_arr = df["branch"].to_numpy(dtype=str)
    idx_all = np.arange(0, len(df) - 1, dtype=int)

    branch_models: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    branch_radii: List[float] = []
    for branch in ("up", "down"):
        idx = idx_all[np.isin(hold_id_arr[idx_all], hold_ids_arr) & (branch_arr[idx_all] == branch)]
        if len(idx) < 20:
            raise ValueError(f"Not enough training samples for branch {branch}.")
        X = state_matrix[idx]
        U = u_arr[idx][:, None]
        Y = state_matrix[idx + 1]
        A_branch, B_branch = fit_linear_state_space_model(X, U, Y)
        branch_models[branch] = (A_branch, B_branch)
        branch_radii.append(float(np.max(np.abs(np.linalg.eigvals(A_branch)))))

    A_nom = 0.5 * (branch_models["up"][0] + branch_models["down"][0])
    B_nom = 0.5 * (branch_models["up"][1] + branch_models["down"][1])
    spectral_radius = float(max(branch_radii))

    return IdentifiedModel(
        variant=variant,
        controller_family=controller_family,
        state_names=list(state_names),
        runtime_observable_names=list(runtime_observable_names),
        A=A_nom,
        B=B_nom,
        x_mean=np.zeros(len(state_names), dtype=float),
        u_mean=u_mean,
        uses_lookup_prefilter=True,
        stable_model=spectral_radius < 1.0,
        spectral_radius=spectral_radius,
        branch_models=branch_models,
        uses_velocity_state=uses_velocity_state,
        velocity_filter_tau_s=velocity_filter_tau_s,
        velocity_slope_window_samples=velocity_slope_window_samples,
        velocity_position_prefilter_tau_s=velocity_position_prefilter_tau_s,
    )


def identify_candidate_models(
    df: pd.DataFrame,
    steady_df: pd.DataFrame,
    train_hold_ids: Iterable[int],
    config: ControllerDesignConfig,
) -> Dict[str, IdentifiedModel]:
    del config
    _, u_bar = build_global_means(steady_df)
    state_matrix, u_arr = _build_state_arrays(
        df=df,
        steady_df=steady_df,
        ts_s=1.0,
        state_names=["s_tilde", "p_tilde"],
    )
    base_model = _identify_model_from_states(
        df=df,
        state_matrix=state_matrix,
        u_arr=u_arr,
        train_hold_ids=train_hold_ids,
        variant="lookup_plus_branchwise_minimal_LQI",
        controller_family="minimal_branchwise_lqi",
        state_names=["s_tilde", "p_tilde"],
        runtime_observable_names=["s_tilde", "p_tilde"],
        u_mean=u_bar,
        uses_velocity_state=False,
        velocity_filter_tau_s=None,
        velocity_slope_window_samples=None,
        velocity_position_prefilter_tau_s=None,
    )
    return {base_model.variant: base_model}


def identify_velocity_model(
    df: pd.DataFrame,
    steady_df: pd.DataFrame,
    train_hold_ids: Iterable[int],
    ts_s: float,
    tau_s: float,
    slope_window_samples: int,
    position_prefilter_tau_s: float,
) -> IdentifiedModel:
    _, u_bar = build_global_means(steady_df)
    state_matrix, u_arr = _build_state_arrays(
        df=df,
        steady_df=steady_df,
        ts_s=ts_s,
        state_names=["s_tilde", "p_tilde", "v_hat_mm_s"],
        velocity_filter_tau_s=tau_s,
        velocity_slope_window_samples=slope_window_samples,
        velocity_position_prefilter_tau_s=position_prefilter_tau_s,
    )
    return _identify_model_from_states(
        df=df,
        state_matrix=state_matrix,
        u_arr=u_arr,
        train_hold_ids=train_hold_ids,
        variant="lookup_plus_branchwise_minimal_LQI_VEL",
        controller_family="minimal_branchwise_lqi_vel",
        state_names=["s_tilde", "p_tilde", "v_hat_mm_s"],
        runtime_observable_names=["s_tilde", "p_tilde", "v_hat_mm_s"],
        u_mean=u_bar,
        uses_velocity_state=True,
        velocity_filter_tau_s=tau_s,
        velocity_slope_window_samples=slope_window_samples,
        velocity_position_prefilter_tau_s=position_prefilter_tau_s,
    )


def evaluate_prediction_metrics(
    model: IdentifiedModel,
    df: pd.DataFrame,
    val_hold_ids: Iterable[int] | None,
    steady_df: pd.DataFrame,
    ts_s: float,
) -> Dict[str, float]:
    state_matrix, u_arr = _build_state_arrays(
        df=df,
        steady_df=steady_df,
        ts_s=ts_s,
        state_names=model.state_names,
        velocity_filter_tau_s=model.velocity_filter_tau_s,
        velocity_slope_window_samples=model.velocity_slope_window_samples,
        velocity_position_prefilter_tau_s=model.velocity_position_prefilter_tau_s,
    )
    hold_id_arr = df["hold_id"].to_numpy(dtype=int)
    branch_arr = df["branch"].to_numpy(dtype=str)
    if val_hold_ids is None:
        val_hold_ids_arr = (
            df.groupby("hold_id", sort=True)["sample_id"]
            .count()
            .loc[lambda s: s > 1]
            .index.to_numpy(dtype=int)
        )
    else:
        val_hold_ids_arr = np.asarray(list(val_hold_ids), dtype=int)
    idx = np.arange(0, len(df) - 11, dtype=int)
    mask = np.isin(hold_id_arr[idx], val_hold_ids_arr) & np.isin(
        branch_arr[idx], np.array(["up", "down"])
    )
    idx = idx[mask]
    if len(idx) == 0:
        raise ValueError("No validation samples available for prediction metrics.")
    if len(idx) > 2500:
        idx = idx[:: max(1, len(idx) // 2500)]

    errs_1: List[Tuple[float, float]] = []
    errs_10: List[Tuple[float, float]] = []
    for sample_idx in idx:
        branch = branch_arr[sample_idx]
        A, B = model.branch_models[branch]
        x_pred = A @ state_matrix[sample_idx] + B[:, 0] * u_arr[sample_idx]
        errs_1.append(
            (
                float(x_pred[0] - state_matrix[sample_idx + 1, 0]),
                float(x_pred[1] - state_matrix[sample_idx + 1, 1]),
            )
        )
        for step in range(1, 10):
            future_idx = int(sample_idx + step)
            A, B = model.branch_models[branch_arr[future_idx]]
            x_pred = A @ x_pred + B[:, 0] * u_arr[future_idx]
        errs_10.append(
            (
                float(x_pred[0] - state_matrix[sample_idx + 10, 0]),
                float(x_pred[1] - state_matrix[sample_idx + 10, 1]),
            )
        )

    e1 = np.asarray(errs_1, dtype=float)
    e10 = np.asarray(errs_10, dtype=float)
    return {
        "one_step_rmse_s_mm": float(np.sqrt(np.mean(np.square(e1[:, 0])))),
        "one_step_rmse_p_bar": float(np.sqrt(np.mean(np.square(e1[:, 1])))),
        "ten_step_rmse_s_mm": float(np.sqrt(np.mean(np.square(e10[:, 0])))),
        "ten_step_rmse_p_bar": float(np.sqrt(np.mean(np.square(e10[:, 1])))),
    }


def build_lqi_matrices(A: np.ndarray, B: np.ndarray, ts_s: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    C_track = np.zeros((1, A.shape[0]), dtype=float)
    C_track[0, 0] = 1.0
    A_aug = np.block(
        [
            [A, np.zeros((A.shape[0], 1), dtype=float)],
            [ts_s * C_track, np.ones((1, 1), dtype=float)],
        ]
    )
    B_aug = np.vstack([B, np.zeros((1, B.shape[1]), dtype=float)])
    return A_aug, B_aug, C_track


def solve_discrete_lqr(A: np.ndarray, B: np.ndarray, Q: np.ndarray, R: np.ndarray) -> np.ndarray:
    P = solve_discrete_are(A, B, Q, R)
    gain_rhs = B.T @ P @ A
    gain_lhs = B.T @ P @ B + R
    return solve(gain_lhs, gain_rhs)


def solve_lqi_for_model(
    model: IdentifiedModel,
    ts_s: float,
    q_diag: np.ndarray,
    r_value: float,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
]:
    Q = np.diag(q_diag)
    R = np.array([[r_value]], dtype=float)
    branch_K_aug: Dict[str, np.ndarray] = {}
    branch_A_aug: Dict[str, np.ndarray] = {}
    branch_B_aug: Dict[str, np.ndarray] = {}
    C_track: np.ndarray | None = None

    for branch, (A_branch, B_branch) in model.branch_models.items():
        A_aug, B_aug, C_track = build_lqi_matrices(A_branch, B_branch, ts_s)
        K_aug = solve_discrete_lqr(A_aug, B_aug, Q, R)
        branch_K_aug[branch] = K_aug
        branch_A_aug[branch] = A_aug
        branch_B_aug[branch] = B_aug

    if C_track is None:
        raise ValueError("Missing augmented tracking matrix.")
    return (
        0.5 * (branch_K_aug["up"] + branch_K_aug["down"]),
        0.5 * (branch_A_aug["up"] + branch_A_aug["down"]),
        0.5 * (branch_B_aug["up"] + branch_B_aug["down"]),
        C_track,
        branch_K_aug,
        branch_A_aug,
        branch_B_aug,
    )


def _gain_stats(branch_K_aug: Dict[str, np.ndarray]) -> Tuple[float, float, float, float]:
    k_up = np.asarray(branch_K_aug["up"]).reshape(-1)
    k_down = np.asarray(branch_K_aug["down"]).reshape(-1)
    max_abs_gain_s = float(max(abs(k_up[0]), abs(k_down[0])))
    max_abs_gain_p = float(max(abs(k_up[1]), abs(k_down[1])))
    gain_l1_norm = float(np.sum(np.abs(k_up)) + np.sum(np.abs(k_down)))
    min_abs_integral_gain = float(min(abs(k_up[-1]), abs(k_down[-1])))
    return gain_l1_norm, max_abs_gain_s, max_abs_gain_p, min_abs_integral_gain


def _alpha_from_tau(ts_s: float, tau_s: float) -> float:
    if tau_s <= 0.0:
        return 0.0
    return float(np.exp(-ts_s / tau_s))


def alpha_from_tau(ts_s: float, tau_s: float) -> float:
    return _alpha_from_tau(ts_s, tau_s)


def _lookup_value(branch: str, s_ref: float, lookup_maps: Dict[str, LookupMap]) -> Tuple[float, float]:
    u_ff = float(lookup_maps[f"{branch}_inverse_u"].evaluate(s_ref))
    p_ff = float(lookup_maps[f"{branch}_inverse_p"].evaluate(s_ref))
    return u_ff, p_ff


def lookup_value(branch: str, s_ref: float, lookup_maps: Dict[str, LookupMap]) -> Tuple[float, float]:
    return _lookup_value(branch, s_ref, lookup_maps)


def build_reference_scenarios(lookup_tables: Dict[str, pd.DataFrame], steps_per_level: int) -> Dict[str, np.ndarray]:
    all_s = np.concatenate(
        [
            lookup_tables["up"]["s_ref_mm"].to_numpy(dtype=float),
            lookup_tables["down"]["s_ref_mm"].to_numpy(dtype=float),
        ]
    )
    lo = float(np.quantile(all_s, 0.20))
    mid = float(np.quantile(all_s, 0.50))
    hi = float(np.quantile(all_s, 0.80))
    q35 = float(np.quantile(all_s, 0.35))
    q65 = float(np.quantile(all_s, 0.65))

    def staircase(levels: List[float]) -> np.ndarray:
        return np.concatenate([np.full(steps_per_level, level, dtype=float) for level in levels])

    return {
        "positive": staircase([lo, mid, hi]),
        "negative": staircase([hi, mid, lo]),
        "mixed": staircase([q35, hi, mid, lo, q65]),
        "interior": staircase([q35, q65, mid, q65, mid]),
    }


def _safe_nanmean(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.all(~np.isfinite(arr)):
        return float("nan")
    return float(np.nanmean(arr))


def _safe_mean(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return float("nan")
    return float(np.mean(arr))


def _segment_bounds(ref_mm: np.ndarray) -> List[Tuple[int, int]]:
    if len(ref_mm) == 0:
        return []
    change_idx = np.flatnonzero(np.abs(np.diff(ref_mm)) > 1e-12) + 1
    starts = np.r_[0, change_idx]
    ends = np.r_[change_idx, len(ref_mm)]
    return [(int(start), int(end)) for start, end in zip(starts, ends)]


def segment_bounds(ref_mm: np.ndarray) -> List[Tuple[int, int]]:
    return _segment_bounds(ref_mm)


def _first_crossing_index(values: np.ndarray, threshold: float) -> int | None:
    idx = np.flatnonzero(values >= threshold)
    return None if len(idx) == 0 else int(idx[0])


def _compute_segment_metrics(
    y0_mm: float,
    ref_level_mm: float,
    y_seg_mm: np.ndarray,
    err_seg_mm: np.ndarray,
    u_cmd_seg_v: np.ndarray,
    u_ff_seg_v: np.ndarray,
    v_seg_mm_s: np.ndarray,
    sat_seg: np.ndarray,
    ts_s: float,
    config: ControllerDesignConfig,
) -> Dict[str, float]:
    tail_n = max(5, int(np.ceil(config.settled_window_frac * len(y_seg_mm))))
    tail_err = err_seg_mm[-tail_n:]
    tail_u = u_cmd_seg_v[-tail_n:]
    tail_v = v_seg_mm_s[-tail_n:]
    control_residual = u_cmd_seg_v - u_ff_seg_v
    local_t = np.arange(len(err_seg_mm), dtype=float) * ts_s

    metrics = {
        "settled_error_rms_mm": float(np.sqrt(np.mean(np.square(tail_err)))),
        "settled_error_p2p_mm": float(np.max(tail_err) - np.min(tail_err)),
        "settled_command_p2p_v": float(np.max(tail_u) - np.min(tail_u)),
        "settled_velocity_rms_mm_s": float(np.sqrt(np.mean(np.square(tail_v)))),
        "steady_state_bias_mm": float(np.mean(tail_err)),
        "IAE_mm_s": float(np.sum(np.abs(err_seg_mm)) * ts_s),
        "ISE_mm2_s": float(np.sum(np.square(err_seg_mm)) * ts_s),
        "ITAE_mm_s2": float(np.sum(local_t * np.abs(err_seg_mm)) * ts_s),
        "control_residual_rms_v": float(np.sqrt(np.mean(np.square(control_residual)))),
        "command_total_variation_v": float(np.sum(np.abs(np.diff(u_cmd_seg_v)))) if len(u_cmd_seg_v) > 1 else 0.0,
        "saturation_fraction": float(np.mean(sat_seg)),
    }

    step_mm = float(ref_level_mm - y0_mm)
    if abs(step_mm) < config.step_metric_min_step_mm:
        metrics.update(
            {
                "rise_time_s_10_90": float("nan"),
                "settling_time_s_2pct": float("nan"),
                "peak_time_s": float("nan"),
                "overshoot_mm": float("nan"),
                "overshoot_pct": float("nan"),
            }
        )
        return metrics

    sign = 1.0 if step_mm >= 0 else -1.0
    progress = sign * (y_seg_mm - y0_mm) / abs(step_mm)
    rise_lo, rise_hi = config.rise_time_bounds
    rise_lo_idx = _first_crossing_index(progress, rise_lo)
    rise_hi_idx = _first_crossing_index(progress, rise_hi)
    if rise_lo_idx is None or rise_hi_idx is None or rise_hi_idx < rise_lo_idx:
        rise_time = float("nan")
    else:
        rise_time = float((rise_hi_idx - rise_lo_idx) * ts_s)

    band_mm = max(config.settling_band_frac * abs(step_mm), config.settling_band_floor_mm)
    deviation = np.abs(y_seg_mm - ref_level_mm)
    settling_time = float("nan")
    for idx in range(len(deviation)):
        if np.all(deviation[idx:] <= band_mm):
            settling_time = float(idx * ts_s)
            break

    directional_response = sign * (y_seg_mm - y0_mm)
    peak_idx = int(np.argmax(directional_response))
    peak_time = float(peak_idx * ts_s)
    overshoot_mm = float(max(np.max(sign * (y_seg_mm - ref_level_mm)), 0.0))
    overshoot_pct = float(100.0 * overshoot_mm / abs(step_mm))
    metrics.update(
        {
            "rise_time_s_10_90": rise_time,
            "settling_time_s_2pct": settling_time,
            "peak_time_s": peak_time,
            "overshoot_mm": overshoot_mm,
            "overshoot_pct": overshoot_pct,
        }
    )
    return metrics


def compute_segment_metrics(
    y0_mm: float,
    ref_level_mm: float,
    y_seg_mm: np.ndarray,
    err_seg_mm: np.ndarray,
    u_cmd_seg_v: np.ndarray,
    u_ff_seg_v: np.ndarray,
    v_seg_mm_s: np.ndarray,
    sat_seg: np.ndarray,
    ts_s: float,
    config: ControllerDesignConfig,
) -> Dict[str, float]:
    return _compute_segment_metrics(
        y0_mm=y0_mm,
        ref_level_mm=ref_level_mm,
        y_seg_mm=y_seg_mm,
        err_seg_mm=err_seg_mm,
        u_cmd_seg_v=u_cmd_seg_v,
        u_ff_seg_v=u_ff_seg_v,
        v_seg_mm_s=v_seg_mm_s,
        sat_seg=sat_seg,
        ts_s=ts_s,
        config=config,
    )


def _nan_controller_metrics() -> Dict[str, float]:
    return {
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


def _simulate_controller(
    model: IdentifiedModel,
    lookup_maps: Dict[str, LookupMap],
    ref_mm: np.ndarray,
    ts_s: float,
    control_limits: Dict[str, float | str | bool],
    steady_df: pd.DataFrame,
    config: ControllerDesignConfig,
    branch_K_aug: Dict[str, np.ndarray],
    transient_profiles: TransientProfiles | None = None,
    use_creep_compensation: bool = False,
    use_dynamic_pressure_reference: bool = False,
    use_pressure_filter: bool = False,
    pressure_filter_tau_s: float | None = None,
) -> Dict[str, float]:
    x_bar_2, u_bar = build_global_means(steady_df)
    x = np.zeros(model.A.shape[0], dtype=float)
    e_int = 0.0
    s_meas_hist: List[float] = [float(x_bar_2[0])]
    s_vel_lp_hist: List[float] = [float(x_bar_2[0])]
    p_hat = float(x_bar_2[1])
    v_hat_fb = 0.0
    s_vel_lp_state = float(x_bar_2[0])
    active_branch = "up"
    branch_deadband = float(control_limits["branch_deadband_mm"])
    position_deadband = float(control_limits["position_deadband_mm"])
    s_lower_bound, s_upper_bound = _practical_displacement_bounds(lookup_maps)
    alpha_v = (
        _alpha_from_tau(ts_s, model.velocity_filter_tau_s)
        if model.uses_velocity_state and model.velocity_filter_tau_s is not None
        else None
    )
    alpha_s_vel = (
        _alpha_from_tau(ts_s, model.velocity_position_prefilter_tau_s)
        if model.uses_velocity_state and model.velocity_position_prefilter_tau_s is not None
        else None
    )
    alpha_p = (
        _alpha_from_tau(ts_s, pressure_filter_tau_s)
        if use_pressure_filter and pressure_filter_tau_s is not None
        else None
    )

    s_hist_mm: List[float] = []
    err_hist_mm: List[float] = []
    u_cmd_hist_v: List[float] = []
    u_ff_hist_v: List[float] = []
    v_slope_raw_hist_mm_s: List[float] = []
    sat_hist: List[bool] = []

    initial_output_mm = float(x_bar_2[0])
    previous_ref_mm: float | None = None
    step_delta_s_mm = 0.0
    step_branch: str | None = None
    time_since_step_s = 0.0
    p_step_bar = float(x_bar_2[1])

    for s_ref_value in ref_mm:
        s_ref = float(s_ref_value)
        s_meas_now = float(x[0] + x_bar_2[0])
        p_meas_now = float(x[1] + x_bar_2[1])
        if alpha_p is not None:
            p_hat = alpha_p * p_hat + (1.0 - alpha_p) * p_meas_now
            p_feedback = p_hat
        else:
            p_feedback = p_meas_now

        s_tilde = s_meas_now - float(s_ref)
        s_tilde_fb = 0.0 if abs(s_tilde) < position_deadband else s_tilde

        if s_ref > s_meas_now + branch_deadband:
            active_branch = "up"
        elif s_ref < s_meas_now - branch_deadband:
            active_branch = "down"

        if previous_ref_mm is None:
            p_step_bar = p_meas_now
            time_since_step_s = 0.0
            step_delta_s_mm = 0.0
            step_branch = None
        elif abs(s_ref - previous_ref_mm) > 1e-12:
            step_delta_s_mm = s_ref - previous_ref_mm
            step_branch = "up" if step_delta_s_mm > 0.0 else "down"
            time_since_step_s = 0.0
            p_step_bar = p_meas_now

        u_ff_static, p_ff_static = _lookup_value(active_branch, s_ref, lookup_maps)
        direction_matches = step_branch == active_branch and abs(step_delta_s_mm) > 1e-12
        u_ff = u_ff_static
        if use_creep_compensation:
            if transient_profiles is None or active_branch not in transient_profiles.creep_by_branch:
                raise ValueError("Creep compensation is enabled but no branch profile is loaded.")
            creep_profile = transient_profiles.creep_by_branch[active_branch]
            if direction_matches and time_since_step_s <= creep_profile.horizon_s + 1e-12:
                creep_fraction = creep_profile.evaluate(time_since_step_s)
                s_virtual = np.clip(
                    s_ref + step_delta_s_mm * creep_fraction,
                    lookup_maps[f"{active_branch}_inverse_u"].grid_x[0],
                    lookup_maps[f"{active_branch}_inverse_u"].grid_x[-1],
                )
                u_ff = float(lookup_maps[f"{active_branch}_inverse_u"].evaluate(s_virtual))

        p_ff = p_ff_static
        if use_dynamic_pressure_reference:
            if transient_profiles is None or active_branch not in transient_profiles.pressure_by_branch:
                raise ValueError("Dynamic pressure reference is enabled but no branch profile is loaded.")
            pressure_profile = transient_profiles.pressure_by_branch[active_branch]
            if direction_matches and time_since_step_s <= pressure_profile.horizon_s + 1e-12:
                pressure_progress = pressure_profile.evaluate(time_since_step_s)
                p_ff = p_step_bar + (p_ff_static - p_step_bar) * pressure_progress
        p_tilde = p_feedback - p_ff
        if alpha_s_vel is not None:
            s_vel_lp_state = update_first_order_lowpass(s_vel_lp_state, s_meas_now, alpha_s_vel)
        else:
            s_vel_lp_state = s_meas_now
        v_slope_raw = compute_latest_linear_fit_slope(
            s_vel_lp_hist + [s_vel_lp_state],
            ts_s,
            model.velocity_slope_window_samples or 2,
        )

        if model.uses_velocity_state:
            if alpha_v is None:
                raise ValueError("Velocity-state controller is missing alpha_v.")
            v_hat_fb = alpha_v * v_hat_fb + (1.0 - alpha_v) * v_slope_raw
            feedback_vec = np.array([s_tilde_fb, p_tilde, v_hat_fb, e_int], dtype=float)
        else:
            feedback_vec = np.array([s_tilde_fb, p_tilde, e_int], dtype=float)

        K_eff = branch_K_aug[active_branch].reshape(-1)
        u_cmd_raw = float(u_ff - (K_eff @ feedback_vec))
        u_cmd = float(np.clip(u_cmd_raw, float(control_limits["u_min_v"]), float(control_limits["u_max_v"])))
        saturated = bool(
            u_cmd >= float(control_limits["u_max_v"]) - 1e-12
            or u_cmd <= float(control_limits["u_min_v"]) + 1e-12
        )
        u_tilde = u_cmd - u_bar

        A_eff, B_eff = model.branch_models[active_branch]
        x_next = A_eff @ x + B_eff[:, 0] * u_tilde
        if not (np.all(np.isfinite(x_next)) and np.isfinite(u_cmd) and np.isfinite(v_hat_fb)):
            return _nan_controller_metrics()
        s_meas_next = float(x_next[0] + x_bar_2[0])
        if not np.isfinite(s_meas_next):
            return _nan_controller_metrics()
        if s_meas_next < s_lower_bound or s_meas_next > s_upper_bound:
            return _nan_controller_metrics()
        error_now = s_meas_now - float(s_ref)
        if u_cmd >= float(control_limits["u_max_v"]) - 1e-12 and error_now < 0:
            pass
        elif u_cmd <= float(control_limits["u_min_v"]) + 1e-12 and error_now > 0:
            pass
        else:
            e_int = float(e_int + (s_meas_next - float(s_ref)) * ts_s)
        e_int = float(
            np.clip(
                e_int,
                -float(control_limits["integral_state_limit"]),
                float(control_limits["integral_state_limit"]),
            )
        )

        s_hist_mm.append(s_meas_next)
        err_hist_mm.append(s_meas_next - float(s_ref))
        u_cmd_hist_v.append(u_cmd)
        u_ff_hist_v.append(u_ff)
        v_slope_raw_hist_mm_s.append(v_slope_raw)
        sat_hist.append(saturated)

        s_meas_hist.append(s_meas_now)
        s_vel_lp_hist.append(s_vel_lp_state)
        x = x_next
        previous_ref_mm = s_ref
        time_since_step_s += ts_s

    s_hist_arr = np.asarray(s_hist_mm, dtype=float)
    err_hist_arr = np.asarray(err_hist_mm, dtype=float)
    u_cmd_hist_arr = np.asarray(u_cmd_hist_v, dtype=float)
    u_ff_hist_arr = np.asarray(u_ff_hist_v, dtype=float)
    v_raw_hist_arr = np.asarray(v_slope_raw_hist_mm_s, dtype=float)
    sat_hist_arr = np.asarray(sat_hist, dtype=bool)

    segment_metrics: List[Dict[str, float]] = []
    bounds = _segment_bounds(ref_mm)
    for start, end in bounds:
        y0_mm = initial_output_mm if start == 0 else float(s_hist_arr[start - 1])
        segment_metrics.append(
            _compute_segment_metrics(
                y0_mm=y0_mm,
                ref_level_mm=float(ref_mm[start]),
                y_seg_mm=s_hist_arr[start:end],
                err_seg_mm=err_hist_arr[start:end],
                u_cmd_seg_v=u_cmd_hist_arr[start:end],
                u_ff_seg_v=u_ff_hist_arr[start:end],
                v_seg_mm_s=v_raw_hist_arr[start:end],
                sat_seg=sat_hist_arr[start:end],
                ts_s=ts_s,
                config=config,
            )
        )

    scenario_itae = float(np.sum([m["ITAE_mm_s2"] for m in segment_metrics]))
    return {
        "tracking_rms_mm": float(np.sqrt(np.mean(np.square(err_hist_arr)))),
        "tracking_overshoot_mm": float(np.max(np.abs(err_hist_arr))),
        "control_effort_mean_abs_v": float(np.mean(np.abs(u_cmd_hist_arr - u_ff_hist_arr))),
        "settled_error_rms_mm": _safe_nanmean([m["settled_error_rms_mm"] for m in segment_metrics]),
        "settled_error_p2p_mm": _safe_nanmean([m["settled_error_p2p_mm"] for m in segment_metrics]),
        "settled_command_p2p_v": _safe_nanmean([m["settled_command_p2p_v"] for m in segment_metrics]),
        "settled_velocity_rms_mm_s": _safe_nanmean([m["settled_velocity_rms_mm_s"] for m in segment_metrics]),
        "rise_time_s_10_90": _safe_nanmean([m["rise_time_s_10_90"] for m in segment_metrics]),
        "settling_time_s_2pct": _safe_nanmean([m["settling_time_s_2pct"] for m in segment_metrics]),
        "peak_time_s": _safe_nanmean([m["peak_time_s"] for m in segment_metrics]),
        "overshoot_mm": _safe_nanmean([m["overshoot_mm"] for m in segment_metrics]),
        "overshoot_pct": _safe_nanmean([m["overshoot_pct"] for m in segment_metrics]),
        "steady_state_bias_mm": _safe_nanmean([m["steady_state_bias_mm"] for m in segment_metrics]),
        "IAE_mm_s": float(np.sum(np.abs(err_hist_arr)) * ts_s),
        "ISE_mm2_s": float(np.sum(np.square(err_hist_arr)) * ts_s),
        "ITAE_mm_s2": scenario_itae,
        "control_residual_rms_v": float(np.sqrt(np.mean(np.square(u_cmd_hist_arr - u_ff_hist_arr)))),
        "command_total_variation_v": float(np.sum(np.abs(np.diff(u_cmd_hist_arr)))) if len(u_cmd_hist_arr) > 1 else 0.0,
        "saturation_fraction": float(np.mean(sat_hist_arr)),
    }


def evaluate_pressure_filter_sensitivity(
    result: ControllerResult,
    lookup_maps: Dict[str, LookupMap],
    ts_s: float,
    control_limits: Dict[str, float | str | bool],
    steady_df: pd.DataFrame,
    config: ControllerDesignConfig,
    transient_profiles: TransientProfiles,
) -> Tuple[bool, float | None]:
    lookup_tables = build_lookup_tables(lookup_maps)
    hold_steps = max(60, int(round(1.1 / ts_s)))
    scenarios = build_reference_scenarios(lookup_tables, hold_steps)

    base_metrics = []
    pressure_control_limits = dict(control_limits)
    pressure_control_limits["branch_deadband_mm"] = float(result.branch_deadband_mm)
    pressure_control_limits["position_deadband_mm"] = float(result.position_deadband_mm)
    for ref_values in scenarios.values():
        base_metrics.append(
            _simulate_controller(
                model=result.evaluation_model,
                lookup_maps=lookup_maps,
                ref_mm=ref_values,
                ts_s=ts_s,
                control_limits=pressure_control_limits,
                steady_df=steady_df,
                config=config,
                branch_K_aug=result.branch_K_aug,
                transient_profiles=transient_profiles,
                use_creep_compensation=True,
                use_dynamic_pressure_reference=True,
                use_pressure_filter=False,
            )
        )
    base_tracking = float(np.mean([m["tracking_rms_mm"] for m in base_metrics]))
    base_settled = float(np.mean([m["settled_error_p2p_mm"] for m in base_metrics]))

    for tau_p_s in (0.02, 0.05):
        metrics = []
        for ref_values in scenarios.values():
            metrics.append(
                _simulate_controller(
                    model=result.evaluation_model,
                    lookup_maps=lookup_maps,
                    ref_mm=ref_values,
                    ts_s=ts_s,
                    control_limits=pressure_control_limits,
                    steady_df=steady_df,
                    config=config,
                    branch_K_aug=result.branch_K_aug,
                    transient_profiles=transient_profiles,
                    use_creep_compensation=True,
                    use_dynamic_pressure_reference=True,
                    use_pressure_filter=True,
                    pressure_filter_tau_s=tau_p_s,
                )
            )
        tracking = float(np.mean([m["tracking_rms_mm"] for m in metrics]))
        settled = float(np.mean([m["settled_error_p2p_mm"] for m in metrics]))
        if settled <= 0.90 * base_settled and tracking <= 1.05 * base_tracking:
            return True, float(tau_p_s)
    return False, None


def _boundary_fields_for_result(result: ControllerResult, config: ControllerDesignConfig) -> str:
    fields: List[str] = []
    if result.q_s in (min(config.q_s_grid), max(config.q_s_grid)):
        fields.append("q_s")
    if result.q_p in (min(config.q_p_grid), max(config.q_p_grid)):
        fields.append("q_p")
    if result.q_i in (min(config.q_i_grid), max(config.q_i_grid)):
        fields.append("q_i")
    if result.r_value in (min(config.r_grid), max(config.r_grid)):
        fields.append("r")
    if result.branch_deadband_mm in (min(config.branch_deadband_grid_mm), max(config.branch_deadband_grid_mm)):
        fields.append("branch_deadband_mm")
    if result.position_deadband_mm in (
        min(config.position_deadband_grid_mm),
        max(config.position_deadband_grid_mm),
    ):
        fields.append("position_deadband_mm")
    if result.uses_velocity_state:
        if result.q_v is not None and result.q_v in (min(config.q_v_grid), max(config.q_v_grid)):
            fields.append("q_v")
        if (
            result.velocity_filter_tau_s is not None
            and result.velocity_filter_tau_s in (min(config.vel_tau_grid_s), max(config.vel_tau_grid_s))
        ):
            fields.append("velocity_filter_tau_s")
        if (
            result.velocity_slope_window_samples is not None
            and result.velocity_slope_window_samples
            in (min(config.vel_slope_window_grid_samples), max(config.vel_slope_window_grid_samples))
        ):
            fields.append("velocity_slope_window_samples")
        if (
            result.velocity_position_prefilter_tau_s is not None
            and result.velocity_position_prefilter_tau_s
            in (
                min(config.vel_position_prefilter_tau_grid_s),
                max(config.vel_position_prefilter_tau_grid_s),
            )
        ):
            fields.append("velocity_position_prefilter_tau_s")
    return "|".join(fields)


def _score_value(value: float) -> float:
    if value is None or not np.isfinite(value):
        return float("inf")
    return float(value)


def _compute_selection_best_tracking(results: Sequence[ControllerResult]) -> float:
    tracking_values = [_score_value(result.tracking_rms_mm) for result in results]
    return min(tracking_values) if tracking_values else float("inf")


def _selection_tracking_gate_passed(result: ControllerResult, best_tracking_rms_mm: float) -> bool:
    if not np.isfinite(best_tracking_rms_mm):
        return False
    tracking_rms_mm = _score_value(result.tracking_rms_mm)
    return tracking_rms_mm <= SELECTION_TRACKING_GATE_FACTOR * best_tracking_rms_mm


def _selection_eligibility_passed(result: ControllerResult, best_tracking_rms_mm: float) -> bool:
    if not result.validation_ok:
        return False
    if not _selection_tracking_gate_passed(result, best_tracking_rms_mm):
        return False
    if _score_value(result.min_abs_integral_gain) < SELECTION_MIN_ABS_INTEGRAL_GAIN:
        return False
    return all(np.isfinite(getattr(result, attr)) for attr in SELECTION_PRIMARY_METRIC_ATTRS)


def _selection_metric_mins(results: Sequence[ControllerResult]) -> Dict[str, float]:
    mins: Dict[str, float] = {}
    for attr in SELECTION_PRIMARY_METRIC_ATTRS:
        values = [float(getattr(result, attr)) for result in results if np.isfinite(getattr(result, attr))]
        mins[attr] = min(values) if values else float("inf")
    return mins


def _selection_score_primary(result: ControllerResult, metric_mins: Dict[str, float]) -> float:
    score = 0.0
    for attr in SELECTION_PRIMARY_METRIC_ATTRS:
        metric_value = float(getattr(result, attr))
        metric_min = float(metric_mins[attr])
        if not np.isfinite(metric_value) or not np.isfinite(metric_min):
            return float("inf")
        metric_floor = float(SELECTION_NORMALIZATION_FLOORS.get(attr, 0.0))
        denom = max(abs(metric_min), metric_floor, 1e-9)
        normalized = metric_value / denom
        if attr == "max_abs_gain_s":
            normalized = normalized ** SELECTION_GAIN_PENALTY_EXPONENT
        score += SELECTION_PRIMARY_WEIGHTS[attr] * normalized
    return float(score)


def _selection_tie_break_key(result: ControllerResult) -> Tuple[float, float, float]:
    return tuple(_score_value(getattr(result, attr)) for attr in SELECTION_TIE_BREAK_ATTRS)


def _rank_controller_candidates(
    results: Sequence[ControllerResult],
) -> Tuple[List[ControllerResult], bool, float]:
    result_list = list(results)
    if not result_list:
        return [], False, float("inf")

    best_tracking_rms_mm = _compute_selection_best_tracking(result_list)
    eligible_results = [result for result in result_list if _selection_eligibility_passed(result, best_tracking_rms_mm)]
    ranked_pool = eligible_results if eligible_results else result_list
    metric_mins = _selection_metric_mins(ranked_pool)

    for result in result_list:
        result.selection_tracking_gate_passed = _selection_tracking_gate_passed(result, best_tracking_rms_mm)
        result.selection_eligibility_passed = _selection_eligibility_passed(result, best_tracking_rms_mm)
        result.selection_score_primary = _selection_score_primary(result, metric_mins)

    ranked_results = sorted(
        ranked_pool,
        key=lambda result: (
            _score_value(result.selection_score_primary),
            *_selection_tie_break_key(result),
        ),
    )
    return ranked_results, bool(eligible_results), best_tracking_rms_mm


def rank_results_for_selection(
    results: Sequence[ControllerResult],
) -> Tuple[List[ControllerResult], bool, float]:
    return _rank_controller_candidates(results)


def _feedback_state_names(result: ControllerResult) -> List[str]:
    return list(result.model.runtime_observable_names) + ["e_I"]


def _make_result(
    config_id: str,
    model: IdentifiedModel,
    evaluation_model: IdentifiedModel,
    q_diag: np.ndarray,
    q_s: float,
    q_p: float,
    q_i: float,
    q_v: float | None,
    r_value: float,
    branch_deadband_mm: float,
    position_deadband_mm: float,
    prediction_metrics: Dict[str, float],
    scenario_metrics: Dict[str, Dict[str, float]],
    K_aug: np.ndarray,
    A_aug: np.ndarray,
    B_aug: np.ndarray,
    C_track: np.ndarray,
    branch_K_aug: Dict[str, np.ndarray],
    branch_A_aug: Dict[str, np.ndarray],
    branch_B_aug: Dict[str, np.ndarray],
    controller_ts_s: float,
    integral_state_limit: float,
    config: ControllerDesignConfig,
) -> ControllerResult:
    all_metrics = scenario_metrics["all"]
    gain_l1_norm, max_abs_gain_s, max_abs_gain_p, min_abs_integral_gain = _gain_stats(branch_K_aug)
    result = ControllerResult(
        config_id=config_id,
        variant=model.variant,
        model=model,
        evaluation_model=evaluation_model,
        K_aug=K_aug,
        A_aug=A_aug,
        B_aug=B_aug,
        C_track=C_track,
        q_diag=q_diag.copy(),
        q_s=q_s,
        q_p=q_p,
        q_i=q_i,
        q_v=q_v,
        r_value=r_value,
        branch_deadband_mm=branch_deadband_mm,
        position_deadband_mm=position_deadband_mm,
        one_step_rmse_s_mm=prediction_metrics["one_step_rmse_s_mm"],
        one_step_rmse_p_bar=prediction_metrics["one_step_rmse_p_bar"],
        ten_step_rmse_s_mm=prediction_metrics["ten_step_rmse_s_mm"],
        ten_step_rmse_p_bar=prediction_metrics["ten_step_rmse_p_bar"],
        tracking_rms_mm=all_metrics["tracking_rms_mm"],
        tracking_overshoot_mm=all_metrics["tracking_overshoot_mm"],
        control_effort_mean_abs_v=all_metrics["control_effort_mean_abs_v"],
        settled_error_rms_mm=all_metrics["settled_error_rms_mm"],
        settled_error_p2p_mm=all_metrics["settled_error_p2p_mm"],
        settled_command_p2p_v=all_metrics["settled_command_p2p_v"],
        settled_velocity_rms_mm_s=all_metrics["settled_velocity_rms_mm_s"],
        rise_time_s_10_90=all_metrics["rise_time_s_10_90"],
        settling_time_s_2pct=all_metrics["settling_time_s_2pct"],
        peak_time_s=all_metrics["peak_time_s"],
        overshoot_mm=all_metrics["overshoot_mm"],
        overshoot_pct=all_metrics["overshoot_pct"],
        steady_state_bias_mm=all_metrics["steady_state_bias_mm"],
        abs_steady_state_bias_mm=abs(all_metrics["steady_state_bias_mm"]),
        iae_mm_s=all_metrics["IAE_mm_s"],
        ise_mm2_s=all_metrics["ISE_mm2_s"],
        itae_mm_s2=all_metrics["ITAE_mm_s2"],
        control_residual_rms_v=all_metrics["control_residual_rms_v"],
        command_total_variation_v=all_metrics["command_total_variation_v"],
        saturation_fraction=all_metrics["saturation_fraction"],
        integral_state_limit=float(integral_state_limit),
        validation_ok=model.stable_model,
        branch_K_aug=branch_K_aug,
        branch_A_aug=branch_A_aug,
        branch_B_aug=branch_B_aug,
        gain_l1_norm=gain_l1_norm,
        max_abs_gain_s=max_abs_gain_s,
        max_abs_gain_p=max_abs_gain_p,
        min_abs_integral_gain=min_abs_integral_gain,
        uses_velocity_state=model.uses_velocity_state,
        velocity_filter_tau_s=model.velocity_filter_tau_s,
        velocity_slope_window_samples=model.velocity_slope_window_samples,
        velocity_position_prefilter_tau_s=model.velocity_position_prefilter_tau_s,
        velocity_filter_alpha=(
            _alpha_from_tau(controller_ts_s, model.velocity_filter_tau_s)
            if model.uses_velocity_state and model.velocity_filter_tau_s is not None
            else None
        ),
        velocity_position_prefilter_alpha=(
            _alpha_from_tau(controller_ts_s, model.velocity_position_prefilter_tau_s)
            if model.uses_velocity_state and model.velocity_position_prefilter_tau_s is not None
            else None
        ),
        uses_creep_compensation=True,
        uses_dynamic_pressure_reference=True,
        transient_profile_horizon_s=3.99,
        transient_lookup_source="open_loop_steady_state_lookup",
        uses_pressure_filter=False,
        pressure_filter_tau_s=None,
        scenario_metrics=scenario_metrics,
        search_grid_boundary_hit=False,
        search_grid_boundary_fields="",
        evaluation_model_stable=evaluation_model.stable_model,
        evaluation_spectral_radius=evaluation_model.spectral_radius,
    )
    scenario_rows = [
        metrics for name, metrics in scenario_metrics.items() if name != "all"
    ]
    scenario_finite = all(
        np.isfinite(float(metrics.get("tracking_rms_mm", float("nan"))))
        for metrics in scenario_rows
    )
    scenario_saturation = [
        float(metrics["saturation_fraction"])
        for metrics in scenario_rows
        if np.isfinite(float(metrics.get("saturation_fraction", float("nan"))))
    ]
    saturation_lockup = bool(
        scenario_saturation
        and len(scenario_saturation) == len(scenario_rows)
        and all(value >= 0.99 for value in scenario_saturation)
    )
    aggregate_primary_finite = all(
        np.isfinite(float(getattr(result, attr)))
        for attr in SELECTION_PRIMARY_METRIC_ATTRS
    )
    result.validation_ok = bool(
        model.stable_model
        and evaluation_model.stable_model
        and scenario_finite
        and aggregate_primary_finite
        and not saturation_lockup
    )
    boundary_fields = _boundary_fields_for_result(result, config)
    result.search_grid_boundary_fields = boundary_fields
    result.search_grid_boundary_hit = boundary_fields != ""
    return result


def _init_lqi_worker(
    model: IdentifiedModel,
    evaluation_model: IdentifiedModel,
    lookup_maps: Dict[str, LookupMap],
    reference_scenarios: Dict[str, np.ndarray],
    ts_s: float,
    eval_ts_s: float,
    control_limits: Dict[str, float | str | bool],
    steady_df: pd.DataFrame,
    config: ControllerDesignConfig,
    prediction_metrics: Dict[str, float],
    transient_profiles: TransientProfiles,
) -> None:
    global _LQI_WORKER_MODEL
    global _LQI_WORKER_EVAL_MODEL
    global _LQI_WORKER_LOOKUP_MAPS
    global _LQI_WORKER_REFERENCE_SCENARIOS
    global _LQI_WORKER_TS_S
    global _LQI_WORKER_EVAL_TS_S
    global _LQI_WORKER_CONTROL_LIMITS
    global _LQI_WORKER_STEADY_DF
    global _LQI_WORKER_CONFIG
    global _LQI_WORKER_PREDICTION_METRICS
    global _LQI_WORKER_TRANSIENT_PROFILES

    _LQI_WORKER_MODEL = model
    _LQI_WORKER_EVAL_MODEL = evaluation_model
    _LQI_WORKER_LOOKUP_MAPS = lookup_maps
    _LQI_WORKER_REFERENCE_SCENARIOS = reference_scenarios
    _LQI_WORKER_TS_S = ts_s
    _LQI_WORKER_EVAL_TS_S = eval_ts_s
    _LQI_WORKER_CONTROL_LIMITS = control_limits
    _LQI_WORKER_STEADY_DF = steady_df
    _LQI_WORKER_CONFIG = config
    _LQI_WORKER_PREDICTION_METRICS = prediction_metrics
    _LQI_WORKER_TRANSIENT_PROFILES = transient_profiles


def _evaluate_lqi_task(
    task: Tuple[str, float, float, float, float | None, float, float, float]
) -> ControllerResult:
    if (
        _LQI_WORKER_MODEL is None
        or _LQI_WORKER_EVAL_MODEL is None
        or _LQI_WORKER_LOOKUP_MAPS is None
        or _LQI_WORKER_REFERENCE_SCENARIOS is None
        or _LQI_WORKER_TS_S is None
        or _LQI_WORKER_EVAL_TS_S is None
        or _LQI_WORKER_CONTROL_LIMITS is None
        or _LQI_WORKER_STEADY_DF is None
        or _LQI_WORKER_CONFIG is None
        or _LQI_WORKER_PREDICTION_METRICS is None
        or _LQI_WORKER_TRANSIENT_PROFILES is None
    ):
        raise RuntimeError("LQI worker context is not initialized.")

    config_id, q_s, q_p, q_i, q_v, r_value, branch_deadband_mm, position_deadband_mm = task
    if _LQI_WORKER_MODEL.uses_velocity_state:
        if q_v is None:
            raise ValueError("Velocity-state LQI task is missing q_v.")
        q_diag = np.array([q_s, q_p, q_v, q_i], dtype=float)
    else:
        q_diag = np.array([q_s, q_p, q_i], dtype=float)

    (
        K_aug,
        A_aug,
        B_aug,
        C_track,
        branch_K_aug,
        branch_A_aug,
        branch_B_aug,
    ) = solve_lqi_for_model(_LQI_WORKER_MODEL, _LQI_WORKER_TS_S, q_diag, r_value)

    eval_control_limits = dict(_LQI_WORKER_CONTROL_LIMITS)
    eval_control_limits["branch_deadband_mm"] = float(branch_deadband_mm)
    eval_control_limits["position_deadband_mm"] = float(position_deadband_mm)
    scenario_metrics = {
        name: _simulate_controller(
            model=_LQI_WORKER_EVAL_MODEL,
            lookup_maps=_LQI_WORKER_LOOKUP_MAPS,
            ref_mm=ref_values,
            ts_s=_LQI_WORKER_EVAL_TS_S,
            control_limits=eval_control_limits,
            steady_df=_LQI_WORKER_STEADY_DF,
            config=_LQI_WORKER_CONFIG,
            branch_K_aug=branch_K_aug,
            transient_profiles=_LQI_WORKER_TRANSIENT_PROFILES,
            use_creep_compensation=True,
            use_dynamic_pressure_reference=_LQI_WORKER_MODEL.controller_family != "branchwise_pid_baseline",
        )
        for name, ref_values in _LQI_WORKER_REFERENCE_SCENARIOS.items()
    }
    scenario_metrics["all"] = _aggregate_scenario_metrics(scenario_metrics)
    return _make_result(
        config_id=config_id,
        model=_LQI_WORKER_MODEL,
        evaluation_model=_LQI_WORKER_EVAL_MODEL,
        q_diag=q_diag,
        q_s=float(q_s),
        q_p=float(q_p),
        q_i=float(q_i),
        q_v=None if q_v is None else float(q_v),
        r_value=float(r_value),
        branch_deadband_mm=float(branch_deadband_mm),
        position_deadband_mm=float(position_deadband_mm),
        prediction_metrics=_LQI_WORKER_PREDICTION_METRICS,
        scenario_metrics=scenario_metrics,
        K_aug=K_aug,
        A_aug=A_aug,
        B_aug=B_aug,
        C_track=C_track,
        branch_K_aug=branch_K_aug,
        branch_A_aug=branch_A_aug,
        branch_B_aug=branch_B_aug,
        controller_ts_s=_LQI_WORKER_TS_S,
        integral_state_limit=float(_LQI_WORKER_CONTROL_LIMITS["integral_state_limit"]),
        config=_LQI_WORKER_CONFIG,
    )


def _should_print_progress(completed: int, total: int, progress_step: int) -> bool:
    early_marks = {1, 2, 5, 10, 25, 50}
    return completed in early_marks or completed % progress_step == 0 or completed == total


def _collect_lqi_results_with_progress(
    futures: List,
    label: str,
    total: int,
    progress_step: int,
    initial_completed: int = 0,
    heartbeat_s: float = 10.0,
) -> List[ControllerResult]:
    results: List[ControllerResult] = []
    pending = set(futures)
    completed = initial_completed
    print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
    while pending:
        done, pending = wait(pending, timeout=heartbeat_s, return_when=FIRST_COMPLETED)
        if not done:
            print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
            continue
        for future in done:
            result = future.result()
            results.append(result)
            completed += 1
            if _should_print_progress(completed, total, progress_step):
                print(f"{label}: {completed}/{total} finished, {total - completed} left", flush=True)
    return results


def evaluate_controllers(
    models: Dict[str, IdentifiedModel],
    design_df: pd.DataFrame,
    design_hold_ids: Iterable[int],
    evaluation_df: pd.DataFrame,
    evaluation_hold_ids: Iterable[int],
    lookup_maps: Dict[str, LookupMap],
    design_ts_s: float,
    evaluation_ts_s: float,
    control_limits: Dict[str, float | str | bool],
    lookup_steady_df: pd.DataFrame,
    config: ControllerDesignConfig,
    transient_profiles: TransientProfiles,
) -> Tuple[Dict[str, ControllerResult], List[ControllerResult]]:
    lookup_tables = build_lookup_tables(lookup_maps)
    hold_steps = max(60, int(round(1.1 / evaluation_ts_s)))
    reference_scenarios = build_reference_scenarios(lookup_tables, hold_steps)
    mp_ctx = mp.get_context("fork")
    max_workers = max(1, (os.cpu_count() or 2) - 1)

    all_results: List[ControllerResult] = []
    family_best: Dict[str, ControllerResult] = {}
    base_total = (
        len(config.q_s_grid)
        * len(config.q_p_grid)
        * len(config.q_i_grid)
        * len(config.r_grid)
        * len(config.branch_deadband_grid_mm)
        * len(config.position_deadband_grid_mm)
    )
    base_progress_step = max(1, min(100, base_total // 10 if base_total > 0 else 1))

    base_model = models["lookup_plus_branchwise_minimal_LQI"]
    base_prediction = evaluate_prediction_metrics(
        base_model,
        evaluation_df,
        evaluation_hold_ids,
        lookup_steady_df,
        evaluation_ts_s,
    )
    if design_df is evaluation_df:
        base_eval_model = base_model
    else:
        base_eval_model = identify_candidate_models(
            evaluation_df,
            lookup_steady_df,
            evaluation_hold_ids,
            config,
        )["lookup_plus_branchwise_minimal_LQI"]
    base_results: List[ControllerResult] = []
    base_tasks = [
        (
            f"base_qs{q_s:g}_qp{q_p:g}_qi{q_i:g}_r{r_value:g}_bd{branch_deadband_mm:g}_pd{position_deadband_mm:g}",
            float(q_s),
            float(q_p),
            float(q_i),
            None,
            float(r_value),
            float(branch_deadband_mm),
            float(position_deadband_mm),
        )
        for q_s in config.q_s_grid
        for q_p in config.q_p_grid
        for q_i in config.q_i_grid
        for r_value in config.r_grid
        for branch_deadband_mm in config.branch_deadband_grid_mm
        for position_deadband_mm in config.position_deadband_grid_mm
    ]
    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=mp_ctx,
        initializer=_init_lqi_worker,
        initargs=(
            base_model,
            base_eval_model,
            lookup_maps,
            reference_scenarios,
            design_ts_s,
            evaluation_ts_s,
            control_limits,
            lookup_steady_df,
            config,
            base_prediction,
            transient_profiles,
        ),
    ) as executor:
        futures = [executor.submit(_evaluate_lqi_task, task) for task in base_tasks]
        for result in _collect_lqi_results_with_progress(
            futures=futures,
            label="LQI base sweep progress",
            total=base_total,
            progress_step=base_progress_step,
        ):
            base_results.append(result)
            all_results.append(result)

    if not base_results:
        raise ValueError("Could not synthesize any base minimal-LQI controllers.")
    ranked_base, _, _ = _rank_controller_candidates(base_results)
    best_base = ranked_base[0]
    family_best[best_base.variant] = best_base

    velocity_results: List[ControllerResult] = []
    top_base_results = ranked_base[: config.top_base_configs_for_vel]
    velocity_total = (
        len(config.vel_tau_grid_s)
        * len(config.vel_slope_window_grid_samples)
        * len(config.vel_position_prefilter_tau_grid_s)
        * len(top_base_results)
        * len(config.q_v_grid)
    )
    velocity_progress_step = max(1, min(50, velocity_total // 10 if velocity_total > 0 else 1))
    velocity_completed = 0
    for position_prefilter_tau_s in config.vel_position_prefilter_tau_grid_s:
        for tau_s in config.vel_tau_grid_s:
            for slope_window_samples in config.vel_slope_window_grid_samples:
                vel_model = identify_velocity_model(
                    design_df,
                    lookup_steady_df,
                    design_hold_ids,
                    design_ts_s,
                    tau_s,
                    slope_window_samples,
                    position_prefilter_tau_s,
                )
                vel_prediction = evaluate_prediction_metrics(
                    vel_model,
                    evaluation_df,
                    evaluation_hold_ids,
                    lookup_steady_df,
                    evaluation_ts_s,
                )
                if design_df is evaluation_df:
                    vel_eval_model = vel_model
                else:
                    vel_eval_model = identify_velocity_model(
                        evaluation_df,
                        lookup_steady_df,
                        evaluation_hold_ids,
                        evaluation_ts_s,
                        tau_s,
                        slope_window_samples,
                        position_prefilter_tau_s,
                    )
                velocity_tasks = [
                    (
                        (
                            f"vel_qs{base_result.q_s:g}_qp{base_result.q_p:g}_qv{q_v:g}"
                            f"_qi{base_result.q_i:g}_r{base_result.r_value:g}"
                            f"_bd{base_result.branch_deadband_mm:g}_pd{base_result.position_deadband_mm:g}"
                            f"_tau{tau_s:g}_n{slope_window_samples}_sp{position_prefilter_tau_s:g}"
                        ),
                        float(base_result.q_s),
                        float(base_result.q_p),
                        float(base_result.q_i),
                        float(q_v),
                        float(base_result.r_value),
                        float(base_result.branch_deadband_mm),
                        float(base_result.position_deadband_mm),
                    )
                    for base_result in top_base_results
                    for q_v in config.q_v_grid
                ]
                with ProcessPoolExecutor(
                    max_workers=max_workers,
                    mp_context=mp_ctx,
                    initializer=_init_lqi_worker,
                    initargs=(
                        vel_model,
                        vel_eval_model,
                        lookup_maps,
                        reference_scenarios,
                        design_ts_s,
                        evaluation_ts_s,
                        control_limits,
                        lookup_steady_df,
                        config,
                        vel_prediction,
                        transient_profiles,
                    ),
                ) as executor:
                    futures = [executor.submit(_evaluate_lqi_task, task) for task in velocity_tasks]
                    for result in _collect_lqi_results_with_progress(
                        futures=futures,
                        label="LQI velocity-state sweep progress",
                        total=velocity_total,
                        progress_step=velocity_progress_step,
                        initial_completed=velocity_completed,
                    ):
                        velocity_results.append(result)
                        all_results.append(result)
                        velocity_completed += 1

    if velocity_results:
        ranked_velocity, _, _ = _rank_controller_candidates(velocity_results)
        if ranked_velocity:
            best_velocity = ranked_velocity[0]
            family_best[best_velocity.variant] = best_velocity

    return family_best, all_results


def _aggregate_scenario_metrics(
    scenario_metrics: Dict[str, Dict[str, float]]
) -> Dict[str, float]:
    metric_names = list(next(iter(scenario_metrics.values())).keys())
    aggregated: Dict[str, float] = {}
    for name in metric_names:
        values = [metrics[name] for metrics in scenario_metrics.values()]
        if name == "tracking_overshoot_mm":
            aggregated[name] = float(np.nanmax(values))
        elif name in {"IAE_mm_s", "ISE_mm2_s", "ITAE_mm_s2", "command_total_variation_v"}:
            aggregated[name] = float(np.mean(values))
        else:
            aggregated[name] = _safe_nanmean(values)
    return aggregated


def aggregate_scenario_metrics(
    scenario_metrics: Dict[str, Dict[str, float]]
) -> Dict[str, float]:
    return _aggregate_scenario_metrics(scenario_metrics)


def pick_final_variant(results: Dict[str, ControllerResult]) -> str:
    ranked_results, _, _ = _rank_controller_candidates(list(results.values()))
    if not ranked_results:
        raise ValueError("No controller-family candidates available for final selection.")
    return ranked_results[0].variant


def build_control_limits(
    control_window_row: pd.Series,
    excitation_summary_row: pd.Series,
    chosen_result: ControllerResult,
) -> Dict[str, float | str | bool]:
    u_min = float(control_window_row["window_low_v"])
    u_max = float(control_window_row["window_high_v"])
    excitation_low = float(excitation_summary_row["excitation_low_v"])
    excitation_high = float(excitation_summary_row["excitation_high_v"])
    integral_limit = float(chosen_result.integral_state_limit)
    return {
        "u_min_v": u_min,
        "u_max_v": u_max,
        "excitation_low_v": excitation_low,
        "excitation_high_v": excitation_high,
        "integral_state_limit": float(integral_limit),
        "branch_selection_rule": "hold_previous_inside_deadband_then_sign_of_position_error",
        "branch_deadband_mm": chosen_result.branch_deadband_mm,
        "anti_windup_mode": "conditional_freeze",
        "position_deadband_mm": chosen_result.position_deadband_mm,
        "velocity_filter_tau_s": chosen_result.velocity_filter_tau_s,
        "velocity_slope_window_samples": chosen_result.velocity_slope_window_samples,
        "velocity_filter_alpha": chosen_result.velocity_filter_alpha,
        "velocity_position_prefilter_tau_s": chosen_result.velocity_position_prefilter_tau_s,
        "velocity_position_prefilter_alpha": chosen_result.velocity_position_prefilter_alpha,
        "min_abs_integral_gain": chosen_result.min_abs_integral_gain,
        "uses_creep_compensation": chosen_result.uses_creep_compensation,
        "uses_dynamic_pressure_reference": chosen_result.uses_dynamic_pressure_reference,
        "transient_profile_horizon_s": chosen_result.transient_profile_horizon_s,
        "transient_lookup_source": chosen_result.transient_lookup_source,
        "uses_pressure_filter": chosen_result.uses_pressure_filter,
    }


def write_lookup_tables(lookup_tables: Dict[str, pd.DataFrame], output_dir: Path) -> None:
    for branch, table_df in lookup_tables.items():
        table_df.to_csv(output_dir / f"lqr_steady_state_lookup_{branch}.csv", index=False)


def write_feedforward_alignment_outputs(
    alignment_points_df: pd.DataFrame,
    alignment_summary_df: pd.DataFrame,
    output_dir: Path,
) -> None:
    alignment_points_df.to_csv(output_dir / "lqr_feedforward_alignment_points.csv", index=False)
    alignment_summary_df.to_csv(output_dir / "lqr_feedforward_alignment_summary.csv", index=False)
    stale_bands_path = output_dir / "lqr_feedforward_correction_bands.csv"
    if stale_bands_path.exists():
        stale_bands_path.unlink()


def formula_for_observable(name: str, chosen_result: ControllerResult) -> str:
    velocity_window_text = (
        str(int(chosen_result.velocity_slope_window_samples))
        if chosen_result.velocity_slope_window_samples is not None
        else "N"
    )
    formulas = {
        "s_tilde": "s_meas_mm - s_ref_mm",
        "p_tilde": "p_meas_bar - p_ff_bar",
        "v_hat_mm_s": (
            f"alpha_v * v_hat_prev_mm_s + (1 - alpha_v) * "
            f"linear_fit_slope_last_{velocity_window_text}_filtered_position_samples_mm_per_s"
        ),
        "e_I": "e_I_prev + (s_meas_mm - s_ref_mm) * Ts_s",
    }
    return formulas[name]


def uses_prev_sample(name: str) -> bool:
    return name in {"v_hat_mm_s", "e_I"}


def write_observable_definition(output_dir: Path, chosen_result: ControllerResult) -> None:
    names = _feedback_state_names(chosen_result)
    rows = [
        {
            "index": idx + 1,
            "name": name,
            "formula": formula_for_observable(name, chosen_result),
            "uses_prev_sample": uses_prev_sample(name),
            "velocity_slope_window_samples": chosen_result.velocity_slope_window_samples,
            "velocity_position_prefilter_tau_s": chosen_result.velocity_position_prefilter_tau_s,
        }
        for idx, name in enumerate(names)
    ]
    pd.DataFrame(rows).to_csv(output_dir / "lqr_observable_definition.csv", index=False)


def _feedback_gain_rows(result: ControllerResult, include_variant: bool = False) -> List[Dict[str, float | str]]:
    gain_rows: List[Dict[str, float | str]] = []
    for branch, K_branch in result.branch_K_aug.items():
        k_branch_row = np.asarray(K_branch).reshape(1, -1)
        for col_idx in range(k_branch_row.shape[1]):
            row: Dict[str, float | str] = {
                "matrix": f"K_aug_{branch}",
                "row": 1,
                "col": col_idx + 1,
                "value": float(k_branch_row[0, col_idx]),
            }
            if include_variant:
                row["variant"] = result.variant
            gain_rows.append(row)
    return gain_rows


def write_matrix_exports(chosen_result: ControllerResult, output_dir: Path) -> None:
    rows = []
    matrices = {
        "A_model": chosen_result.model.A,
        "B_model": chosen_result.model.B,
        "C_track": chosen_result.C_track,
        "A_lqi_aug": chosen_result.A_aug,
        "B_lqi_aug": chosen_result.B_aug,
    }
    for branch, (A_branch, B_branch) in chosen_result.model.branch_models.items():
        matrices[f"A_model_{branch}"] = A_branch
        matrices[f"B_model_{branch}"] = B_branch
    for branch, A_aug_branch in chosen_result.branch_A_aug.items():
        matrices[f"A_lqi_aug_{branch}"] = A_aug_branch
    for branch, B_aug_branch in chosen_result.branch_B_aug.items():
        matrices[f"B_lqi_aug_{branch}"] = B_aug_branch

    for name, matrix in matrices.items():
        for row_idx in range(matrix.shape[0]):
            for col_idx in range(matrix.shape[1]):
                rows.append(
                    {"matrix": name, "row": row_idx + 1, "col": col_idx + 1, "value": float(matrix[row_idx, col_idx])}
                )
    pd.DataFrame(rows).to_csv(output_dir / "lqr_lifted_model_matrices.csv", index=False)

    pd.DataFrame(_feedback_gain_rows(chosen_result)).to_csv(output_dir / "lqr_feedback_gains.csv", index=False)


def write_nonselected_feedback_gains(
    output_dir: Path,
    controller_results: Dict[str, ControllerResult],
    chosen_variant: str,
) -> None:
    rows: List[Dict[str, float | str]] = []
    for variant, result in controller_results.items():
        if variant == chosen_variant:
            continue
        rows.extend(_feedback_gain_rows(result, include_variant=True))
    pd.DataFrame(rows, columns=["variant", "matrix", "row", "col", "value"]).to_csv(
        output_dir / "lqr_feedback_gains_nonselected.csv", index=False
    )


def write_controller_limits(output_dir: Path, control_limits: Dict[str, float | str | bool]) -> None:
    pd.DataFrame([control_limits]).to_csv(output_dir / "lqr_controller_limits.csv", index=False)


def build_design_summary(
    controller_results: Dict[str, ControllerResult],
    chosen_variant: str,
    ts_s: float,
    steady_df: pd.DataFrame,
    data_df: pd.DataFrame,
    schedule_check: Dict[str, float | int | bool],
    tail_fraction: float,
    input_dataset: str,
) -> pd.DataFrame:
    rows = []
    accepted_steady = int(steady_df["steady_enough"].sum())
    recommended_order = "lookup_plus_branchwise_minimal_LQI>lookup_plus_branchwise_minimal_LQI_VEL"
    ranked_family_results, selection_gate_used, selection_best_tracking_rms_mm = _rank_controller_candidates(
        list(controller_results.values())
    )
    ranked_variants = {result.variant for result in ranked_family_results}
    for variant, result in controller_results.items():
        rows.append(
            {
                "variant": variant,
                "chosen_variant": variant == chosen_variant,
                "input_dataset": input_dataset,
                "data_split_policy": "single_dataset_no_hold_split",
                "controller_family": result.model.controller_family,
                "sampling_time_s": ts_s,
                "n_samples": int(len(data_df)),
                "n_holds_total": int(data_df["hold_id"].nunique()),
                "n_holds_full": int((data_df.groupby("hold_id")["sample_id"].count() > 1).sum()),
                "steady_tail_fraction": float(tail_fraction),
                "n_steady_points_total": int(len(steady_df)),
                "n_steady_points_accepted": accepted_steady,
                "schedule_available": schedule_check["schedule_available"],
                "schedule_match_holds": schedule_check["schedule_match_holds"],
                "schedule_match_commands": schedule_check["schedule_match_commands"],
                "uses_lookup_prefilter": result.model.uses_lookup_prefilter,
                "uses_branchwise_dynamics": True,
                "uses_piecewise_feedforward_correction": False,
                "uses_velocity_state": result.uses_velocity_state,
                "uses_creep_compensation": result.uses_creep_compensation,
                "uses_dynamic_pressure_reference": result.uses_dynamic_pressure_reference,
                "transient_profile_horizon_s": result.transient_profile_horizon_s,
                "transient_lookup_source": result.transient_lookup_source,
                "uses_pressure_filter": result.uses_pressure_filter,
                "stable_model": result.model.stable_model,
                "spectral_radius": result.model.spectral_radius,
                "evaluation_stable_model": result.evaluation_model_stable,
                "evaluation_spectral_radius": result.evaluation_spectral_radius,
                "prediction_dataset_role": "single_dataset_in_sample",
                "selection_dataset_role": "single_dataset_simulation",
                "one_step_rmse_s_mm": result.one_step_rmse_s_mm,
                "one_step_rmse_p_bar": result.one_step_rmse_p_bar,
                "ten_step_rmse_s_mm": result.ten_step_rmse_s_mm,
                "ten_step_rmse_p_bar": result.ten_step_rmse_p_bar,
                "tracking_rms_mm": result.tracking_rms_mm,
                "tracking_overshoot_mm": result.tracking_overshoot_mm,
                "control_effort_mean_abs_v": result.control_effort_mean_abs_v,
                "settled_error_rms_mm": result.settled_error_rms_mm,
                "settled_error_p2p_mm": result.settled_error_p2p_mm,
                "settled_command_p2p_v": result.settled_command_p2p_v,
                "settled_velocity_rms_mm_s": result.settled_velocity_rms_mm_s,
                "rise_time_s_10_90": result.rise_time_s_10_90,
                "settling_time_s_2pct": result.settling_time_s_2pct,
                "peak_time_s": result.peak_time_s,
                "overshoot_mm": result.overshoot_mm,
                "overshoot_pct": result.overshoot_pct,
                "steady_state_bias_mm": result.steady_state_bias_mm,
                "abs_steady_state_bias_mm": result.abs_steady_state_bias_mm,
                "IAE_mm_s": result.iae_mm_s,
                "ISE_mm2_s": result.ise_mm2_s,
                "ITAE_mm_s2": result.itae_mm_s2,
                "control_residual_rms_v": result.control_residual_rms_v,
                "command_total_variation_v": result.command_total_variation_v,
                "saturation_fraction": result.saturation_fraction,
                "integral_state_limit": result.integral_state_limit,
                "validation_ok": result.validation_ok,
                "q_s": result.q_s,
                "q_p": result.q_p,
                "q_v": result.q_v,
                "q_i": result.q_i,
                "branch_deadband_mm": result.branch_deadband_mm,
                "position_deadband_mm": result.position_deadband_mm,
                "q_diag": "|".join(f"{x:.6g}" for x in result.q_diag),
                "r_value": result.r_value,
                "state_names": "|".join(_feedback_state_names(result)),
                "n_runtime_observables": len(_feedback_state_names(result)),
                "gain_l1_norm": result.gain_l1_norm,
                "max_abs_gain_s": result.max_abs_gain_s,
                "max_abs_gain_p": result.max_abs_gain_p,
                "min_abs_integral_gain": result.min_abs_integral_gain,
                "velocity_filter_tau_s": result.velocity_filter_tau_s,
                "velocity_slope_window_samples": result.velocity_slope_window_samples,
                "velocity_position_prefilter_tau_s": result.velocity_position_prefilter_tau_s,
                "velocity_filter_alpha": result.velocity_filter_alpha,
                "velocity_position_prefilter_alpha": result.velocity_position_prefilter_alpha,
                "pressure_filter_tau_s": result.pressure_filter_tau_s,
                "search_grid_boundary_hit": result.search_grid_boundary_hit,
                "search_grid_boundary_fields": result.search_grid_boundary_fields,
                "selection_gate_tracking_factor": SELECTION_TRACKING_GATE_FACTOR,
                "selection_gain_penalty_exponent": SELECTION_GAIN_PENALTY_EXPONENT,
                "selection_score_formula": SELECTION_SCORE_FORMULA,
                "selection_primary_metrics": "|".join(SELECTION_PRIMARY_METRIC_EXPORT_NAMES),
                "selection_tie_break_metrics": "|".join(SELECTION_TIE_BREAK_EXPORT_NAMES),
                "selection_best_tracking_rms_mm": selection_best_tracking_rms_mm,
                "selection_tracking_gate_passed": result.selection_tracking_gate_passed,
                "selection_eligibility_passed": result.selection_eligibility_passed,
                "selection_gate_used": selection_gate_used,
                "selection_score_primary": result.selection_score_primary,
                "selection_ranked_pool_member": variant in ranked_variants,
                "recommended_variant_order": recommended_order,
            }
        )
    return pd.DataFrame(rows)


def build_controller_sweep_metrics(
    all_results: List[ControllerResult],
    family_best_results: Dict[str, ControllerResult],
    chosen_variant: str,
    input_dataset: str,
) -> pd.DataFrame:
    rows: List[Dict[str, float | str | bool]] = []
    selected_config_id = family_best_results[chosen_variant].config_id
    for result in all_results:
        family_best_config_id = family_best_results[result.variant].config_id if result.variant in family_best_results else ""
        for scenario_name, metrics in result.scenario_metrics.items():
            row = {
                "variant": result.variant,
                "controller_family": result.model.controller_family,
                "config_id": result.config_id,
                "scenario": scenario_name,
                "chosen_variant": result.variant == chosen_variant,
                "family_best": result.config_id == family_best_config_id,
                "global_selected": result.config_id == selected_config_id,
                "input_dataset": input_dataset,
                "data_split_policy": "single_dataset_no_hold_split",
                "prediction_dataset_role": "single_dataset_in_sample",
                "selection_dataset_role": "single_dataset_simulation",
                "uses_velocity_state": result.uses_velocity_state,
                "uses_creep_compensation": result.uses_creep_compensation,
                "uses_dynamic_pressure_reference": result.uses_dynamic_pressure_reference,
                "transient_profile_horizon_s": result.transient_profile_horizon_s,
                "transient_lookup_source": result.transient_lookup_source,
                "velocity_filter_tau_s": result.velocity_filter_tau_s,
                "velocity_slope_window_samples": result.velocity_slope_window_samples,
                "velocity_position_prefilter_tau_s": result.velocity_position_prefilter_tau_s,
                "velocity_position_prefilter_alpha": result.velocity_position_prefilter_alpha,
                "q_s": result.q_s,
                "q_p": result.q_p,
                "q_v": result.q_v,
                "q_i": result.q_i,
                "r_value": result.r_value,
                "branch_deadband_mm": result.branch_deadband_mm,
                "position_deadband_mm": result.position_deadband_mm,
                "n_runtime_observables": len(_feedback_state_names(result)),
                "stable_model": result.model.stable_model,
                "validation_ok": result.validation_ok,
                "one_step_rmse_s_mm": result.one_step_rmse_s_mm,
                "one_step_rmse_p_bar": result.one_step_rmse_p_bar,
                "ten_step_rmse_s_mm": result.ten_step_rmse_s_mm,
                "ten_step_rmse_p_bar": result.ten_step_rmse_p_bar,
                "gain_l1_norm": result.gain_l1_norm,
                "max_abs_gain_s": result.max_abs_gain_s,
                "max_abs_gain_p": result.max_abs_gain_p,
                "min_abs_integral_gain": result.min_abs_integral_gain,
                "integral_state_limit": result.integral_state_limit,
                "abs_steady_state_bias_mm": result.abs_steady_state_bias_mm,
                "search_grid_boundary_hit": result.search_grid_boundary_hit,
                "search_grid_boundary_fields": result.search_grid_boundary_fields,
                "selection_gain_penalty_exponent": SELECTION_GAIN_PENALTY_EXPONENT,
                "selection_score_formula": SELECTION_SCORE_FORMULA,
                "selection_primary_metrics": "|".join(SELECTION_PRIMARY_METRIC_EXPORT_NAMES),
                "selection_tie_break_metrics": "|".join(SELECTION_TIE_BREAK_EXPORT_NAMES),
                "selection_score_primary": result.selection_score_primary,
                "selection_tracking_gate_passed": result.selection_tracking_gate_passed,
                "selection_eligibility_passed": result.selection_eligibility_passed,
            }
            row.update(metrics)
            rows.append(row)
    return pd.DataFrame(rows)
