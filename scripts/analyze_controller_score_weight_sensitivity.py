#!/home/tomislav/bin/micromamba/prefix/envs/pykmd_mkl/bin/python

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TRACKING_GATE_FACTOR = 1.12
MIN_ABS_INTEGRAL_GAIN = 0.001
GAIN_PENALTY_EXPONENT = 3.0
BIAS_NORMALIZATION_FLOOR = 0.01
WEIGHT_PERTURBATION_FACTORS = (0.75, 1.25)
LEGACY_GAIN_WEIGHTS = (0.10, 0.30)
ALTERNATIVE_GAIN_EXPONENTS = (1.0, 2.0, 4.0)
CONTROLLER_LABELS = {
    "lookup_plus_branchwise_minimal_LQI": "Base LQI",
    "lookup_plus_branchwise_minimal_LQI_VEL": "Velocity-state LQI",
    "lookup_plus_branchwise_PID": "PID",
    "lookup_only_feedforward": "Feedforward-only",
}
COMMON_GAIN_METRICS = {
    "max_abs_gain_s": "Maximum displacement gain",
    "gain_l1_norm": "Gain L1 norm",
    "min_abs_integral_gain": "Minimum integral-gain magnitude",
}

PRIMARY_METRICS = {
    "abs_steady_state_bias_mm": "abs_steady_state_bias_mm",
    "settled_error_p2p_mm": "settled_error_p2p_mm",
    "itae_mm_s2": "ITAE_mm_s2",
    "tracking_rms_mm": "tracking_rms_mm",
    "overshoot_pct": "overshoot_pct",
    "settling_time_s_2pct": "settling_time_s_2pct",
    "settled_command_p2p_v": "settled_command_p2p_v",
    "settled_velocity_rms_mm_s": "settled_velocity_rms_mm_s",
    "max_abs_gain_s": "max_abs_gain_s",
}
BASE_WEIGHTS = {
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
IMPLEMENTED_WEIGHT_SUM = sum(BASE_WEIGHTS.values())
IMPLEMENTED_NON_GAIN_SUM = IMPLEMENTED_WEIGHT_SUM - BASE_WEIGHTS["max_abs_gain_s"]
TIE_BREAK_COLUMNS = [
    "saturation_fraction",
    "command_total_variation_v",
    "gain_l1_norm",
]
REQUIRED_COLUMNS = [
    "variant",
    "config_id",
    "scenario",
    "validation_ok",
    "min_abs_integral_gain",
    *PRIMARY_METRICS.values(),
    *TIE_BREAK_COLUMNS,
]
OPTIONAL_COLUMNS = [
    "controller_family",
    "stage",
    "global_selected",
    "family_best",
    "q_s",
    "q_p",
    "q_v",
    "q_i",
    "r_value",
    "branch_deadband_mm",
    "position_deadband_mm",
    "integral_state_limit",
    "max_abs_gain_p",
    "kp_up_v_per_mm",
    "ki_up_v_per_mm_s",
    "kd_up_v_s_per_mm",
    "kp_down_v_per_mm",
    "ki_down_v_per_mm_s",
    "kd_down_v_s_per_mm",
    "tau_d_s",
    "slope_window_samples",
    "position_prefilter_tau_s",
    "velocity_filter_tau_s",
    "velocity_slope_window_samples",
    "velocity_position_prefilter_tau_s",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Re-score existing aggregate controller sweeps at alternative gain weights "
            "without rerunning simulations."
        )
    )
    parser.add_argument(
        "--lqr_csv",
        default="analysis_outputs/lqr_controller_sweep_metrics.csv",
        help="Saved LQI sweep metrics CSV.",
    )
    parser.add_argument(
        "--pid_csv",
        default="analysis_outputs/pid_controller_sweep_metrics.csv",
        help="Saved PID sweep metrics CSV.",
    )
    parser.add_argument(
        "--output_csv",
        default="analysis_outputs/controller_score_weight_sensitivity.csv",
        help="Sensitivity-result CSV to write.",
    )
    parser.add_argument(
        "--gain_stats_csv",
        default="analysis_outputs/controller_score_weight_sensitivity_gain_statistics.csv",
        help="Selected-gain statistics CSV to write.",
    )
    parser.add_argument(
        "--selection_summary_csv",
        default="analysis_outputs/controller_score_weight_sensitivity_summary.csv",
        help="Selection-sensitivity summary CSV to write.",
    )
    parser.add_argument(
        "--plot_path",
        default="analysis_outputs/controller_score_weight_sensitivity_boxplots.png",
        help="Reviewer-only boxplot to write in analysis_outputs.",
    )
    parser.add_argument(
        "--reviewer_plot",
        default="article/revision1/controller_score_weight_sensitivity_boxplots.png",
        help="Reviewer-only copy of the boxplot.",
    )
    return parser.parse_args()


def _as_bool(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.strip().str.lower()
    parsed = values.map(
        {
            "true": True,
            "false": False,
            "1": True,
            "0": False,
            "yes": True,
            "no": False,
        }
    )
    if parsed.isna().any():
        bad = sorted(values.loc[parsed.isna()].unique())
        raise ValueError(f"Unrecognized boolean values: {bad}")
    return parsed.astype(bool)


def load_aggregate_rows(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing sweep file: {path}")

    header = pd.read_csv(path, nrows=0)
    missing = [column for column in REQUIRED_COLUMNS if column not in header.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    usecols = [
        column
        for column in [*REQUIRED_COLUMNS, *OPTIONAL_COLUMNS]
        if column in header.columns
    ]

    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=50_000):
        aggregate = chunk.loc[chunk["scenario"].astype(str).eq("all")].copy()
        if not aggregate.empty:
            chunks.append(aggregate)
    if not chunks:
        raise ValueError(f"No scenario=all rows found in {path}")

    result = pd.concat(chunks, ignore_index=True)
    result["config_id"] = result["config_id"].astype(str)
    result["validation_ok"] = _as_bool(result["validation_ok"])
    for column in [*PRIMARY_METRICS.values(), "min_abs_integral_gain", *TIE_BREAK_COLUMNS]:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    for column in ["global_selected", "family_best"]:
        if column in result:
            result[column] = _as_bool(result[column])
        else:
            result[column] = False

    duplicate_ids = result.loc[result["config_id"].duplicated(), "config_id"]
    if not duplicate_ids.empty:
        raise ValueError(
            f"{path} contains duplicate aggregate configuration IDs, including "
            f"{duplicate_ids.iloc[0]!r}"
        )
    return result


def build_weights(
    perturbed_metric: str | None = None,
    factor: float = 1.0,
    gain_weight: float | None = None,
) -> dict[str, float]:
    if perturbed_metric is None and gain_weight is None:
        weights = dict(BASE_WEIGHTS)
    else:
        if gain_weight is not None:
            perturbed_metric = "max_abs_gain_s"
            target_weight = gain_weight
        else:
            if perturbed_metric not in BASE_WEIGHTS:
                raise ValueError(f"Unknown score weight: {perturbed_metric}")
            target_weight = BASE_WEIGHTS[perturbed_metric] * factor
        if not 0.0 < target_weight < IMPLEMENTED_WEIGHT_SUM:
            raise ValueError("Perturbed score weight must be positive and below the total scale.")
        other_sum = IMPLEMENTED_WEIGHT_SUM - BASE_WEIGHTS[perturbed_metric]
        other_scale = (IMPLEMENTED_WEIGHT_SUM - target_weight) / other_sum
        weights = {
            name: value * other_scale
            for name, value in BASE_WEIGHTS.items()
            if name != perturbed_metric
        }
        weights[perturbed_metric] = target_weight
    if not np.isclose(sum(weights.values()), IMPLEMENTED_WEIGHT_SUM):
        raise AssertionError("Sensitivity coefficients changed their total scale.")
    return weights


def build_profiles() -> list[dict[str, object]]:
    profiles: list[dict[str, object]] = [
        {
            "profile_id": "baseline",
            "profile_type": "baseline",
            "perturbed_metric": "",
            "perturbation_factor": 1.0,
            "gain_weight": BASE_WEIGHTS["max_abs_gain_s"],
            "gain_penalty_exponent": GAIN_PENALTY_EXPONENT,
            "weights": build_weights(),
        }
    ]
    for metric in BASE_WEIGHTS:
        for factor in WEIGHT_PERTURBATION_FACTORS:
            profiles.append(
                {
                    "profile_id": f"weight_{metric}_x{factor:.2f}",
                    "profile_type": "one_factor_at_a_time",
                    "perturbed_metric": metric,
                    "perturbation_factor": factor,
                    "gain_weight": np.nan,
                    "gain_penalty_exponent": GAIN_PENALTY_EXPONENT,
                    "weights": build_weights(metric, factor=factor),
                }
            )
    for gain_weight in LEGACY_GAIN_WEIGHTS:
        profiles.append(
            {
                "profile_id": f"gain_weight_{gain_weight:.2f}",
                "profile_type": "legacy_gain_stress",
                "perturbed_metric": "max_abs_gain_s",
                "perturbation_factor": gain_weight / BASE_WEIGHTS["max_abs_gain_s"],
                "gain_weight": gain_weight,
                "gain_penalty_exponent": GAIN_PENALTY_EXPONENT,
                "weights": build_weights(gain_weight=gain_weight),
            }
        )
    for exponent in ALTERNATIVE_GAIN_EXPONENTS:
        profiles.append(
            {
                "profile_id": f"gain_exponent_{exponent:.0f}",
                "profile_type": "gain_exponent",
                "perturbed_metric": "max_abs_gain_s",
                "perturbation_factor": np.nan,
                "gain_weight": BASE_WEIGHTS["max_abs_gain_s"],
                "gain_penalty_exponent": exponent,
                "weights": build_weights(),
            }
        )
    if len(profiles) != 24:
        raise AssertionError(f"Expected 24 sensitivity profiles, found {len(profiles)}")
    return profiles


def _score_rows(
    rows: pd.DataFrame,
    weights: dict[str, float],
    gain_penalty_exponent: float,
) -> tuple[pd.DataFrame, bool, float]:
    if gain_penalty_exponent <= 0.0:
        raise ValueError("Gain penalty exponent must be positive.")
    rows = rows.copy()
    finite_primary = rows[list(PRIMARY_METRICS.values())].notna().all(axis=1)
    finite_tracking = rows["tracking_rms_mm"].replace([np.inf, -np.inf], np.nan).notna()
    finite_tracking_values = rows.loc[finite_tracking, "tracking_rms_mm"]
    if finite_tracking_values.empty:
        best_tracking = float("inf")
    else:
        best_tracking = float(finite_tracking_values.min())
    tracking_passed = rows["tracking_rms_mm"] <= TRACKING_GATE_FACTOR * best_tracking
    eligible = (
        rows["validation_ok"]
        & finite_primary
        & tracking_passed
        & (rows["min_abs_integral_gain"] >= MIN_ABS_INTEGRAL_GAIN)
    )
    gate_used = bool(eligible.any())
    ranked_pool = rows.loc[eligible].copy() if gate_used else rows.copy()

    metric_mins = {
        name: float(ranked_pool[column].dropna().min())
        if ranked_pool[column].notna().any()
        else float("inf")
        for name, column in PRIMARY_METRICS.items()
    }
    scores = np.zeros(len(ranked_pool), dtype=float)
    score_valid = np.ones(len(ranked_pool), dtype=bool)
    for name, column in PRIMARY_METRICS.items():
        values = ranked_pool[column].to_numpy(dtype=float)
        minimum = metric_mins[name]
        floor = BIAS_NORMALIZATION_FLOOR if name == "abs_steady_state_bias_mm" else 0.0
        denominator = max(abs(minimum), floor, 1e-9)
        normalized = values / denominator
        if name == "max_abs_gain_s":
            normalized = normalized**gain_penalty_exponent
        score_valid &= np.isfinite(normalized)
        scores += weights[name] * np.where(np.isfinite(normalized), normalized, 0.0)
    scores[~score_valid] = np.inf
    ranked_pool["_score"] = scores
    ranked_pool = ranked_pool.sort_values(
        ["_score", *TIE_BREAK_COLUMNS],
        ascending=True,
        kind="mergesort",
    ).reset_index(drop=True)
    ranked_pool["_rank"] = np.arange(1, len(ranked_pool) + 1)
    return ranked_pool, gate_used, best_tracking


def _current_id(rows: pd.DataFrame, marker: str) -> str:
    marked = rows.loc[rows[marker], "config_id"]
    if len(marked) != 1:
        raise ValueError(f"Expected exactly one {marker}=True aggregate row, found {len(marked)}")
    return str(marked.iloc[0])


SELECTED_FIELD_MAP = {
    "selected_q_s": "q_s",
    "selected_q_p": "q_p",
    "selected_q_v": "q_v",
    "selected_q_i": "q_i",
    "selected_r_value": "r_value",
    "selected_branch_deadband_mm": "branch_deadband_mm",
    "selected_position_deadband_mm": "position_deadband_mm",
    "selected_integral_state_limit": "integral_state_limit",
    "selected_gain_l1_norm": "gain_l1_norm",
    "selected_max_abs_gain_s": "max_abs_gain_s",
    "selected_max_abs_gain_p": "max_abs_gain_p",
    "selected_min_abs_integral_gain": "min_abs_integral_gain",
    "selected_kp_up_v_per_mm": "kp_up_v_per_mm",
    "selected_ki_up_v_per_mm_s": "ki_up_v_per_mm_s",
    "selected_kd_up_v_s_per_mm": "kd_up_v_s_per_mm",
    "selected_kp_down_v_per_mm": "kp_down_v_per_mm",
    "selected_ki_down_v_per_mm_s": "ki_down_v_per_mm_s",
    "selected_kd_down_v_s_per_mm": "kd_down_v_s_per_mm",
    "selected_tau_d_s": "tau_d_s",
    "selected_slope_window_samples": "slope_window_samples",
    "selected_position_prefilter_tau_s": "position_prefilter_tau_s",
    "selected_velocity_filter_tau_s": "velocity_filter_tau_s",
    "selected_velocity_slope_window_samples": "velocity_slope_window_samples",
    "selected_velocity_position_prefilter_tau_s": "velocity_position_prefilter_tau_s",
}


def _profile_fields(profile: dict[str, object]) -> dict[str, object]:
    weights = profile["weights"]
    if not isinstance(weights, dict):
        raise TypeError("Sensitivity profile weights must be a dictionary.")
    output: dict[str, object] = {
        "profile_id": profile["profile_id"],
        "profile_type": profile["profile_type"],
        "perturbed_metric": profile["perturbed_metric"],
        "perturbation_factor": profile["perturbation_factor"],
        "gain_weight": profile["gain_weight"],
        "gain_penalty_exponent": profile["gain_penalty_exponent"],
        "weight_sum": round(sum(weights.values()), 12),
    }
    for name, value in weights.items():
        output[f"weight_{name}"] = value
    return output


def sensitivity_row(
    pool_name: str,
    source_name: str,
    rows: pd.DataFrame,
    profile: dict[str, object],
    baseline_config_id: str,
    family_base_config_id: str = "",
    family_velocity_config_id: str = "",
) -> dict[str, object]:
    exponent = float(profile["gain_penalty_exponent"])
    weights = profile["weights"]
    if not isinstance(weights, dict):
        raise TypeError("Sensitivity profile weights must be a dictionary.")
    ranked, gate_used, best_tracking = _score_rows(rows, weights, exponent)
    selected = ranked.iloc[0]
    baseline_ranked = ranked.loc[ranked["config_id"].eq(baseline_config_id), "_rank"]
    eligible_mask = (
        rows["validation_ok"]
        & rows[list(PRIMARY_METRICS.values())].notna().all(axis=1)
        & (rows["tracking_rms_mm"] <= TRACKING_GATE_FACTOR * best_tracking)
        & (rows["min_abs_integral_gain"] >= MIN_ABS_INTEGRAL_GAIN)
    )
    output: dict[str, object] = {
        **_profile_fields(profile),
        "pool": pool_name,
        "source_csv": source_name,
        "analysis_status": "scored_feedback_pool",
        "candidate_count": len(rows),
        "eligible_count": int(eligible_mask.sum()),
        "selection_gate_used": gate_used,
        "best_tracking_rms_mm": best_tracking,
        "selected_controller_label": CONTROLLER_LABELS.get(
            str(selected["variant"]), str(selected["variant"])
        ),
        "selected_variant": selected["variant"],
        "selected_config_id": selected["config_id"],
        "selected_config_rank": 1,
        "selected_score": float(selected["_score"]),
        "selected_tracking_rms_mm": float(selected["tracking_rms_mm"]),
        "baseline_config_id": baseline_config_id,
        "baseline_config_rank": (
            int(baseline_ranked.iloc[0]) if not baseline_ranked.empty else np.nan
        ),
        "selection_same_as_baseline": selected["config_id"] == baseline_config_id,
        "family_base_config_id": family_base_config_id,
        "family_velocity_config_id": family_velocity_config_id,
        "selected_gain_data_available": True,
    }
    for output_name, source_name in SELECTED_FIELD_MAP.items():
        output[output_name] = (
            selected[source_name]
            if source_name in selected.index
            else np.nan
        )
    return output


def feedforward_row(profile: dict[str, object]) -> dict[str, object]:
    return {
        **_profile_fields(profile),
        "pool": "feedforward_only_reference",
        "source_csv": "not_applicable_no_feedback_sweep",
        "analysis_status": "not_scored_no_feedback_gains",
        "candidate_count": np.nan,
        "eligible_count": np.nan,
        "selection_gate_used": np.nan,
        "best_tracking_rms_mm": np.nan,
        "selected_controller_label": CONTROLLER_LABELS["lookup_only_feedforward"],
        "selected_variant": "lookup_only_feedforward",
        "selected_config_id": "",
        "selected_config_rank": np.nan,
        "selected_score": np.nan,
        "selected_tracking_rms_mm": np.nan,
        "baseline_config_id": "",
        "baseline_config_rank": np.nan,
        "selection_same_as_baseline": np.nan,
        "family_base_config_id": "",
        "family_velocity_config_id": "",
        "selected_gain_data_available": False,
    }

def build_gain_observations(result_df: pd.DataFrame) -> pd.DataFrame:
    pool_labels = {
        "lqi_base_pool": "Base LQI",
        "lqi_velocity_pool": "Velocity-state LQI",
        "pid_pool": "PID",
    }
    feedback_rows = result_df.loc[result_df["pool"].isin(pool_labels)].copy()
    observations: list[dict[str, object]] = []
    for pool, controller_label in pool_labels.items():
        pool_rows = feedback_rows.loc[feedback_rows["pool"].eq(pool)]
        baseline_rows = pool_rows.loc[pool_rows["profile_id"].eq("baseline")]
        if len(baseline_rows) != 1:
            raise ValueError(f"Expected one baseline row for {pool}.")
        baseline = baseline_rows.iloc[0]
        for _, row in pool_rows.iterrows():
            for metric, metric_label in COMMON_GAIN_METRICS.items():
                value = pd.to_numeric(row[f"selected_{metric}"], errors="coerce")
                baseline_value = pd.to_numeric(
                    baseline[f"selected_{metric}"], errors="coerce"
                )
                if not np.isfinite(value) or not np.isfinite(baseline_value):
                    continue
                if baseline_value == 0.0:
                    raise ValueError(f"Zero baseline for gain metric {metric} in {pool}.")
                observations.append(
                    {
                        "profile_id": row["profile_id"],
                        "profile_type": row["profile_type"],
                        "perturbed_metric": row["perturbed_metric"],
                        "pool": pool,
                        "controller": controller_label,
                        "gain_metric": metric,
                        "gain_metric_label": metric_label,
                        "selected_value": float(value),
                        "baseline_value": float(baseline_value),
                        "relative_to_baseline": float(value / baseline_value),
                    }
                )
    result = pd.DataFrame(observations)
    if result.empty:
        raise ValueError("No selected feedback gain observations were produced.")
    return result


def build_gain_statistics(observations: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (controller, metric, metric_label), group in observations.groupby(
        ["controller", "gain_metric", "gain_metric_label"], sort=False
    ):
        ratios = group["relative_to_baseline"].to_numpy(dtype=float)
        baseline_values = group["baseline_value"].to_numpy(dtype=float)
        rows.append(
            {
                "controller": controller,
                "gain_metric": metric,
                "gain_metric_label": metric_label,
                "n_profiles": len(group),
                "baseline_value": float(baseline_values[0]),
                "relative_min": float(np.min(ratios)),
                "relative_q1": float(np.quantile(ratios, 0.25)),
                "relative_median": float(np.median(ratios)),
                "relative_q3": float(np.quantile(ratios, 0.75)),
                "relative_max": float(np.max(ratios)),
                "relative_mean": float(np.mean(ratios)),
            }
        )
    return pd.DataFrame(rows)


def build_selection_summary(result_df: pd.DataFrame) -> pd.DataFrame:
    pool_labels = {
        "lqi_base_pool": "Base LQI",
        "lqi_velocity_pool": "Velocity-state LQI",
        "pid_pool": "PID",
    }
    summary_rows: list[dict[str, object]] = []
    for pool, controller_label in pool_labels.items():
        group = result_df.loc[result_df["pool"].eq(pool)].copy()
        baseline = group.loc[group["profile_id"].eq("baseline")].iloc[0]
        summary_rows.append(
            {
                "controller": controller_label,
                "analysis_scope": pool,
                "analysis_status": "scored_feedback_pool",
                "profile_count": len(group),
                "baseline_config_id": baseline["selected_config_id"],
                "baseline_selection_same_count": int(group["selection_same_as_baseline"].sum()),
                "selection_change_count": int(group["selection_same_as_baseline"].eq(False).sum()),
                "unique_selected_config_count": int(group["selected_config_id"].nunique()),
                "max_abs_gain_s_ratio_min": float(
                    group["selected_max_abs_gain_s"].min()
                    / baseline["selected_max_abs_gain_s"]
                ),
                "max_abs_gain_s_ratio_max": float(
                    group["selected_max_abs_gain_s"].max()
                    / baseline["selected_max_abs_gain_s"]
                ),
                "gain_l1_ratio_min": float(
                    group["selected_gain_l1_norm"].min()
                    / baseline["selected_gain_l1_norm"]
                ),
                "gain_l1_ratio_max": float(
                    group["selected_gain_l1_norm"].max()
                    / baseline["selected_gain_l1_norm"]
                ),
            }
        )
    family = result_df.loc[result_df["pool"].eq("lqi_family_comparison")].copy()
    summary_rows.append(
        {
            "controller": "LQI family comparison",
            "analysis_scope": "lqi_family_comparison",
            "analysis_status": "scored_family_comparison",
            "profile_count": len(family),
            "baseline_config_id": family.loc[family["profile_id"].eq("baseline"), "selected_config_id"].iloc[0],
            "baseline_selection_same_count": int(family["selection_same_as_baseline"].sum()),
            "selection_change_count": int(family["selection_same_as_baseline"].eq(False).sum()),
            "unique_selected_config_count": int(family["selected_config_id"].nunique()),
            "lqi_base_selection_count": int(
                family["selected_controller_label"].eq("Base LQI").sum()
            ),
            "lqi_velocity_selection_count": int(
                family["selected_controller_label"].eq("Velocity-state LQI").sum()
            ),
        }
    )
    summary_rows.append(
        {
            "controller": "Feedforward-only",
            "analysis_scope": "feedforward_only_reference",
            "analysis_status": "not_scored_no_feedback_gains",
            "profile_count": result_df["profile_id"].nunique(),
            "baseline_config_id": "",
            "baseline_selection_same_count": np.nan,
            "selection_change_count": np.nan,
            "unique_selected_config_count": np.nan,
            "note": "No feedback gains or feedback-family selection score.",
        }
    )
    return pd.DataFrame(summary_rows)


def plot_gain_boxplots(observations: pd.DataFrame, output_path: Path) -> None:
    controller_order = ["Base LQI", "Velocity-state LQI", "PID"]
    metric_order = list(COMMON_GAIN_METRICS)
    colors = ["#4477AA", "#CC6677", "#228833"]
    fig, axes = plt.subplots(
        nrows=len(metric_order),
        ncols=1,
        figsize=(7.2, 8.8),
        sharex=True,
    )
    for panel_index, metric in enumerate(metric_order):
        ax = axes[panel_index]
        groups = [
            observations.loc[
                observations["controller"].eq(controller)
                & observations["gain_metric"].eq(metric),
                "relative_to_baseline",
            ].to_numpy(dtype=float)
            for controller in controller_order
        ]
        boxplot = ax.boxplot(
            groups,
            positions=np.arange(1, len(controller_order) + 1),
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.0},
            whiskerprops={"color": "#555555", "linewidth": 0.8},
            capprops={"color": "#555555", "linewidth": 0.8},
        )
        for patch, color in zip(boxplot["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.65)
            patch.set_edgecolor("#333333")
            patch.set_linewidth(0.8)
        for position, (values, color) in enumerate(zip(groups, colors), start=1):
            offsets = np.linspace(-0.18, 0.18, len(values))
            ax.scatter(
                position + offsets,
                values,
                s=10,
                color=color,
                alpha=0.45,
                edgecolors="none",
                zorder=2,
            )
        for position in range(1, len(controller_order) + 1):
            ax.plot(
                position,
                1.0,
                marker="D",
                markersize=4.5,
                color="black",
                markerfacecolor="white",
                linestyle="none",
                zorder=3,
            )
        ax.axhline(1.0, color="#555555", linewidth=0.8, linestyle="--")
        ax.set_ylabel(f"{COMMON_GAIN_METRICS[metric]}\n(relative to baseline)", fontsize=9)
        ax.tick_params(axis="both", labelsize=9)
        ax.grid(axis="y", linewidth=0.5, alpha=0.35)
        ax.text(
            0.01,
            0.92,
            f"({chr(ord('a') + panel_index)})",
            transform=ax.transAxes,
            fontsize=10,
            va="top",
        )
        finite_values = np.concatenate(groups)
        lower = float(np.min(finite_values))
        upper = float(np.max(finite_values))
        margin = max(0.03, 0.08 * (upper - lower))
        ax.set_ylim(lower - margin, upper + margin)
    axes[-1].set_xticks(np.arange(1, len(controller_order) + 1))
    axes[-1].set_xticklabels(controller_order, fontsize=9)
    axes[-1].set_xlabel("Controller variant", fontsize=10)
    fig.text(
        0.5,
        0.01,
        "Feedforward-only has no feedback gains and is not included in the boxes. "
        "Boxes summarize deterministic sensitivity profiles, not uncertainty intervals.",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.19, right=0.98, top=0.99, bottom=0.12, hspace=0.18)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    lqr_path = Path(args.lqr_csv)
    pid_path = Path(args.pid_csv)
    lqr = load_aggregate_rows(lqr_path)
    pid = load_aggregate_rows(pid_path)

    base = lqr.loc[lqr["variant"].eq("lookup_plus_branchwise_minimal_LQI")].copy()
    velocity = lqr.loc[
        lqr["variant"].eq("lookup_plus_branchwise_minimal_LQI_VEL")
    ].copy()
    if base.empty or velocity.empty:
        raise ValueError("Both LQI variants must be present in the aggregate sweep.")

    current_base_id = _current_id(base, "family_best")
    current_velocity_id = _current_id(velocity, "family_best")
    current_lqi_id = _current_id(lqr, "global_selected")
    current_pid_id = _current_id(pid, "global_selected")

    rows: list[dict[str, object]] = []
    for profile in build_profiles():
        rows.append(
            sensitivity_row(
                "lqi_base_pool",
                str(lqr_path),
                base,
                profile,
                current_base_id,
            )
        )
        rows.append(
            sensitivity_row(
                "lqi_velocity_pool",
                str(lqr_path),
                velocity,
                profile,
                current_velocity_id,
            )
        )
        base_winner = _score_rows(
            base,
            profile["weights"],
            float(profile["gain_penalty_exponent"]),
        )[0].iloc[0]
        velocity_winner = _score_rows(
            velocity,
            profile["weights"],
            float(profile["gain_penalty_exponent"]),
        )[0].iloc[0]
        family_rows = pd.DataFrame([base_winner, velocity_winner])
        rows.append(
            sensitivity_row(
                "lqi_family_comparison",
                str(lqr_path),
                family_rows,
                profile,
                current_lqi_id,
                family_base_config_id=str(base_winner["config_id"]),
                family_velocity_config_id=str(velocity_winner["config_id"]),
            )
        )
        rows.append(
            sensitivity_row(
                "pid_pool",
                str(pid_path),
                pid,
                profile,
                current_pid_id,
            )
        )
        rows.append(feedforward_row(profile))

    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df = pd.DataFrame(rows)
    result_df.to_csv(output_path, index=False)

    observations = build_gain_observations(result_df)
    gain_stats = build_gain_statistics(observations)
    gain_stats_path = Path(args.gain_stats_csv)
    gain_stats_path.parent.mkdir(parents=True, exist_ok=True)
    gain_stats.to_csv(gain_stats_path, index=False)

    selection_summary = build_selection_summary(result_df)
    summary_path = Path(args.selection_summary_csv)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    selection_summary.to_csv(summary_path, index=False)

    plot_path = Path(args.plot_path)
    plot_gain_boxplots(observations, plot_path)
    reviewer_plot_path = Path(args.reviewer_plot)
    reviewer_plot_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(plot_path, reviewer_plot_path)

    print(f"Wrote {len(result_df)} sensitivity rows to {output_path}")
    print(f"Wrote gain statistics to {gain_stats_path}")
    print(f"Wrote selection summary to {summary_path}")
    print(f"Wrote reviewer boxplot to {plot_path} and {reviewer_plot_path}")


if __name__ == "__main__":
    main()
