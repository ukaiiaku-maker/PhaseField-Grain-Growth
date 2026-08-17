from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _front_ranks(objectives: np.ndarray) -> np.ndarray:
    """Return one-based nondominated ranks for minimization objectives."""
    ranks = np.zeros(len(objectives), dtype=int)
    remaining = np.arange(len(objectives))
    rank = 1
    while len(remaining):
        front = []
        for candidate in remaining:
            others = remaining[remaining != candidate]
            dominated = any(
                np.all(objectives[other] <= objectives[candidate])
                and np.any(objectives[other] < objectives[candidate])
                for other in others
            )
            if not dominated:
                front.append(candidate)
        front_array = np.asarray(front, dtype=int)
        ranks[front_array] = rank
        remaining = remaining[~np.isin(remaining, front_array)]
        rank += 1
    return ranks


def pareto_summary(summary: pd.DataFrame, target_exponent: float = 2.0,
                   exponent_tolerance: float = 0.5) -> pd.DataFrame:
    """Rank jerkiness candidates subject to an ensemble-scaling constraint.

    Scaling error is minimized while trajectory CV, event-count Fano factor,
    and waiting-time burstiness are maximized. Missing intermittency objectives
    are assigned the worst observed value rather than silently treated as zero.
    """
    required = {
        "regime", "temperature", "n", "n_ci_low", "n_ci_high", "K",
        "jerkiness_CV", "Fano", "burstiness",
    }
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"summary lacks Pareto columns: {sorted(missing)}")
    result = summary.copy()
    result["target_exponent"] = float(target_exponent)
    result["scaling_error"] = np.abs(result["n"] - target_exponent)
    result["ci_contains_target"] = (
        (result["n_ci_low"] <= target_exponent)
        & (result["n_ci_high"] >= target_exponent)
    )
    result["physically_admissible"] = (
        np.isfinite(result["n"])
        & np.isfinite(result["K"])
        & (result["K"] > 0)
        & (
            (result["scaling_error"] <= exponent_tolerance)
            | result["ci_contains_target"]
        )
    )
    result["pareto_rank"] = pd.array([pd.NA] * len(result), dtype="Int64")
    eligible = result.index[result["physically_admissible"]].to_numpy()
    if len(eligible):
        objective_columns = ["scaling_error", "jerkiness_CV", "Fano", "burstiness"]
        objective_data = result.loc[eligible, objective_columns].to_numpy(float)
        for column in range(1, objective_data.shape[1]):
            finite = np.isfinite(objective_data[:, column])
            worst = np.min(objective_data[finite, column]) - 1.0 if np.any(finite) else -1.0
            objective_data[~finite, column] = worst
            objective_data[:, column] *= -1.0
        result.loc[eligible, "pareto_rank"] = _front_ranks(objective_data)
    result["pareto_front"] = result["pareto_rank"] == 1
    return result.sort_values(
        ["physically_admissible", "pareto_rank", "scaling_error", "jerkiness_CV"],
        ascending=[False, True, True, False], na_position="last",
    ).reset_index(drop=True)


def write_pareto_summary(summary_path: str | Path, output_path: str | Path,
                         target_exponent: float = 2.0,
                         exponent_tolerance: float = 0.5) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    ranked = pareto_summary(
        pd.read_csv(summary_path), target_exponent, exponent_tolerance
    )
    ranked.to_csv(output, index=False)
    return output
