from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import fit_growth_law
from grain_growth_pf.analysis.jerkiness import jerkiness_metrics


SUMMARY_COLUMNS = [
    "regime", "temperature", "n", "n_ci_low", "n_ci_high", "K", "K_ci",
    "Q_app", "Q_app_ci", "jerkiness_CV", "Fano", "burstiness",
    "reverse_motion_fraction", "velocity_curvature_R2", "pinned_fraction",
    "number_of_events", "number_of_realizations", "Git_SHA",
]


def analyze_run(run_dir: str | Path) -> dict[str, object]:
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    config = manifest["config"]
    tracks = ensemble_radius(load_tracks(run_dir / "grain_tracks.csv"))
    fit = fit_growth_law(tracks["time"].to_numpy(), tracks["R_A"].to_numpy(), transient_fraction=0.15)
    per_grain = load_tracks(run_dir / "grain_tracks.csv")
    metrics = []
    for _, grain in per_grain.groupby("grain_id"):
        grain = grain.sort_values("time")
        if len(grain) >= 3:
            metrics.append(jerkiness_metrics(grain["time"].to_numpy(), grain["area"].to_numpy()))
    jerk = float(np.mean([m["jerkiness_CV"] for m in metrics])) if metrics else np.nan
    burst = float(np.mean([m["burstiness"] for m in metrics])) if metrics else np.nan
    with (run_dir / "events.csv").open() as handle:
        events = max(sum(1 for _ in handle) - 1, 0)
    return {
        "regime": config["regime"], "temperature": config["pf"]["temperature"],
        "n": fit.exponent, "n_ci_low": np.nan, "n_ci_high": np.nan,
        "K": fit.coefficient, "K_ci": np.nan, "Q_app": np.nan, "Q_app_ci": np.nan,
        "jerkiness_CV": jerk, "Fano": np.nan, "burstiness": burst,
        "reverse_motion_fraction": np.nan, "velocity_curvature_R2": np.nan,
        "pinned_fraction": np.nan, "number_of_events": events,
        "number_of_realizations": 1, "Git_SHA": manifest["git_sha"],
    }


def analyze_campaign(campaign_dir: str | Path, output: str | Path | None = None) -> pd.DataFrame:
    campaign_dir = Path(campaign_dir)
    manifest = json.loads((campaign_dir / "campaign_manifest.json").read_text())
    raw = pd.DataFrame([analyze_run(path) for path in manifest["runs"]])
    grouped_rows = []
    for (regime, temperature), group in raw.groupby(["regime", "temperature"]):
        row = group.iloc[0].to_dict()
        row["n"] = group["n"].mean(); row["K"] = group["K"].mean()
        row["n_ci_low"] = group["n"].quantile(0.025); row["n_ci_high"] = group["n"].quantile(0.975)
        row["K_ci"] = group["K"].std(ddof=1) if len(group) > 1 else np.nan
        row["number_of_events"] = group["number_of_events"].sum()
        row["number_of_realizations"] = len(group)
        grouped_rows.append(row)
    summary = pd.DataFrame(grouped_rows, columns=SUMMARY_COLUMNS)
    for regime, indices in summary.groupby("regime").groups.items():
        subset = summary.loc[indices]
        if len(subset) >= 4 and np.all(subset["K"] > 0):
            activation = fit_activation_energy(subset["temperature"].to_numpy(), subset["K"].to_numpy())
            summary.loc[indices, "Q_app"] = activation.activation_energy_ev
            summary.loc[indices, "Q_app_ci"] = 1.96 * activation.standard_error_ev
    target = Path(output) if output else campaign_dir / "mechanism_summary.csv"
    summary.to_csv(target, index=False)
    return summary

