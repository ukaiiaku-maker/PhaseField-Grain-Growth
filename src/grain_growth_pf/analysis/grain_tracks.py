from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_tracks(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def ensemble_radius(tracks: pd.DataFrame) -> pd.DataFrame:
    enriched = tracks.assign(radius_squared=tracks["radius"].to_numpy(float) ** 2)
    # A resumed run receives a new provenance run_id but remains one physical
    # realization. Group on physical clocks so copied checkpoint trajectories
    # are not split into lexicographically ordered run-id blocks.
    grouped = enriched.groupby(["time", "step"], as_index=False).agg(
        run_id=("run_id", "last"),
        mean_area=("area", "mean"),
        mean_radius=("radius", "mean"),
        median_radius=("radius", "median"),
        mean_radius_squared=("radius_squared", "mean"),
        mean_perimeter=("perimeter", "mean"),
        grain_count=("grain_id", "nunique"),
    )
    grouped["R_A"] = np.sqrt(grouped["mean_area"] / np.pi)
    grouped["R_mean"] = grouped["mean_radius"]
    grouped["R_median"] = grouped["median_radius"]
    grouped["R_rms"] = np.sqrt(grouped["mean_radius_squared"])
    grouped["R_perimeter"] = grouped["mean_perimeter"] / (2.0 * np.pi)
    return grouped.sort_values(["time", "step"]).reset_index(drop=True)
