from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_tracks(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)


def ensemble_radius(tracks: pd.DataFrame) -> pd.DataFrame:
    grouped = tracks.groupby(["run_id", "time", "step"], as_index=False).agg(
        mean_area=("area", "mean"), grain_count=("grain_id", "nunique")
    )
    grouped["R_A"] = np.sqrt(grouped["mean_area"] / np.pi)
    return grouped

