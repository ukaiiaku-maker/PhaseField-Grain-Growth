from __future__ import annotations

from pathlib import Path

import pandas as pd

from grain_growth_pf.analysis.campaign import SUMMARY_COLUMNS


def aggregate_summaries(summary_paths: list[str | Path], output: str | Path) -> Path:
    """Combine finalized summaries, with later inputs replacing matched rows."""
    if not summary_paths:
        raise ValueError("at least one summary is required")
    frames = []
    for priority, raw_path in enumerate(summary_paths):
        path = Path(raw_path)
        frame = pd.read_csv(path)
        missing = set(SUMMARY_COLUMNS).difference(frame.columns)
        if missing:
            raise ValueError(f"{path} lacks required columns: {sorted(missing)}")
        frame = frame[SUMMARY_COLUMNS].copy()
        frame["_priority"] = priority
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values("_priority").drop_duplicates(
        ["regime", "temperature"], keep="last"
    )
    combined = combined.drop(columns="_priority").sort_values(
        ["regime", "temperature"]
    ).reset_index(drop=True)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(target, index=False)
    return target
