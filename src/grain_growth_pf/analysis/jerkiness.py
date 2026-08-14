from __future__ import annotations

import numpy as np


def jerkiness_metrics(time: np.ndarray, trajectory: np.ndarray, events: np.ndarray | None = None,
                      stationary_tolerance: float = 1e-10) -> dict[str, float]:
    time, trajectory = np.asarray(time, float), np.asarray(trajectory, float)
    velocity = np.diff(trajectory) / np.diff(time)
    magnitude = np.abs(velocity)
    mean = magnitude.mean() if len(magnitude) else 0.0
    cv = float(magnitude.std() / mean) if mean else 0.0
    stationary = float(np.mean(magnitude <= stationary_tolerance)) if len(magnitude) else 1.0
    ordered = np.sort(magnitude)[::-1]
    total = ordered.sum()
    concentration = {}
    for fraction in (0.01, 0.05, 0.10):
        count = max(1, int(np.ceil(len(ordered) * fraction)))
        concentration[f"motion_top_{int(fraction*100)}pct"] = float(ordered[:count].sum() / total) if total else 0.0
    result = {"jerkiness_CV": cv, "stationary_fraction": stationary, **concentration}
    if events is not None:
        counts = np.asarray(events, float)
        result["Fano"] = float(counts.var() / counts.mean()) if counts.mean() else 0.0
    waiting_indices = np.flatnonzero(magnitude > stationary_tolerance)
    if len(waiting_indices) > 1:
        waits = np.diff(time[1:][waiting_indices])
        result["burstiness"] = float((waits.std() - waits.mean()) / (waits.std() + waits.mean())) if waits.std() + waits.mean() else 0.0
    else:
        result["burstiness"] = -1.0
    return result

