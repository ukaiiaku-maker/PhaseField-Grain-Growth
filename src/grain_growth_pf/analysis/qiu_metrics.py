from __future__ import annotations

import numpy as np


def correlation_metrics(velocity: np.ndarray, curvature: np.ndarray,
                        shear: np.ndarray | None = None, deficit: np.ndarray | None = None) -> dict[str, float]:
    def corr2(a: np.ndarray, b: np.ndarray) -> float:
        c = np.corrcoef(a, b)[0, 1]
        return float(c * c) if np.isfinite(c) else 0.0
    result = {"velocity_curvature_R2": corr2(np.asarray(velocity), np.asarray(curvature))}
    result["reverse_motion_fraction"] = float(np.mean(np.asarray(velocity) * np.asarray(curvature) < 0))
    if shear is not None:
        result["velocity_shear_R2"] = corr2(np.asarray(velocity), np.asarray(shear))
    if deficit is not None:
        result["velocity_deficit_R2"] = corr2(np.asarray(velocity), np.asarray(deficit))
    return result

