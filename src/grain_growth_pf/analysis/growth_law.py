from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


@dataclass(frozen=True)
class GrowthFit:
    exponent: float
    coefficient: float
    intercept: float
    r_squared: float
    residual_autocorrelation: float
    fit_start: float
    fit_end: float


def _linear_for_n(time: np.ndarray, radius: np.ndarray, exponent: float) -> tuple[float, float, np.ndarray]:
    y = radius**exponent
    design = np.column_stack((time, np.ones_like(time)))
    coefficient, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    residual = y - (coefficient * time + intercept)
    return float(coefficient), float(intercept), residual


def fit_growth_law(time: np.ndarray, radius: np.ndarray, n_bounds: tuple[float, float] = (1.0, 6.0),
                   transient_fraction: float = 0.2) -> GrowthFit:
    time, radius = np.asarray(time, float), np.asarray(radius, float)
    valid = np.isfinite(time) & np.isfinite(radius) & (radius > 0)
    time, radius = time[valid], radius[valid]
    start = min(max(int(len(time) * transient_fraction), 0), max(len(time) - 3, 0))
    time, radius = time[start:], radius[start:]
    if len(time) < 3:
        raise ValueError("at least three post-transient samples are required")

    def objective(x: np.ndarray) -> np.ndarray:
        transformed = radius ** float(x[0])
        # Normalize by signal variation, not its mean. Mean normalization
        # spuriously rewards n->1 whenever R is large compared with Delta R.
        return _linear_for_n(time, radius, float(x[0]))[2] / max(np.std(transformed), 1e-15)

    result = least_squares(objective, x0=np.array([2.0]), bounds=n_bounds)
    exponent = float(result.x[0])
    coefficient, intercept, residual = _linear_for_n(time, radius, exponent)
    y = radius**exponent
    sst = float(np.sum((y - y.mean())**2))
    r2 = 1.0 - float(np.sum(residual**2)) / sst if sst else 1.0
    ac = float(np.corrcoef(residual[:-1], residual[1:])[0, 1]) if len(residual) > 2 and np.std(residual) else 0.0
    return GrowthFit(exponent, coefficient, intercept, r2, ac, float(time[0]), float(time[-1]))


def bootstrap_exponent(time: np.ndarray, radii_by_realization: np.ndarray, samples: int,
                       seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    radii = np.asarray(radii_by_realization, float)
    estimates = []
    for _ in range(samples):
        selection = rng.integers(0, len(radii), len(radii))
        estimates.append(fit_growth_law(time, radii[selection].mean(axis=0)).exponent)
    return tuple(np.quantile(estimates, [0.025, 0.975]))
