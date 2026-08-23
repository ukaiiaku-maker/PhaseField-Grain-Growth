from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class InverseGrowthRateFit:
    """Linear diagnostic fit for ``1 / Rdot = a R + b``."""

    a: float
    b: float
    r_squared: float


def fit_inverse_growth_rate(
    radius: np.ndarray,
    growth_rate: np.ndarray,
) -> InverseGrowthRateFit:
    """Fit the explicitly labelled inverse-growth-rate resistance law."""
    x = np.asarray(radius, dtype=float)
    rate = np.asarray(growth_rate, dtype=float)
    if x.shape != rate.shape or x.ndim != 1 or len(x) < 3:
        raise ValueError("radius and growth_rate must be matching 1-D arrays of length >= 3")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(rate)) or np.any(rate <= 0.0):
        raise ValueError("inverse-growth-rate fitting requires finite positive growth rates")
    inverse_rate = 1.0 / rate
    slope, intercept = np.polyfit(x, inverse_rate, 1)
    prediction = slope * x + intercept
    total = float(np.sum((inverse_rate - inverse_rate.mean()) ** 2))
    residual = float(np.sum((inverse_rate - prediction) ** 2))
    r_squared = 1.0 - residual / total if total else np.nan
    return InverseGrowthRateFit(float(slope), float(intercept), float(r_squared))
