from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import savgol_filter


@dataclass(frozen=True)
class GrainSizeMoments:
    number_mean: float
    population: float
    area_weighted: float


@dataclass(frozen=True)
class WindowGrowthFit:
    n_best: float
    coefficient: float
    normalized_rmse: float
    profile_width: float
    hit_search_bound: bool
    k2: float
    k3: float
    mean_growth_rate: float
    normalized_growth_rate: float
    exponents: np.ndarray
    residual_profile: np.ndarray


def grain_size_moments(areas: np.ndarray) -> GrainSizeMoments:
    values = np.asarray(areas, dtype=float)
    values = values[np.isfinite(values) & (values > 0.0)]
    if not len(values):
        raise ValueError("at least one positive finite grain area is required")
    diameters = 2.0 * np.sqrt(values / np.pi)
    return GrainSizeMoments(
        number_mean=float(diameters.mean()),
        population=float(2.0 * np.sqrt(values.sum() / (np.pi * len(values)))),
        area_weighted=float(np.average(diameters, weights=values)),
    )


def _through_origin_slope(x: np.ndarray, y: np.ndarray) -> float:
    denominator = float(x @ x)
    return max(float(x @ y) / denominator, 0.0) if denominator > 0.0 else 0.0


def profile_growth_window(
    time: np.ndarray,
    size: np.ndarray,
    *,
    n_min: float = 1.0,
    n_max: float = 20.0,
    n_step: float = 0.1,
    extended_n_max: float = 50.0,
) -> WindowGrowthFit:
    """Fit ``G**n-Ga**n=K(t-ta)`` in measured-size error space."""
    time = np.asarray(time, dtype=float)
    size = np.asarray(size, dtype=float)
    valid = np.isfinite(time) & np.isfinite(size) & (size > 0.0)
    time, size = time[valid], size[valid]
    if len(time) < 5 or np.any(np.diff(time) <= 0.0):
        raise ValueError("a growth window requires five samples at increasing times")
    elapsed = time - time[0]
    scale = max(float(np.ptp(size)), 1e-15)

    def scan(upper: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        exponents = np.arange(n_min, upper + 0.5 * n_step, n_step)
        coefficients = np.empty_like(exponents)
        residuals = np.empty_like(exponents)
        for index, exponent in enumerate(exponents):
            transformed = size**exponent - size[0]**exponent
            coefficient = _through_origin_slope(elapsed, transformed)
            prediction = np.maximum(
                size[0]**exponent + coefficient * elapsed, 1e-300
            ) ** (1.0 / exponent)
            coefficients[index] = coefficient
            residuals[index] = np.sqrt(np.mean((prediction - size) ** 2)) / scale
        return exponents, coefficients, residuals

    exponents, coefficients, residuals = scan(n_max)
    best = int(np.argmin(residuals))
    if best == len(exponents) - 1 and extended_n_max > n_max:
        exponents, coefficients, residuals = scan(extended_n_max)
        best = int(np.argmin(residuals))
    minimum = float(residuals[best])
    tolerance = minimum * 1.05 + 1e-12
    supported = exponents[residuals <= tolerance]
    width = float(supported[-1] - supported[0]) if len(supported) else 0.0
    k2 = _through_origin_slope(elapsed, size**2 - size[0]**2)
    k3 = _through_origin_slope(elapsed, size**3 - size[0]**3)
    mean_rate = float((size[-1] - size[0]) / elapsed[-1])
    return WindowGrowthFit(
        n_best=float(exponents[best]), coefficient=float(coefficients[best]),
        normalized_rmse=minimum, profile_width=width,
        hit_search_bound=best == len(exponents) - 1,
        k2=k2, k3=k3, mean_growth_rate=mean_rate,
        normalized_growth_rate=float(mean_rate / size[0]),
        exponents=exponents, residual_profile=residuals,
    )


def differential_effective_exponent(
    time: np.ndarray, size: np.ndarray, *, window_length: int = 11,
) -> tuple[np.ndarray, np.ndarray]:
    """Return smoothed ``Gdot`` and ``1-d ln(Gdot)/d ln(G)`` cross-checks."""
    time = np.asarray(time, dtype=float)
    size = np.asarray(size, dtype=float)
    if len(time) < 5 or len(time) != len(size):
        raise ValueError("at least five paired samples are required")
    window = min(window_length, len(time) if len(time) % 2 else len(time) - 1)
    window = max(window, 5)
    smooth = savgol_filter(size, window, min(3, window - 2), mode="interp")
    rate = np.gradient(smooth, time)
    positive = (rate > 0.0) & (smooth > 0.0)
    exponent = np.full(len(rate), np.nan)
    if np.count_nonzero(positive) >= 3:
        indices = np.flatnonzero(positive)
        # Digitized grain populations can give identical consecutive G values.
        # A derivative with respect to log(G) is undefined on such plateaus, so
        # retain one sample per strictly distinct smoothed size and leave the
        # other locations explicitly unavailable instead of emitting infinities.
        _, unique_offsets = np.unique(smooth[indices], return_index=True)
        distinct = indices[np.sort(unique_offsets)]
        if len(distinct) >= 3:
            exponent[distinct] = 1.0 - np.gradient(
                np.log(rate[distinct]), np.log(smooth[distinct])
            )
    return rate, exponent
