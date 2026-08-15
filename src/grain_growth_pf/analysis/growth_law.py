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


@dataclass(frozen=True)
class GrowthProfile:
    exponents: np.ndarray
    normalized_rmse: np.ndarray
    residual_autocorrelation: np.ndarray


@dataclass(frozen=True)
class CommonExponentFit:
    exponent: float
    coefficients: np.ndarray
    initial_radii: np.ndarray
    normalized_rmse: float


def scan_growth_exponent(time: np.ndarray, radius: np.ndarray,
                         exponents: np.ndarray | None = None) -> GrowthProfile:
    """Profile the generalized growth law in the measured-radius error space."""
    time, radius = np.asarray(time, float), np.asarray(radius, float)
    valid = np.isfinite(time) & np.isfinite(radius) & (radius > 0)
    time, radius = time[valid], radius[valid]
    if len(time) < 3:
        raise ValueError("at least three samples are required")
    grid = np.asarray(exponents if exponents is not None else np.linspace(1.0, 6.0, 251), float)
    errors, autocorrelations = [], []
    shifted_time = time - time[0]
    radius_scale = max(float(np.std(radius)), 1e-15)
    for exponent in grid:
        coefficient = max(float(np.polyfit(shifted_time, radius**exponent, 1)[0]), 1e-15)

        def fixed_exponent_residual(parameters: np.ndarray) -> np.ndarray:
            prediction = np.maximum(
                parameters[1] ** exponent + parameters[0] * shifted_time, 1e-30
            ) ** (1.0 / exponent)
            return (prediction - radius) / radius_scale

        result = least_squares(
            fixed_exponent_residual,
            x0=np.array([coefficient, radius[0]]),
            bounds=([0.0, 1e-15], [np.inf, np.inf]),
            x_scale="jac",
        )
        residual = fixed_exponent_residual(result.x) * radius_scale
        errors.append(float(np.sqrt(np.mean(residual**2)) / radius_scale))
        autocorrelations.append(
            float(np.corrcoef(residual[:-1], residual[1:])[0, 1])
            if len(residual) > 2 and np.std(residual) else 0.0
        )
    return GrowthProfile(grid, np.asarray(errors), np.asarray(autocorrelations))


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

    shifted_time = time - time[0]
    radius_scale = max(float(np.std(radius)), 1e-15)

    def objective(parameters: np.ndarray) -> np.ndarray:
        exponent, coefficient, initial_radius = parameters
        prediction = np.maximum(
            initial_radius**exponent + coefficient * shifted_time, 1e-30
        ) ** (1.0 / exponent)
        return (prediction - radius) / radius_scale

    candidates = []
    for initial_exponent in (1.01, 2.0, 3.0, 5.5):
        initial_exponent = float(np.clip(initial_exponent, *n_bounds))
        initial_coefficient = max(
            float(np.polyfit(shifted_time, radius**initial_exponent, 1)[0]), 1e-15
        )
        candidates.append(least_squares(
            objective,
            x0=np.array([initial_exponent, initial_coefficient, radius[0]]),
            bounds=([n_bounds[0], 0.0, 1e-15], [n_bounds[1], np.inf, np.inf]),
            x_scale="jac",
            max_nfev=10_000,
        ))
    result = min(candidates, key=lambda candidate: candidate.cost)
    exponent, coefficient, initial_radius = map(float, result.x)
    prediction = np.maximum(
        initial_radius**exponent + coefficient * shifted_time, 1e-30
    ) ** (1.0 / exponent)
    residual = radius - prediction
    intercept = initial_radius**exponent - coefficient * time[0]
    sst = float(np.sum((radius - radius.mean())**2))
    r2 = 1.0 - float(np.sum(residual**2)) / sst if sst else 1.0
    ac = float(np.corrcoef(residual[:-1], residual[1:])[0, 1]) if len(residual) > 2 and np.std(residual) else 0.0
    return GrowthFit(exponent, coefficient, intercept, r2, ac, float(time[0]), float(time[-1]))


def fit_growth_law_fixed_exponent(time: np.ndarray, radius: np.ndarray, exponent: float,
                                  transient_fraction: float = 0.2) -> GrowthFit:
    """Fit ``K`` and ``R0`` in radius space for a prescribed common exponent."""
    if exponent <= 0:
        raise ValueError("exponent must be positive")
    time, radius = np.asarray(time, float), np.asarray(radius, float)
    valid = np.isfinite(time) & np.isfinite(radius) & (radius > 0)
    time, radius = time[valid], radius[valid]
    start = min(max(int(len(time) * transient_fraction), 0), max(len(time) - 3, 0))
    time, radius = time[start:], radius[start:]
    if len(time) < 3:
        raise ValueError("at least three post-transient samples are required")
    shifted_time = time - time[0]
    radius_scale = max(float(np.std(radius)), 1e-15)
    initial_coefficient = max(
        float(np.polyfit(shifted_time, radius**exponent, 1)[0]), 1e-15
    )

    def objective(parameters: np.ndarray) -> np.ndarray:
        prediction = np.maximum(
            parameters[1] ** exponent + parameters[0] * shifted_time, 1e-30
        ) ** (1.0 / exponent)
        return (prediction - radius) / radius_scale

    result = least_squares(
        objective,
        x0=np.array([initial_coefficient, radius[0]]),
        bounds=([0.0, 1e-15], [np.inf, np.inf]),
        x_scale="jac",
    )
    coefficient, initial_radius = map(float, result.x)
    prediction = np.maximum(
        initial_radius**exponent + coefficient * shifted_time, 1e-30
    ) ** (1.0 / exponent)
    residual = radius - prediction
    intercept = initial_radius**exponent - coefficient * time[0]
    sst = float(np.sum((radius - radius.mean())**2))
    r2 = 1.0 - float(np.sum(residual**2)) / sst if sst else 1.0
    ac = float(np.corrcoef(residual[:-1], residual[1:])[0, 1]) if len(residual) > 2 and np.std(residual) else 0.0
    return GrowthFit(float(exponent), coefficient, intercept, r2, ac,
                     float(time[0]), float(time[-1]))


def fit_common_exponent(time_series: list[np.ndarray], radius_series: list[np.ndarray],
                        n_bounds: tuple[float, float] = (1.0, 6.0)) -> CommonExponentFit:
    """Jointly fit one exponent and a separate coefficient at each condition."""
    if len(time_series) != len(radius_series) or not time_series:
        raise ValueError("time and radius series must be nonempty and paired")
    times, radii, scales = [], [], []
    for raw_time, raw_radius in zip(time_series, radius_series):
        time = np.asarray(raw_time, float)
        radius = np.asarray(raw_radius, float)
        valid = np.isfinite(time) & np.isfinite(radius) & (radius > 0)
        time, radius = time[valid], radius[valid]
        if len(time) < 3:
            raise ValueError("each condition requires at least three samples")
        times.append(time - time[0])
        radii.append(radius)
        scales.append(max(float(np.std(radius)), 1e-15) * np.sqrt(len(radius)))

    def objective(parameters: np.ndarray) -> np.ndarray:
        exponent = parameters[0]
        coefficients = np.exp(parameters[1:1 + len(radii)])
        initial_radii = np.exp(parameters[1 + len(radii):])
        residuals = []
        for time, radius, scale, coefficient, initial_radius in zip(
            times, radii, scales, coefficients, initial_radii
        ):
            prediction = np.maximum(
                initial_radius**exponent + coefficient * time, 1e-30
            ) ** (1.0 / exponent)
            residuals.append((prediction - radius) / scale)
        return np.concatenate(residuals)

    candidates = []
    for initial_exponent in (1.01, 2.0, 3.0, 5.5):
        initial_exponent = float(np.clip(initial_exponent, *n_bounds))
        initial_coefficients = [
            max(float(np.polyfit(time, radius**initial_exponent, 1)[0]), 1e-15)
            for time, radius in zip(times, radii)
        ]
        x0 = np.array([
            initial_exponent,
            *np.log(initial_coefficients),
            *np.log([radius[0] for radius in radii]),
        ])
        lower = np.array([n_bounds[0], *([-50.0] * (2 * len(radii)))])
        upper = np.array([n_bounds[1], *([50.0] * (2 * len(radii)))])
        candidates.append(least_squares(
            objective, x0=x0, bounds=(lower, upper), x_scale="jac", max_nfev=10_000
        ))
    result = min(candidates, key=lambda candidate: candidate.cost)
    exponent = float(result.x[0])
    coefficients = np.exp(result.x[1:1 + len(radii)])
    initial_radii = np.exp(result.x[1 + len(radii):])
    normalized_rmse = float(np.sqrt(np.mean(objective(result.x) ** 2)))
    return CommonExponentFit(exponent, coefficients, initial_radii, normalized_rmse)


def bootstrap_exponent(time: np.ndarray, radii_by_realization: np.ndarray, samples: int,
                       seed: int, transient_fraction: float = 0.2) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    radii = np.asarray(radii_by_realization, float)
    estimates = []
    for _ in range(samples):
        selection = rng.integers(0, len(radii), len(radii))
        estimates.append(fit_growth_law(
            time, radii[selection].mean(axis=0), transient_fraction=transient_fraction
        ).exponent)
    return tuple(np.quantile(estimates, [0.025, 0.975]))
