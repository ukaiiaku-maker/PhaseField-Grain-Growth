from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
from scipy.special import gammainc, lambertw


@dataclass(frozen=True)
class CrossoverGrowthFit:
    model: str
    intrinsic_constant: float
    crossover_strength: float
    size_exponent: float
    r_squared: float
    normalized_rmse: float
    parameter_at_bound: bool


def crossover_radius_prediction(time: np.ndarray, initial_radius: float,
                                intrinsic_constant: float,
                                crossover_strength: float,
                                size_exponent: float) -> np.ndarray:
    """Invert the exact Class-B additive growth-time expression.

    K dt = (R^2-R0^2)/2 + c (R^(p+2)-R0^(p+2))/(p+2).
    The Class-C exchange crossover is the special case p=1 and c=1/Rx.
    """
    elapsed = np.asarray(time, dtype=float) - float(np.asarray(time, dtype=float)[0])
    if intrinsic_constant <= 0 or crossover_strength < 0 or size_exponent <= 0:
        raise ValueError("crossover growth parameters must be positive")
    target = intrinsic_constant * elapsed
    lower = np.full_like(target, initial_radius)
    upper = np.full_like(target, max(initial_radius * 1.1, initial_radius + 1.0))

    def integrated(radius: np.ndarray) -> np.ndarray:
        return (
            0.5 * (radius**2 - initial_radius**2)
            + crossover_strength / (size_exponent + 2.0)
            * (radius ** (size_exponent + 2.0) - initial_radius ** (size_exponent + 2.0))
        )

    for _ in range(64):
        needs_growth = integrated(upper) < target
        if not np.any(needs_growth):
            break
        upper[needs_growth] *= 2.0
    for _ in range(64):
        middle = 0.5 * (lower + upper)
        below = integrated(middle) < target
        lower[below] = middle[below]
        upper[~below] = middle[~below]
    return 0.5 * (lower + upper)


def fit_crossover_growth(time: np.ndarray, radius: np.ndarray,
                         size_exponent: float | None = None) -> CrossoverGrowthFit:
    """Fit manuscript Class-B/Class-C kinetics directly in radius space."""
    time = np.asarray(time, dtype=float)
    radius = np.asarray(radius, dtype=float)
    if len(time) < 4 or len(time) != len(radius) or np.any(np.diff(time) <= 0):
        raise ValueError("fit requires at least four radius samples at increasing times")
    r0 = float(radius[0])
    scale = max(float(np.ptp(radius)), np.finfo(float).eps)
    elapsed = float(time[-1] - time[0])
    initial_k = max(float(radius[-1] ** 2 - r0**2) / (2.0 * elapsed), 1e-8)
    log_k_bounds = (np.log(1e-10), np.log(max(1e6, initial_k * 1e4)))
    log_c_bounds = (np.log(1e-12), np.log(1e4))

    def residual(parameters: np.ndarray) -> np.ndarray:
        exponent = float(size_exponent if size_exponent is not None else parameters[2])
        predicted = crossover_radius_prediction(
            time, r0, np.exp(parameters[0]), np.exp(parameters[1]), exponent
        )
        return (predicted - radius) / scale

    starts = (0.5, 1.0, 2.0, 3.0) if size_exponent is None else (size_exponent,)
    fits = []
    for exponent_start in starts:
        x0 = [np.log(initial_k), np.log(max(r0 ** (-exponent_start), 1e-10))]
        lower = [log_k_bounds[0], log_c_bounds[0]]
        upper = [log_k_bounds[1], log_c_bounds[1]]
        if size_exponent is None:
            x0.append(exponent_start)
            lower.append(0.1)
            upper.append(4.0)
        fits.append(least_squares(residual, x0, bounds=(lower, upper), max_nfev=1000))
    fit = min(fits, key=lambda item: np.sum(item.fun**2))
    exponent = float(size_exponent if size_exponent is not None else fit.x[2])
    predicted = crossover_radius_prediction(
        time, r0, np.exp(fit.x[0]), np.exp(fit.x[1]), exponent
    )
    residual_sum = float(np.sum((radius - predicted) ** 2))
    total_sum = float(np.sum((radius - radius.mean()) ** 2))
    at_bound = bool(
        np.isclose(fit.x[0], log_k_bounds, atol=1e-4).any()
        or np.isclose(fit.x[1], log_c_bounds, atol=1e-4).any()
        or (size_exponent is None and (fit.x[2] <= 0.1001 or fit.x[2] >= 3.9999))
    )
    return CrossoverGrowthFit(
        model="class_c_exchange" if size_exponent == 1.0 else "class_b_additive",
        intrinsic_constant=float(np.exp(fit.x[0])),
        crossover_strength=float(np.exp(fit.x[1])),
        size_exponent=exponent,
        r_squared=1.0 - residual_sum / total_sum if total_sum else 1.0,
        normalized_rmse=float(np.sqrt(np.mean((radius - predicted) ** 2)) / scale),
        parameter_at_bound=at_bound,
    )


def intrinsic_radius(time: np.ndarray, initial_radius: float, growth_constant: float) -> np.ndarray:
    """R^2=R0^2+2 K t for dR/dt=K/R."""
    return np.sqrt(np.maximum(initial_radius**2 + 2.0 * growth_constant * np.asarray(time), 0.0))


def poisson_activity(required_hits: int, hazard: np.ndarray | float) -> np.ndarray:
    return gammainc(required_hits, np.asarray(hazard, dtype=float))


def power_hazard(radius: np.ndarray | float, amplitude: float, exponent: float) -> np.ndarray:
    return amplitude * np.asarray(radius, dtype=float) ** (-exponent)


def exchange_activity(radius: np.ndarray | float, crossover_radius: float) -> np.ndarray:
    return 1.0 / (1.0 + np.asarray(radius, dtype=float) / crossover_radius)


def drag_activity(drag_number: np.ndarray | float) -> np.ndarray:
    return 1.0 / (1.0 + np.asarray(drag_number, dtype=float))


def series_activity(*activities: np.ndarray | float) -> np.ndarray:
    values = np.asarray(activities, dtype=float)
    return 1.0 / np.sum(1.0 / np.maximum(values, np.finfo(float).tiny), axis=0)


def parallel_activity(*activities: np.ndarray | float) -> np.ndarray:
    values = np.clip(np.asarray(activities, dtype=float), 0.0, 1.0)
    return 1.0 - np.prod(1.0 - values, axis=0)


def work_limited_activity(required_hits: int, raw_hazard: np.ndarray | float,
                          force_gate: np.ndarray | float = 1.0,
                          work_gate: np.ndarray | float = 1.0) -> np.ndarray:
    return np.asarray(force_gate) * np.asarray(work_gate) * poisson_activity(required_hits, raw_hazard)


def integrate_growth_law(radius_grid: np.ndarray, intrinsic_constant: float,
                         activity: np.ndarray) -> np.ndarray:
    radius = np.asarray(radius_grid, float)
    gamma = np.asarray(activity, float)
    dt_dR = radius / (intrinsic_constant * np.maximum(gamma, np.finfo(float).tiny))
    return np.concatenate(([0.0], cumulative_trapezoid(dt_dR, radius)))


def lambert_hazard_solution(z: np.ndarray | float) -> np.ndarray:
    """Positive real solution of Y exp(Y)=Z."""
    return lambertw(np.asarray(z, dtype=float)).real


def asymptotic_exponent(hazard_size_exponent: float, hits: int = 1) -> float:
    """Rare-event Class-B closure: Gamma~Lambda^K gives n=2+pK."""
    return 2.0 + hazard_size_exponent * hits
