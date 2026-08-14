from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import gammainc, lambertw


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

