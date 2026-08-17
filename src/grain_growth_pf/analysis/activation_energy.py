from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from grain_growth_pf.disconnections.mode import K_B_EV


@dataclass(frozen=True)
class ActivationFit:
    activation_energy_ev: float
    prefactor: float
    r_squared: float
    standard_error_ev: float


def fit_activation_energy(temperatures: np.ndarray, coefficients: np.ndarray) -> ActivationFit:
    t, k = np.asarray(temperatures, float), np.asarray(coefficients, float)
    if len(t) < 4 or np.any(t <= 0) or np.any(k <= 0):
        raise ValueError("at least four positive temperatures and coefficients are required")
    x, y = 1.0 / t, np.log(k)
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    residual = y - predicted
    sst = np.sum((y - y.mean())**2)
    r2 = 1.0 - np.sum(residual**2) / sst if sst else 1.0
    slope_se = np.sqrt(np.sum(residual**2) / (len(x) - 2) / np.sum((x - x.mean())**2))
    return ActivationFit(float(-K_B_EV * slope), float(np.exp(intercept)), float(r2), float(K_B_EV * slope_se))


def local_activation_energies(temperatures: np.ndarray,
                              coefficients: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return adjacent-temperature apparent barriers and harmonic-midpoint T."""
    temperature = np.asarray(temperatures, float)
    coefficient = np.asarray(coefficients, float)
    if len(temperature) < 2 or np.any(temperature <= 0) or np.any(coefficient <= 0):
        raise ValueError("positive temperatures and coefficients are required")
    order = np.argsort(temperature)
    temperature, coefficient = temperature[order], coefficient[order]
    inverse = 1.0 / temperature
    local_q = -K_B_EV * np.diff(np.log(coefficient)) / np.diff(inverse)
    midpoint = 2.0 / (inverse[:-1] + inverse[1:])
    return midpoint, local_q
