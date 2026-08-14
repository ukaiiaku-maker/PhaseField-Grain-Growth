from __future__ import annotations

import numpy as np

from .mode import DisconnectionMode


def isotropic_surrogate_library(
    b_shells: tuple[float, ...] = (0.25, 0.5, 1.0),
    directions: int = 8,
    step_heights: tuple[float, ...] = (0.25, 0.5),
    barrier_core_ev: float = 0.35,
    b_coefficient_ev: float = 0.25,
    h_coefficient_ev: float = 0.10,
    b_power: float = 2.0,
    attempt_frequency: float = 1e9,
    site_multiplicity: float = 1.0,
    seed: int = 0,
    disorder_std_ev: float = 0.0,
) -> list[DisconnectionMode]:
    if directions < 2 or not b_shells or min(b_shells) <= 0:
        raise ValueError("finite positive Burgers shells and at least two directions are required")
    rng = np.random.default_rng(seed)
    modes: list[DisconnectionMode] = []
    b0, h0 = min(b_shells), min(map(abs, step_heights))
    for si, magnitude in enumerate(b_shells):
        family = "easy" if si == 0 else ("intermediate" if si < len(b_shells) - 1 else "high")
        for di, angle in enumerate(np.linspace(0, 2 * np.pi, directions, endpoint=False)):
            b = (magnitude * np.cos(angle), magnitude * np.sin(angle))
            for height in step_heights:
                for sign in (-1.0, 1.0):
                    h = sign * abs(height)
                    barrier = (
                        barrier_core_ev
                        + b_coefficient_ev * (magnitude / b0) ** b_power
                        + h_coefficient_ev * (abs(h) / h0) ** 2
                        + float(rng.normal(0, disorder_std_ev))
                    )
                    modes.append(DisconnectionMode(
                        f"{family}:b{si}:d{di}:h{h:g}", b, h, 0.0,
                        max(barrier, 0.0), attempt_frequency, site_multiplicity,
                        activation_volume_normal=h,
                        activation_volume_shear=magnitude,
                        delta_s=magnitude,
                        family=family,
                    ))
    return modes

