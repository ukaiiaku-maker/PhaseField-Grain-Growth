from __future__ import annotations

from itertools import combinations_with_replacement, product
from typing import Iterable

import numpy as np

from .mode import DisconnectionMode, ModeDriving


def feasible_combinations(modes: Iterable[DisconnectionMode], target_burgers: tuple[float, float],
                          target_step: float, max_events: int = 3,
                          tolerance: float = 1e-8) -> list[tuple[DisconnectionMode, ...]]:
    candidates = list(modes)
    target_b = np.asarray(target_burgers, dtype=float)
    feasible: list[tuple[DisconnectionMode, ...]] = []
    for count in range(1, max_events + 1):
        for combo in combinations_with_replacement(candidates, count):
            b = sum((np.asarray(m.burgers) for m in combo), start=np.zeros(2))
            h = sum(m.step_height for m in combo)
            if np.linalg.norm(b - target_b) <= tolerance and abs(h - target_step) <= tolerance:
                feasible.append(combo)
    return feasible


def select_admissible_modes(modes: Iterable[DisconnectionMode], compatibility_required: bool,
                            allowed_families: set[str] | None = None) -> list[DisconnectionMode]:
    modes = list(modes)
    if allowed_families is not None:
        modes = [m for m in modes if m.family in allowed_families]
    if compatibility_required:
        secondary = [m for m in modes if m.family != "easy"]
        return secondary
    return modes


def combination_rates(combinations: Iterable[tuple[DisconnectionMode, ...]], temperature: float,
                      driving: ModeDriving) -> np.ndarray:
    # A necessary combination is serial: residence times add. Its effective
    # rate is the reciprocal total mean residence time, never a sum of rates.
    rates = []
    for combo in combinations:
        individual = np.array([m.rate(temperature, driving) for m in combo])
        rates.append(0.0 if np.any(individual == 0) else 1.0 / np.sum(1.0 / individual))
    return np.asarray(rates)

