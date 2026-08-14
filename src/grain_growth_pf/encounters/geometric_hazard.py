from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class GeometricEncounter:
    measure: float
    threshold: float
    overshoot: float


class GeometricEncounterClock:
    """Poisson process in a monotone physical reaction coordinate Q."""

    def __init__(self, density: float, rng: np.random.Generator):
        if density < 0:
            raise ValueError("encounter density must be nonnegative")
        self.density = density
        self.rng = rng
        self.cumulative_hazard = 0.0
        self.threshold = float(rng.exponential())
        self.total_measure = 0.0

    def advance(self, delta_measure: float) -> list[GeometricEncounter]:
        if delta_measure < 0:
            raise ValueError("pass the magnitude of physical geometric change")
        self.total_measure += delta_measure
        end = self.cumulative_hazard + self.density * delta_measure
        result: list[GeometricEncounter] = []
        while end >= self.threshold and self.density > 0:
            q_event = self.total_measure - (end - self.threshold) / self.density
            result.append(GeometricEncounter(q_event, self.threshold, end - self.threshold))
            self.threshold += float(self.rng.exponential())
        self.cumulative_hazard = end
        return result

