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

    def advance(self, delta_measure: float,
                maximum_events: int | None = None) -> list[GeometricEncounter]:
        if delta_measure < 0:
            raise ValueError("pass the magnitude of physical geometric change")
        if maximum_events is not None and maximum_events < 1:
            raise ValueError("maximum_events must be positive when specified")
        start_measure = self.total_measure
        start_hazard = self.cumulative_hazard
        end = self.cumulative_hazard + self.density * delta_measure
        result: list[GeometricEncounter] = []
        while end >= self.threshold and self.density > 0:
            event_threshold = self.threshold
            q_event = start_measure + (event_threshold - start_hazard) / self.density
            result.append(GeometricEncounter(q_event, event_threshold, end - event_threshold))
            self.threshold += float(self.rng.exponential())
            if maximum_events is not None and len(result) >= maximum_events:
                # Encounter changes the physical mobility state. Discard the
                # geometry remainder that would occur only after that change.
                result[-1].overshoot = 0.0
                self.cumulative_hazard = event_threshold
                self.total_measure = q_event
                return result
        self.cumulative_hazard = end
        self.total_measure = start_measure + delta_measure
        return result
