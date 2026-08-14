from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


def exponential_threshold(rng: np.random.Generator) -> float:
    # Generator.exponential avoids log(0) while being exactly Exp(1).
    return float(rng.exponential())


@dataclass
class HazardEvent:
    event_time: float
    overshoot: float
    threshold: float
    channel: int | None = None


@dataclass
class CumulativeHazardClock:
    rng: np.random.Generator
    cumulative_hazard: float = 0.0
    threshold: float | None = None
    last_rate: float | None = None

    def __post_init__(self) -> None:
        if self.threshold is None:
            self.threshold = exponential_threshold(self.rng)

    def advance(self, rate: float, dt: float, time: float,
                previous_rate: float | None = None) -> list[HazardEvent]:
        """Integrate a piecewise-linear rate and return all first passages.

        The event time is interpolated in cumulative-hazard space. Carrying an
        overshoot permits more than one event per PF step without a Bernoulli
        approximation, giving timestep-invariant integrated statistics.
        """
        if rate < 0 or dt < 0 or not np.isfinite(rate):
            raise ValueError("rate must be finite/nonnegative and dt nonnegative")
        r0 = rate if previous_rate is None and self.last_rate is None else (
            self.last_rate if previous_rate is None else previous_rate
        )
        r0 = max(float(r0), 0.0)
        increment = 0.5 * (r0 + rate) * dt
        start_h = self.cumulative_hazard
        end_h = start_h + increment
        events: list[HazardEvent] = []
        elapsed_fraction = 0.0
        while end_h >= float(self.threshold) and increment > 0:
            target = float(self.threshold) - start_h
            fraction = min(max(target / increment, elapsed_fraction), 1.0)
            event_time = time + fraction * dt
            overshoot = end_h - float(self.threshold)
            events.append(HazardEvent(event_time, overshoot, float(self.threshold)))
            start_h = float(self.threshold)
            self.threshold = start_h + exponential_threshold(self.rng)
            elapsed_fraction = fraction
        self.cumulative_hazard = end_h
        self.last_rate = float(rate)
        return events

    def reset(self) -> None:
        self.cumulative_hazard = 0.0
        self.threshold = exponential_threshold(self.rng)
        self.last_rate = None


@dataclass
class ParallelHazardClock:
    rng: np.random.Generator
    clock: CumulativeHazardClock = field(init=False)

    def __post_init__(self) -> None:
        self.clock = CumulativeHazardClock(self.rng)

    def advance(self, rates: Sequence[float], dt: float, time: float) -> list[HazardEvent]:
        rates_array = np.asarray(rates, dtype=float)
        if np.any(rates_array < 0) or np.any(~np.isfinite(rates_array)):
            raise ValueError("all channel rates must be finite and nonnegative")
        total = float(rates_array.sum())
        base = self.clock.advance(total, dt, time)
        if total == 0:
            return base
        probabilities = rates_array / total
        for event in base:
            event.channel = int(self.rng.choice(len(rates_array), p=probabilities))
        return base

