from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammainc

from .hazard import CumulativeHazardClock, HazardEvent


def poisson_completion_probability(required_hits: int, window_hazard: float) -> float:
    if required_hits < 1 or window_hazard < 0:
        raise ValueError("K >= 1 and Lambda >= 0 are required")
    # P[N >= K] = gammainc(K, Lambda); stable at both tails.
    return float(gammainc(required_hits, window_hazard))


@dataclass
class CompletionEvent:
    time: float
    hits: int


class MultiHitProcess:
    def __init__(self, required_hits: int, rng: np.random.Generator,
                 interpretation: str = "persistent_hits"):
        if required_hits < 1:
            raise ValueError("required_hits must be at least one")
        if interpretation not in {"persistent_hits", "packet_reset"}:
            raise ValueError("unknown multihit interpretation")
        self.required_hits = required_hits
        self.interpretation = interpretation
        self.hit_count = 0
        self.clock = CumulativeHazardClock(rng)

    def advance(self, rate: float, dt: float, time: float) -> list[CompletionEvent]:
        completions: list[CompletionEvent] = []
        for event in self.clock.advance(rate, dt, time):
            self.hit_count += 1
            if self.hit_count >= self.required_hits:
                completions.append(CompletionEvent(event.event_time, self.hit_count))
                self.hit_count = 0
        return completions

    def close_window(self) -> bool:
        complete = self.hit_count >= self.required_hits
        if self.interpretation == "packet_reset":
            self.hit_count = 0
            self.clock.reset()
        return complete

