from __future__ import annotations

import numpy as np

from grain_growth_pf.stochastic.multihit import MultiHitProcess


class EntityPin:
    def __init__(self, entity_id: str, required_hits: int, rng: np.random.Generator,
                 interpretation: str = "persistent_hits"):
        self.entity_id = entity_id
        self.pinned = True
        self.release = MultiHitProcess(required_hits, rng, interpretation)

    def advance(self, rate: float, dt: float, time: float) -> bool:
        if self.pinned and self.release.advance(
            rate, dt, time, stop_after_completion=True
        ):
            self.pinned = False
        return not self.pinned
