from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class GBSegment:
    grain_i: int
    grain_j: int
    segment_id: int
    length: float = 0.0
    curvature: float = 0.0
    normal: tuple[float, float] = (0.0, 0.0)
    velocity: float = 0.0
    misorientation: float = 0.0
    inclination: float = 0.0
    current_mode: str | None = None
    shear_incompatibility: float = 0.0
    free_volume_deficit: float = 0.0
    encounter_state: dict[str, object] = field(default_factory=dict)
    activation_state: dict[str, object] = field(default_factory=dict)
    remaining_release_quota: float = 0.0
    points: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))
    age: int = 0

    @property
    def entity_id(self) -> str:
        return f"gb:{min(self.grain_i,self.grain_j)}-{max(self.grain_i,self.grain_j)}:{self.segment_id}"

