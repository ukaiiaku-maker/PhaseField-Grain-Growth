from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ParticleField:
    positions: np.ndarray
    radii: np.ndarray
    shape: tuple[float, float]

    @classmethod
    def random(cls, count: int, radius: float, shape: tuple[float, float], seed: int) -> "ParticleField":
        rng = np.random.default_rng(seed)
        return cls(rng.uniform([0, 0], shape, size=(count, 2)), np.full(count, radius), shape)

    def contacts(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, float)
        delta = points[:, None, :] - self.positions[None, :, :]
        box = np.asarray(self.shape)
        delta -= np.round(delta / box) * box
        return np.any(np.linalg.norm(delta, axis=-1) <= self.radii, axis=1)

    def zener_pressure_3d(self, volume_fraction: float, gb_energy: float) -> float:
        mean_radius = float(np.mean(self.radii))
        return 3.0 * volume_fraction * gb_energy / (2.0 * mean_radius)

