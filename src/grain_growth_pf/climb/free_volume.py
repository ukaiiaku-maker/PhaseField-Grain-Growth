from __future__ import annotations

from dataclasses import dataclass


@dataclass
class FreeVolumeState:
    excess_volume_per_area: float
    point_defect_formation_volume: float
    stiffness: float
    required_total: float = 0.0
    accommodated_total: float = 0.0
    dissipated_energy: float = 0.0

    @property
    def deficit(self) -> float:
        return self.required_total - self.accommodated_total

    @property
    def energy(self) -> float:
        return 0.5 * self.stiffness * self.deficit**2

    @property
    def chemical_potential(self) -> float:
        return self.stiffness * self.deficit

    def require_for_area_change(self, area_change: float) -> float:
        if self.point_defect_formation_volume <= 0:
            raise ValueError("point-defect formation volume must be positive")
        quota = self.excess_volume_per_area * abs(area_change) / self.point_defect_formation_volume
        self.required_total += quota
        return quota

    def accommodate(self, quota: float) -> float:
        before = self.energy
        accepted = min(max(quota, 0.0), max(self.deficit, 0.0))
        self.accommodated_total += accepted
        self.dissipated_energy += max(0.0, before - self.energy)
        return accepted

    def check_balance(self, tolerance: float = 1e-12) -> bool:
        return abs(self.required_total - self.accommodated_total - self.deficit) <= tolerance
