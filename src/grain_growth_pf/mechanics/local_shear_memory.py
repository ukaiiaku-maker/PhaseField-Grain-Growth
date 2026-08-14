from __future__ import annotations

from dataclasses import dataclass


@dataclass
class LocalShearMemory:
    stiffness: float
    state: float = 0.0
    relaxation_time: float | None = None
    dissipated_energy: float = 0.0

    @property
    def energy(self) -> float:
        return 0.5 * self.stiffness * self.state**2

    @property
    def internal_shear_stress(self) -> float:
        # tau = -dE/ds under the convention ds=beta dx_n.
        return -self.stiffness * self.state

    def migrate(self, beta: float, normal_displacement: float, dt: float = 0.0) -> None:
        old = self.energy
        self.state += beta * normal_displacement
        if self.relaxation_time is not None and self.relaxation_time > 0 and dt > 0:
            before = self.energy
            self.state *= max(0.0, 1.0 - dt / self.relaxation_time)
            self.dissipated_energy += max(0.0, before - self.energy)
        # Stored work may rise during forced migration; this is supplied by the
        # capillary/mechanical driving and tracked rather than called dissipation.
        _ = old

    def release(self, delta_s: float) -> float:
        before = self.energy
        direction = 1.0 if self.state >= 0 else -1.0
        self.state -= direction * min(abs(delta_s), abs(self.state))
        released = max(0.0, before - self.energy)
        self.dissipated_energy += released
        return released

    def normal_velocity(self, mobility: float, capillary_pressure: float,
                        beta: float, chemical_pressure: float = 0.0) -> float:
        return mobility * (capillary_pressure + beta * self.internal_shear_stress + chemical_pressure)

