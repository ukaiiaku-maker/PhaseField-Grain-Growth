from __future__ import annotations

from typing import Any

import numpy as np

from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.simulation import DomainPhysics, EventResolvedSimulation


class MigrationClosureSimulation(EventResolvedSimulation):
    """Event-resolved PF simulation with explicit migration/gating closures.

    ``hybrid`` preserves the legacy behavior: a successful disconnection event
    both changes compatibility state and adds its normal step to the hidden PF
    displacement ledger.

    ``gate_only`` is the corrected reduced model for compatibility-limited
    migration. A high-barrier event restores an admissible state and may relax
    shear/free-volume variables, but it does *not* add an extra normal PF
    displacement on top of the ordinary curvature-driven migration law.

    The event kinematics are still recorded in the ledger, so ``h``, ``b`` and
    strain increments remain available for diagnostics and later construction
    of a fully event-resolved continuum limit.

    A second independent control, ``blocked_gate_profile``, addresses the
    pixel-scale roughness seen in strongly climb-limited movies. ``line`` keeps
    the legacy one-pixel centerline arrest. ``diffuse`` extends the mobility
    suppression smoothly across the physical diffuse-interface width, so a
    blocked GB domain arrests the interface coherently instead of pinning a
    jagged chain of label-boundary pixels.
    """

    VALID_MIGRATION_CLOSURES = {"hybrid", "gate_only"}
    VALID_GATE_PROFILES = {"line", "diffuse"}

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        closure = str(self.config.parameters.get("migration_closure", "gate_only"))
        if closure not in self.VALID_MIGRATION_CLOSURES:
            raise ValueError(
                f"migration_closure must be one of {sorted(self.VALID_MIGRATION_CLOSURES)}, "
                f"got {closure!r}"
            )
        profile = str(self.config.parameters.get("blocked_gate_profile", "line"))
        if profile not in self.VALID_GATE_PROFILES:
            raise ValueError(
                f"blocked_gate_profile must be one of {sorted(self.VALID_GATE_PROFILES)}, "
                f"got {profile!r}"
            )
        self.migration_closure = closure
        self.blocked_gate_profile = profile

    def _record_event(
        self,
        domain: DomainPhysics,
        segment: GBSegment,
        mode: DisconnectionMode,
        rate: float,
        driving: ModeDriving,
        event_type: str,
        delta_length: float,
        event_time: float | None = None,
        *,
        ledger_position: str | None = None,
        field_position: tuple[int, int] | None = None,
        effective_barrier_ev: float | None = None,
    ) -> None:
        prior_hidden_displacement = domain.normal_displacement_ledger
        prior_release_remaining = domain.normal_release_remaining

        super()._record_event(
            domain,
            segment,
            mode,
            rate,
            driving,
            event_type,
            delta_length,
            event_time,
            ledger_position=ledger_position,
            field_position=field_position,
            effective_barrier_ev=effective_barrier_ev,
        )

        if self.migration_closure == "gate_only":
            # The event may change admissibility, Burgers residual, shear state,
            # free-volume state, and the full-field eigenstrain. It must not
            # also inject a second normal-migration channel into the PF model.
            domain.normal_displacement_ledger = prior_hidden_displacement
            domain.normal_release_remaining = prior_release_remaining

    def _apply_diffuse_blocked_gate(self) -> None:
        if self.blocked_gate_profile != "diffuse":
            return

        cfg = self.config
        dx = float(cfg.pf.grid_spacing)
        halfwidth = float(
            cfg.parameters.get("blocked_gate_halfwidth", 0.5 * cfg.pf.interface_width)
        )
        if not np.isfinite(halfwidth) or halfwidth <= 0:
            raise ValueError("blocked_gate_halfwidth must be finite and positive")
        floor = float(cfg.parameters.get("pinned_mobility_fraction", 0.0))
        if not np.isfinite(floor) or not 0.0 <= floor <= 1.0:
            raise ValueError("pinned_mobility_fraction must lie in [0,1]")

        radius = max(1, int(np.ceil(halfwidth / dx)))
        mobility = self.solver.mobility_scale.copy()
        height, width = cfg.pf.shape
        periodic = cfg.pf.boundary_conditions == "periodic"

        offsets: list[tuple[int, int, float]] = []
        for oy in range(-radius, radius + 1):
            for ox in range(-radius, radius + 1):
                distance = float(np.hypot(oy * dx, ox * dx))
                if distance > halfwidth:
                    continue
                # Raised-cosine taper: exact floor at the centerline and a
                # continuous return to the intrinsic mobility at the edge.
                weight = 0.5 * (1.0 + np.cos(np.pi * distance / halfwidth))
                target = 1.0 - (1.0 - floor) * weight
                offsets.append((oy, ox, float(target)))

        for segment in self.snapshot.boundaries.values():
            domain = self.domains.get(segment.entity_id)
            if domain is None or not domain.blocked or not len(segment.points):
                continue
            points = np.unique(segment.points.astype(int), axis=0)
            for y0, x0 in points:
                for oy, ox, target in offsets:
                    y, x = int(y0 + oy), int(x0 + ox)
                    if periodic:
                        y %= height
                        x %= width
                    elif y < 0 or y >= height or x < 0 or x >= width:
                        continue
                    mobility[y, x] = min(mobility[y, x], target)

        self.solver.set_mobility_scale(mobility)

    def _update_physics(self) -> None:
        super()._update_physics()
        self._apply_diffuse_blocked_gate()
