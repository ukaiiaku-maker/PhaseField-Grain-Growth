from __future__ import annotations

from typing import Any

from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.simulation import DomainPhysics, EventResolvedSimulation


class MigrationClosureSimulation(EventResolvedSimulation):
    """Event-resolved PF simulation with an explicit migration-closure choice.

    ``hybrid`` preserves the legacy behavior: a successful disconnection event
    both changes compatibility state and adds its normal step to the hidden PF
    displacement ledger.

    ``gate_only`` is the corrected reduced model for compatibility-limited
    migration.  A high-barrier event restores an admissible state and may relax
    shear/free-volume variables, but it does *not* add an extra normal PF
    displacement on top of the ordinary curvature-driven migration law.

    The event kinematics are still recorded in the ledger, so ``h``, ``b`` and
    strain increments remain available for diagnostics and later construction
    of a fully event-resolved continuum limit.
    """

    VALID_MIGRATION_CLOSURES = {"hybrid", "gate_only"}

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        closure = str(self.config.parameters.get("migration_closure", "gate_only"))
        if closure not in self.VALID_MIGRATION_CLOSURES:
            raise ValueError(
                f"migration_closure must be one of {sorted(self.VALID_MIGRATION_CLOSURES)}, "
                f"got {closure!r}"
            )
        self.migration_closure = closure

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
            # free-volume state, and the full-field eigenstrain.  It must not
            # also inject a second normal-migration channel into the PF model.
            domain.normal_displacement_ledger = prior_hidden_displacement
            domain.normal_release_remaining = prior_release_remaining
