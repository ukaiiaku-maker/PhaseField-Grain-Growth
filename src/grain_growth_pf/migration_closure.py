from __future__ import annotations

import csv
import json
from dataclasses import replace
from typing import Any

import numpy as np

from grain_growth_pf.disconnections.barriers import assign_barriers
from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving, K_B_EV
from grain_growth_pf.entities.arclength_tracker import ArclengthEntityTracker
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.entities.triple_junction import TripleJunction
from grain_growth_pf.pf.kinematics import interface_kinematics
from grain_growth_pf.simulation import DomainPhysics, EventResolvedSimulation


class MigrationClosureSimulation(EventResolvedSimulation):
    """Event-resolved PF simulation with corrected jerky-growth closures.

    ``hybrid`` preserves the legacy behavior in which successful events add a
    second normal PF migration channel. ``gate_only`` is the corrected reduced
    model: events change kinetic admissibility and internal state, while actual
    GB migration remains capillary/mechanical PF motion.

    The corrected branch also supports connected physical-arclength GB domains,
    local swept-area encounter measures, local shear accumulation, a physical TJ
    gate length, and an explicit activation-work diagnostic ledger.
    """

    VALID_MIGRATION_CLOSURES = {"hybrid", "gate_only"}
    VALID_GATE_PROFILES = {"line", "diffuse"}
    ACTIVATION_WORK_FIELDS = (
        "time", "step", "event_type", "entity_id", "grain_i", "grain_j",
        "DeltaG0", "effective_DeltaG", "capillary_pressure",
        "resolved_shear", "free_volume_chemical_potential",
        "activation_volume_normal", "activation_volume_shear",
        "activation_vacancies", "work_capillary", "work_shear",
        "work_free_volume", "work_total_without_tj_residual",
        "shear_state_before_release", "free_volume_deficit_before_release",
    )

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        config = kwargs.get("config", args[0] if args else None)
        if config is None:
            raise TypeError("MigrationClosureSimulation requires a ModelConfig")
        closure = str(config.parameters.get("migration_closure", "gate_only"))
        if closure not in self.VALID_MIGRATION_CLOSURES:
            raise ValueError(
                f"migration_closure must be one of {sorted(self.VALID_MIGRATION_CLOSURES)}, "
                f"got {closure!r}"
            )
        profile = str(config.parameters.get("blocked_gate_profile", "line"))
        if profile not in self.VALID_GATE_PROFILES:
            raise ValueError(
                f"blocked_gate_profile must be one of {sorted(self.VALID_GATE_PROFILES)}, "
                f"got {profile!r}"
            )
        self.migration_closure = closure
        self.blocked_gate_profile = profile
        super().__init__(*args, **kwargs)

        if bool(config.parameters.get("arclength_domains", True)):
            old_domains = dict(self.domains)
            old_tj_domains = dict(self.tj_domains)
            self.tracker = ArclengthEntityTracker(
                self.orientations,
                config.pf.grid_spacing,
                float(config.parameters.get("event_domain_length", 8.0)),
                periodic=config.pf.boundary_conditions == "periodic",
            )
            self.snapshot = self.tracker.update(self.solver.labels)
            self.domains = {
                key: value for key, value in old_domains.items()
                if key in self.snapshot.boundaries
            }
            self.tj_domains = {
                key: value for key, value in old_tj_domains.items()
                if key in self.snapshot.triple_junctions
            }
            self._boundary_to_tjs = self._index_boundary_tjs()
            # Refresh the mobility/driving fields at the same physical time on
            # the connected domains. At initialization the measured displacement
            # is zero, so this does not advance any physical reaction coordinate.
            self._update_physics()

        work_path = self.output_dir / "activation_work.csv"
        work_exists = work_path.exists() and work_path.stat().st_size > 0
        self._activation_work_handle = work_path.open(
            "a", newline="", encoding="utf-8"
        )
        self._activation_work_writer = csv.DictWriter(
            self._activation_work_handle, fieldnames=self.ACTIVATION_WORK_FIELDS
        )
        if not work_exists:
            self._activation_work_writer.writeheader()
            self._activation_work_handle.flush()

    def run(self):
        try:
            return super().run()
        finally:
            if hasattr(self, "_activation_work_handle"):
                self._activation_work_handle.flush()
                self._activation_work_handle.close()

    def _activation_rates(self, domain: DomainPhysics, segment: GBSegment) -> tuple[
        list[DisconnectionMode], np.ndarray, float, np.ndarray, float, np.ndarray
    ]:
        """Evaluate TST rates with separated capillary, shear, and climb work.

        The legacy reduced model added free-volume chemical potential to the
        normal pressure while leaving ``activation_vacancies`` at zero. That
        mixed climb chemistry into the normal step work. Here the activation
        work is the intended decomposition

            p_cap V_n + tau V_tau + Delta_mu_v N_v.
        """
        capillary = self.config.pf.gb_energy * segment.curvature
        candidates = [
            mode for mode in self.modes
            if (mode.family != "easy" if domain.blocked else True)
        ]
        if self.config.parameters.get("barrier_distribution") == "gb_character":
            candidates = assign_barriers(
                candidates,
                "gb_character",
                self.config.seed,
                float(self.config.parameters.get("barrier_mean_ev", 0.5)),
                misorientation=segment.misorientation,
                character_coefficient_ev=float(
                    self.config.parameters.get("barrier_character_coefficient_ev", 0.1)
                ),
            )
        if "mixed_shear_climb_event" in self.config.active_modules:
            candidates = [mode for mode in candidates if mode.delta_s > 0 and mode.delta_q > 0]

        normal = np.asarray(segment.normal, dtype=float)
        tangent = np.asarray((-normal[1], normal[0]))
        stress = None
        if self.full_field is not None and len(segment.points):
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            y = position[0] % self.config.pf.shape[0]
            x = position[1] % self.config.pf.shape[1]
            stress = self.full_field.stress[:, :, y, x]

        burgers = np.asarray([mode.burgers for mode in candidates], dtype=float)
        magnitudes = np.linalg.norm(burgers, axis=1)
        directions = burgers / np.maximum(magnitudes[:, None], np.finfo(float).tiny)
        resolved_shear = domain.shear.internal_shear_stress * (directions @ tangent)
        if stress is not None:
            resolved_shear += np.einsum("mi,ij,j->m", directions, stress, normal)

        vacancy_mu = float(domain.free_volume.chemical_potential)
        barriers = np.fromiter((mode.barrier_ev for mode in candidates), dtype=float)
        prefactors = np.fromiter(
            (mode.site_multiplicity * mode.attempt_frequency for mode in candidates),
            dtype=float,
        )
        work = (
            float(capillary) * np.fromiter(
                (mode.activation_volume_normal for mode in candidates), dtype=float
            )
            + resolved_shear * np.fromiter(
                (mode.activation_volume_shear for mode in candidates), dtype=float
            )
            + vacancy_mu * np.fromiter(
                (mode.activation_vacancies for mode in candidates), dtype=float
            )
        )

        if "tj_burgers_residual" in self.config.active_modules:
            stiffness = float(self.config.parameters.get("tj_residual_stiffness_ev", 1.0))
            if stiffness <= 0 or not np.isfinite(stiffness):
                raise ValueError("tj_residual_stiffness_ev must be finite and positive")
            packet = float(self.config.parameters.get("packet_size", 1.0))
            event_burgers = packet * burgers
            residual_energy_change = np.zeros(len(candidates), dtype=float)
            for tj, sign in self._signed_boundary_tjs(segment):
                residual = np.asarray(tj.residual_burgers, dtype=float)
                updated = residual[None, :] + sign * event_burgers
                residual_energy_change += 0.5 * stiffness * (
                    np.sum(updated * updated, axis=1) - float(residual @ residual)
                )
            work -= residual_energy_change

        effective_barriers = np.maximum(0.0, barriers - work)
        rates = prefactors * np.exp(
            -effective_barriers / (K_B_EV * self.config.pf.temperature)
        )
        return (
            candidates,
            rates,
            float(capillary),
            resolved_shear,
            vacancy_mu,
            effective_barriers,
        )

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
        prior_compatibility_pending = domain.compatibility_pending

        capillary_pressure = float(self.config.pf.gb_energy * segment.curvature)
        work_capillary = capillary_pressure * mode.activation_volume_normal
        work_shear = driving.resolved_shear * mode.activation_volume_shear
        work_free_volume = (
            driving.vacancy_chemical_potential * mode.activation_vacancies
        )
        if hasattr(self, "_activation_work_writer"):
            effective = (
                mode.effective_barrier_ev(driving)
                if effective_barrier_ev is None else float(effective_barrier_ev)
            )
            self._activation_work_writer.writerow({
                "time": self.solver.time if event_time is None else event_time,
                "step": self.solver.step_number,
                "event_type": event_type,
                "entity_id": domain.entity_id,
                "grain_i": segment.grain_i,
                "grain_j": segment.grain_j,
                "DeltaG0": mode.barrier_ev,
                "effective_DeltaG": effective,
                "capillary_pressure": capillary_pressure,
                "resolved_shear": driving.resolved_shear,
                "free_volume_chemical_potential": driving.vacancy_chemical_potential,
                "activation_volume_normal": mode.activation_volume_normal,
                "activation_volume_shear": mode.activation_volume_shear,
                "activation_vacancies": mode.activation_vacancies,
                "work_capillary": work_capillary,
                "work_shear": work_shear,
                "work_free_volume": work_free_volume,
                "work_total_without_tj_residual": (
                    work_capillary + work_shear + work_free_volume
                ),
                "shear_state_before_release": domain.shear.state,
                "free_volume_deficit_before_release": domain.free_volume.deficit,
            })

        record_mode = mode
        if self.migration_closure == "gate_only" and event_type == "climb_quota_completion":
            # Serial climb already accommodated its quota immediately before
            # this summary event. Prevent a second release in the legacy base
            # event writer while retaining the physical activation work above.
            record_mode = replace(
                mode, point_defect_quota=0.0, delta_q=0.0,
                activation_vacancies=0.0,
            )

        super()._record_event(
            domain,
            segment,
            record_mode,
            rate,
            driving,
            event_type,
            delta_length,
            event_time,
            ledger_position=ledger_position,
            field_position=field_position,
            effective_barrier_ev=effective_barrier_ev,
        )

        if event_type == "compatibility_release":
            domain.compatibility_pending = False
        else:
            domain.compatibility_pending = prior_compatibility_pending
        climb_enabled = bool(set(self.config.active_modules).intersection({
            "free_volume", "serial_climb", "nucleation_limited",
            "multihit_nucleation", "exchange_limited", "transport_limited",
            "mixed_shear_climb_event", "independent_and",
        }))
        climb_blocked = (
            climb_enabled
            and domain.free_volume.deficit
            >= float(self.config.parameters.get("climb_trigger_quota", 0.25))
        )
        domain.blocked = domain.compatibility_pending or climb_blocked

        if self.migration_closure == "gate_only":
            domain.normal_displacement_ledger = prior_hidden_displacement
            domain.normal_release_remaining = prior_release_remaining

    def _tj_mode_driving(
        self, tj: TripleJunction, mode: DisconnectionMode,
    ) -> ModeDriving:
        """Average the conjugate GB forces meeting at a triple junction.

        A TJ release samples the connected physical-arclength domains that meet
        at the junction.  Averaging supplies one local transition-state force
        without making the result depend on pixel count or domain subdivision.
        """
        burgers = np.asarray(mode.burgers, dtype=float)
        b_direction = burgers / max(
            float(np.linalg.norm(burgers)), np.finfo(float).tiny
        )
        capillary: list[float] = []
        shear: list[float] = []
        vacancy_mu: list[float] = []

        stress = None
        if self.full_field is not None:
            y, x = np.rint(tj.position).astype(int) % np.asarray(self.config.pf.shape)
            stress = self.full_field.stress[:, :, y, x]

        for boundary_id in sorted(tj.adjoining_boundaries):
            segment = self.snapshot.boundaries.get(boundary_id)
            boundary_domain = self.domains.get(boundary_id)
            if segment is None or boundary_domain is None:
                continue
            normal = np.asarray(segment.normal, dtype=float)
            tangent = np.asarray((-normal[1], normal[0]))
            resolved = (
                boundary_domain.shear.internal_shear_stress
                * float(b_direction @ tangent)
            )
            if stress is not None:
                resolved += float(b_direction @ stress @ normal)
            capillary.append(float(self.config.pf.gb_energy * segment.curvature))
            shear.append(resolved)
            vacancy_mu.append(float(boundary_domain.free_volume.chemical_potential))

        return ModeDriving(
            float(np.mean(capillary)) if capillary else 0.0,
            float(np.mean(shear)) if shear else 0.0,
            float(np.mean(vacancy_mu)) if vacancy_mu else 0.0,
        )

    def _record_tj_activation_work(
        self,
        domain: DomainPhysics,
        tj: TripleJunction,
        mode: DisconnectionMode,
        driving: ModeDriving,
        bare_barrier_ev: float,
        effective_barrier_ev: float,
        event_time: float,
    ) -> None:
        if not hasattr(self, "_activation_work_writer"):
            return
        boundary_domains = [
            self.domains[boundary_id]
            for boundary_id in sorted(tj.adjoining_boundaries)
            if boundary_id in self.domains
        ]
        shear_state = (
            float(np.mean([item.shear.state for item in boundary_domains]))
            if boundary_domains else 0.0
        )
        free_volume_deficit = (
            float(np.mean([item.free_volume.deficit for item in boundary_domains]))
            if boundary_domains else 0.0
        )
        work_capillary = driving.normal_pressure * mode.activation_volume_normal
        work_shear = driving.resolved_shear * mode.activation_volume_shear
        work_free_volume = (
            driving.vacancy_chemical_potential * mode.activation_vacancies
        )
        self._activation_work_writer.writerow({
            "time": event_time,
            "step": self.solver.step_number,
            "event_type": "tj_compatibility_release",
            "entity_id": domain.entity_id,
            "grain_i": tj.grain_ids[0],
            "grain_j": tj.grain_ids[1],
            "DeltaG0": bare_barrier_ev,
            "effective_DeltaG": effective_barrier_ev,
            "capillary_pressure": driving.normal_pressure,
            "resolved_shear": driving.resolved_shear,
            "free_volume_chemical_potential": driving.vacancy_chemical_potential,
            "activation_volume_normal": mode.activation_volume_normal,
            "activation_volume_shear": mode.activation_volume_shear,
            "activation_vacancies": mode.activation_vacancies,
            "work_capillary": work_capillary,
            "work_shear": work_shear,
            "work_free_volume": work_free_volume,
            "work_total_without_tj_residual": (
                work_capillary + work_shear + work_free_volume
            ),
            "shear_state_before_release": shear_state,
            "free_volume_deficit_before_release": free_volume_deficit,
        })

    def _tj_gate_radius_pixels(self) -> int:
        # The corrected closure applies its TJ gate exactly once, below, using
        # tj_correlation_length/grid_spacing.  Ignore the legacy pixel radius.
        return 0

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

    def _apply_physical_tj_gate(self, mobility: np.ndarray) -> None:
        modules = set(self.config.active_modules)
        if not modules.intersection({
            "tj_compatibility", "tj_pinning", "tj_burgers_strict",
            "tj_burgers_residual", "tj_geometric_surrogate",
        }):
            return
        dx = float(self.config.pf.grid_spacing)
        halfwidth = float(self.config.parameters.get(
            "tj_correlation_length", 0.5 * self.config.pf.interface_width
        ))
        if not np.isfinite(halfwidth) or halfwidth < 0:
            raise ValueError("tj_correlation_length must be finite and nonnegative")
        radius = int(np.ceil(halfwidth / dx - 1e-12)) if halfwidth > 0 else 0
        floor = float(self.config.parameters.get("pinned_mobility_fraction", 0.0))
        for key, tj in self.snapshot.triple_junctions.items():
            domain = self.tj_domains.get(key)
            if domain is None or not domain.blocked:
                continue
            y, x = np.rint(tj.position).astype(int) % np.asarray(self.config.pf.shape)
            for oy in range(-radius, radius + 1):
                for ox in range(-radius, radius + 1):
                    mobility[(y + oy) % mobility.shape[0], (x + ox) % mobility.shape[1]] = min(
                        mobility[(y + oy) % mobility.shape[0], (x + ox) % mobility.shape[1]],
                        floor,
                    )

    def _update_physics(self) -> None:
        """Update corrected local kinematics, barriers, and internal stresses."""
        cfg, modules = self.config, set(self.config.active_modules)
        self._boundary_to_tjs = self._index_boundary_tjs()
        mobility = np.ones(cfg.pf.shape)
        self.driving_field.fill(0.0)
        entity_elapsed = self.solver.time - self.previous_entity_time
        current_ids = set(self.snapshot.boundaries)
        self.domains = {
            key: state for key, state in self.domains.items() if key in current_ids
        }

        for key, segment in self.snapshot.boundaries.items():
            if key not in self.domains:
                self.domains[key] = self._new_domain(segment)
            domain = self.domains[key]
            domain.previous_length = segment.length

            ledger_position = (
                json.dumps(segment.points.mean(axis=0).tolist())
                if len(segment.points) else ""
            )
            field_position = (
                tuple(segment.points[len(segment.points) // 2].astype(int))
                if len(segment.points) else None
            )

            measured_curvature, measured_velocity, measured_normal = interface_kinematics(
                self.solver.eta[segment.grain_i],
                self.previous_entity_eta[segment.grain_i],
                segment.points,
                entity_elapsed,
                cfg.pf.grid_spacing,
                periodic=cfg.pf.boundary_conditions == "periodic",
                partner_phase=self.solver.eta[segment.grain_j],
            )
            segment.curvature = measured_curvature
            segment.velocity = measured_velocity
            local_normal_displacement = (
                float(measured_velocity * entity_elapsed)
                if entity_elapsed > 0 and np.isfinite(measured_velocity) else 0.0
            )
            swept_measure = abs(local_normal_displacement) * max(segment.length, cfg.pf.grid_spacing)

            # Preserve legacy bookkeeping fields for checkpoint compatibility,
            # but do not use whole-grain area changes as a local GB displacement.
            grain_i = self.snapshot.grains[segment.grain_i]
            grain_j = self.snapshot.grains[segment.grain_j]
            domain.previous_area_i = grain_i.area
            domain.previous_area_j = grain_j.area
            domain.previous_time = self.solver.time

            if domain.normal_release_remaining * local_normal_displacement > 0.0:
                consumed = min(
                    abs(domain.normal_release_remaining), abs(local_normal_displacement)
                )
                domain.normal_release_remaining -= np.sign(
                    domain.normal_release_remaining
                ) * consumed
                if abs(domain.normal_release_remaining) < 1e-12:
                    domain.normal_release_remaining = 0.0

            release_distance = float(
                cfg.parameters.get("pf_release_displacement", cfg.pf.grid_spacing)
            )
            if release_distance <= 0:
                raise ValueError("pf_release_displacement must be positive")
            if (
                domain.normal_release_remaining == 0.0
                and abs(domain.normal_displacement_ledger) >= release_distance
            ):
                direction = np.sign(domain.normal_displacement_ledger)
                domain.normal_release_remaining = direction * release_distance
                domain.normal_displacement_ledger -= direction * release_distance

            ci = np.asarray(grain_i.centroid)
            cj = np.asarray(grain_j.centroid)
            normal = cj - ci
            if cfg.pf.boundary_conditions == "periodic":
                box = np.asarray(cfg.pf.shape, dtype=float)
                normal -= np.round(normal / box) * box
            normal /= max(np.linalg.norm(normal), np.finfo(float).tiny)
            if np.linalg.norm(measured_normal):
                normal = np.asarray(measured_normal)
            segment.normal = tuple(normal)

            if "shear_memory" in modules or "shear_feedback" in modules:
                beta = float(cfg.parameters.get("easy_beta", 0.35))
                domain.shear.migrate(beta, local_normal_displacement, cfg.pf.time_step)

            if (
                "qiu_reference_shear" in modules
                and self.full_field is not None
                and local_normal_displacement
                and len(segment.points)
            ):
                beta = float(cfg.parameters.get("easy_beta", 0.35))
                tangent = np.asarray((-normal[1], normal[0]))
                displacement = beta * local_normal_displacement * tangent
                strain = 0.5 * (
                    np.outer(displacement, normal) + np.outer(normal, displacement)
                )
                position = tuple(segment.points[len(segment.points) // 2].astype(int))
                self.full_field.add_event(position, strain)

            if modules.intersection({
                "free_volume", "serial_climb", "nucleation_limited",
                "multihit_nucleation", "exchange_limited", "transport_limited",
                "mixed_shear_climb_event", "independent_and",
            }):
                domain.free_volume.require_for_area_change(swept_measure)

            compatibility_trigger = (
                abs(domain.shear.state) >= float(cfg.parameters.get("shear_trigger", 0.25))
                and domain.free_volume.deficit
                >= float(cfg.parameters.get("climb_trigger_quota", 0.25))
            )
            if (
                modules.intersection({"mixed_shear_climb_event", "independent_and"})
                and compatibility_trigger
                and not domain.blocked
            ):
                domain.blocked = True
                domain.compatibility_pending = True
                self._begin_activation_window(domain)
                domain.climb.activate(self.solver.time)

            # A geometric-surrogate *type* no longer implicitly activates GB
            # barriers. The locus is selected explicitly by the module list so
            # a pure TJ case does not secretly contain GB barriers.
            gb_encounter_enabled = bool(modules.intersection({
                "gb_compatibility", "gb_area_point_defect_pinning", "gb_pinning"
            }))
            if (
                self.solver.step_number > 0
                and gb_encounter_enabled
                and not domain.blocked
                and domain.encounter.advance(swept_measure, maximum_events=1)
            ):
                domain.blocked = True
                domain.compatibility_pending = True
                self._begin_activation_window(domain)
                domain.climb.activate(self.solver.time)

            if domain.compatibility_pending and self.solver.step_number > 0:
                (
                    candidates, rates, normal_pressure, resolved_shear,
                    vacancy_mu, effective,
                ) = self._activation_rates(domain, segment)
                total_rate = float(rates.sum())
                completions, hits = self._advance_activation(
                    domain,
                    total_rate,
                    cfg.pf.time_step,
                    self.solver.time - cfg.pf.time_step,
                    stop_after_completion=True,
                )
                self._record_activation_hits(
                    domain, total_rate, hits, segment=segment, position=ledger_position
                )
                if completions:
                    selected = (
                        int(domain.rng.choice(len(candidates), p=rates / total_rate))
                        if total_rate else 0
                    )
                    mode = candidates[selected]
                    driving = self._mode_driving(
                        normal_pressure, resolved_shear, vacancy_mu, selected
                    )
                    self._record_event(
                        domain,
                        segment,
                        mode,
                        total_rate,
                        driving,
                        "compatibility_release",
                        swept_measure,
                        completions[0].time,
                        ledger_position=ledger_position,
                        field_position=field_position,
                        effective_barrier_ev=float(effective[selected]),
                    )
            elif (
                cfg.compatibility_model == "explicit_modes"
                and not gb_encounter_enabled
                and self.solver.step_number > 0
            ):
                if modules.intersection({"tj_burgers_strict", "tj_burgers_residual"}):
                    self._advance_tj_coupled_mode_flux(
                        domain,
                        segment,
                        swept_measure,
                        ledger_position,
                        field_position,
                    )
                else:
                    (
                        candidates, rates, normal_pressure, resolved_shear,
                        vacancy_mu, effective,
                    ) = self._activation_rates(domain, segment)
                    total_rate = float(rates.sum())
                    completions, hits = self._advance_activation(
                        domain,
                        total_rate,
                        cfg.pf.time_step,
                        self.solver.time - cfg.pf.time_step,
                    )
                    self._record_activation_hits(
                        domain, total_rate, hits, segment=segment, position=ledger_position
                    )
                    for completion in completions:
                        selected = (
                            int(domain.rng.choice(len(candidates), p=rates / total_rate))
                            if total_rate else 0
                        )
                        mode = candidates[selected]
                        driving = self._mode_driving(
                            normal_pressure, resolved_shear, vacancy_mu, selected
                        )
                        self._record_event(
                            domain,
                            segment,
                            mode,
                            total_rate,
                            driving,
                            "disconnection_mode",
                            swept_measure,
                            completion.time,
                            ledger_position=ledger_position,
                            field_position=field_position,
                            effective_barrier_ev=float(effective[selected]),
                        )

            self._advance_climb(domain, segment, swept_measure)

            pair_force = float(cfg.parameters.get("easy_beta", 0.35)) * self._boundary_resolved_shear(
                domain, segment
            )
            if domain.normal_release_remaining:
                pair_force += np.sign(domain.normal_release_remaining) * float(
                    cfg.parameters.get("event_normal_pressure", 1.0)
                )
            for yx in segment.points.astype(int):
                y, x = yx % np.asarray(cfg.pf.shape)
                if domain.blocked:
                    mobility[y, x] = min(
                        mobility[y, x],
                        float(cfg.parameters.get("pinned_mobility_fraction", 0.0)),
                    )
                self.driving_field[segment.grain_i, y, x] += pair_force
                self.driving_field[segment.grain_j, y, x] -= pair_force

        if self.particles is not None:
            for segment in self.snapshot.boundaries.values():
                if len(segment.points):
                    contact = self.particles.contacts(segment.points)
                    for yx in segment.points[contact].astype(int):
                        y, x = yx % np.asarray(cfg.pf.shape)
                        mobility[y, x] = 0.0

        self._update_tj_physics(mobility)
        self._apply_physical_tj_gate(mobility)
        self.solver.set_mobility_scale(mobility)
        if self.full_field is not None:
            self.full_field.solve()
        self.previous_entity_eta = self.solver.eta.copy()
        self.previous_entity_time = self.solver.time
        self._apply_diffuse_blocked_gate()
