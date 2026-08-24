from __future__ import annotations

import csv
import json
from dataclasses import replace
from typing import Any

import numpy as np

from grain_growth_pf.climb.area_loss import (
    ConservedDefectInventory,
    best_compatible_tj_velocity,
    gb_excess_volume_density,
    released_point_defect_quota,
    validate_tj_sink_candidate,
)
from grain_growth_pf.climb.exchange import butler_volmer_flux
from grain_growth_pf.climb.serial_cycle import ClimbStage
from grain_growth_pf.climb.transport import diffusivity, transport_time
from grain_growth_pf.disconnections.barriers import assign_barriers
from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving, K_B_EV
from grain_growth_pf.entities.arclength_tracker import ArclengthEntityTracker
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.entities.triple_junction import TripleJunction
from grain_growth_pf.io.event_trace import EventTraceRecorder, trace_entity_key
from grain_growth_pf.mechanics.force_balance import NormalForceBalance, normal_force_balance
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
    DEFECT_HISTORY_FIELDS = (
        "time", "step", "gb_measure", "released_excess_volume",
        "released_defect_quota", "N_required", "N_required_signed",
        "N_accommodated_GB", "N_accommodated_TJ", "N_external",
        "N_active_deficit", "N_retired", "N_stored_signed",
        "conservation_residual", "vacancy_conservation_residual",
        "interstitial_conservation_residual", "max_abs_conservation_residual",
        "GB_accommodation_fraction", "TJ_accommodation_fraction",
        "external_accommodation_fraction", "active_deficit_fraction",
        "retired_inventory_fraction",
    )
    MECHANISM_STATE_FIELDS = (
        "run_id", "time", "step", "entity_type", "entity_id", "grain_ids",
        "G_pending", "T_pending", "C_pending", "blocked", "climb_stage",
        "local_defect_inventory", "GB_sink_activity", "TJ_sink_activity",
    )
    EVENT_TRACE_RELEASE_TYPES = {
        "compatibility_release", "tj_compatibility_release",
        "gb_sink_completion", "tj_sink_completion", "climb_quota_completion",
    }

    def _boundary_force_balance(
        self, domain: DomainPhysics, segment: GBSegment,
    ) -> NormalForceBalance:
        """Return the executed external force plus sharp-interface diagnostics.

        Capillarity is evolved by the diffuse PF free-energy functional.  The
        scalar ``gamma*kappa`` is recorded as its sharp-interface estimate but
        is not added to ``driving_field`` a second time.  Corrected climb
        chemistry changes activation/pinning and has no continuous pressure,
        so ``p_chem`` is exactly zero here.  The event-release impulse is
        retained separately from chemistry.
        """

        event_pressure = 0.0
        if domain.normal_release_remaining:
            event_pressure = np.sign(domain.normal_release_remaining) * float(
                self.config.parameters.get("event_normal_pressure", 1.0)
            )
        return normal_force_balance(
            capillary_pressure=float(self.config.pf.gb_energy * segment.curvature),
            chemical_pressure=0.0,
            beta=float(self.config.parameters.get("easy_beta", 0.35)),
            resolved_shear=self._boundary_resolved_shear(domain, segment),
            event_pressure=float(event_pressure),
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
        modules = set(config.active_modules)
        self.area_loss_enabled = "area_loss_climb" in modules
        if modules.intersection({"gb_defect_sink", "tj_defect_sink"}) and not self.area_loss_enabled:
            raise ValueError("GB/TJ defect sinks require the area_loss_climb source")
        if self.area_loss_enabled and "gb_defect_sink" not in modules:
            raise ValueError("area_loss_climb requires an explicit gb_defect_sink alternative")
        self.defect_inventory = ConservedDefectInventory()
        self.previous_total_gb_measure: float | None = None
        self.last_released_excess_volume = 0.0
        self.last_released_defect_quota = 0.0
        self.max_abs_conservation_residual = 0.0
        self._area_loss_initializing = True
        self._defect_history_handle = None
        self._defect_history_writer = None
        self._mechanism_state_handle = None
        self._mechanism_state_writer = None
        self._last_gb_sink_entities: set[str] = set()
        self._last_tj_sink_entities: set[str] = set()
        self._pending_sink_completions: list[dict[str, Any]] = []
        self.event_trace_enabled = bool(
            config.parameters.get("event_trace_enabled", False)
        )
        self.event_trace_recorder: EventTraceRecorder | None = None
        self._event_trace_checkpoint_state: dict[str, Any] | None = None
        self._event_trace_local_displacement: dict[str, float] = {}
        self._event_trace_records_by_entity: dict[str, list[dict[str, Any]]] = {}
        super().__init__(*args, **kwargs)

        if self.event_trace_enabled:
            resume = bool(kwargs.get("resume", args[2] if len(args) > 2 else False))
            self.event_trace_recorder = EventTraceRecorder(
                self.output_dir,
                run_id=self.run_id,
                pre_steps=int(config.parameters.get("event_trace_pre_steps", 25)),
                post_steps=int(config.parameters.get("event_trace_post_steps", 50)),
                stride=int(config.parameters.get("event_trace_stride", 1)),
                output_format=str(config.parameters.get("event_trace_format", "csv")),
                resume=resume,
                checkpoint_state=self._event_trace_checkpoint_state,
            )
            self.ledger.observer = self._observe_event_record

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

        if self.area_loss_enabled:
            self.previous_total_gb_measure = self._total_gb_measure()
            history_path = self.output_dir / "defect_inventory.csv"
            history_exists = history_path.exists() and history_path.stat().st_size > 0
            self._defect_history_handle = history_path.open(
                "a", newline="", encoding="utf-8"
            )
            self._defect_history_writer = csv.DictWriter(
                self._defect_history_handle, fieldnames=self.DEFECT_HISTORY_FIELDS
            )
            if not history_exists:
                self._defect_history_writer.writeheader()
            self._area_loss_initializing = False
            self._write_defect_history()
        else:
            self._area_loss_initializing = False

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

        state_path = self.output_dir / "mechanism_state.csv"
        state_exists = state_path.exists() and state_path.stat().st_size > 0
        self._mechanism_state_handle = state_path.open(
            "a", newline="", encoding="utf-8"
        )
        self._mechanism_state_writer = csv.DictWriter(
            self._mechanism_state_handle, fieldnames=self.MECHANISM_STATE_FIELDS
        )
        if not state_exists:
            self._mechanism_state_writer.writeheader()
        self._write_mechanism_state()

    def run(self):
        try:
            return super().run()
        finally:
            if hasattr(self, "_activation_work_handle"):
                self._activation_work_handle.flush()
                self._activation_work_handle.close()
            if self._defect_history_handle is not None:
                self._defect_history_handle.flush()
                self._defect_history_handle.close()
            if self._mechanism_state_handle is not None:
                self._mechanism_state_handle.flush()
                self._mechanism_state_handle.close()
            if self.event_trace_recorder is not None:
                self.event_trace_recorder.close()

    def _total_gb_measure(self) -> float:
        return float(sum(segment.length for segment in self.snapshot.boundaries.values()))

    def _extra_checkpoint_state(self) -> dict[str, Any]:
        state = super()._extra_checkpoint_state()
        if self.area_loss_enabled:
            state["area_loss_climb"] = {
                "inventory": self.defect_inventory.state_dict(),
                "previous_total_gb_measure": self.previous_total_gb_measure,
                "last_released_excess_volume": self.last_released_excess_volume,
                "last_released_defect_quota": self.last_released_defect_quota,
                "max_abs_conservation_residual": self.max_abs_conservation_residual,
            }
        if self.event_trace_recorder is not None:
            state["event_trace"] = self.event_trace_recorder.checkpoint_state()
        return state

    def _load_extra_checkpoint_state(self, state: dict[str, Any]) -> None:
        super()._load_extra_checkpoint_state(state)
        area_loss = state.get("area_loss_climb")
        event_trace = state.get("event_trace")
        self._event_trace_checkpoint_state = (
            dict(event_trace) if event_trace is not None else None
        )
        if not self.area_loss_enabled:
            return
        if area_loss is None:
            raise ValueError("area-loss climb checkpoint lacks conserved inventory")
        self.defect_inventory = ConservedDefectInventory.from_state_dict(
            dict(area_loss["inventory"])
        )
        previous = area_loss.get("previous_total_gb_measure")
        self.previous_total_gb_measure = None if previous is None else float(previous)
        self.last_released_excess_volume = float(
            area_loss.get("last_released_excess_volume", 0.0)
        )
        self.last_released_defect_quota = float(
            area_loss.get("last_released_defect_quota", 0.0)
        )
        self.max_abs_conservation_residual = float(
            area_loss.get("max_abs_conservation_residual", 0.0)
        )

    def _gb_arclength_coordinate(self, entity_id: str) -> float | None:
        segment = self.snapshot.boundaries.get(entity_id)
        if segment is None:
            return None
        pair = tuple(sorted((segment.grain_i, segment.grain_j)))
        candidates = sorted(
            (
                other for other in self.snapshot.boundaries.values()
                if tuple(sorted((other.grain_i, other.grain_j))) == pair
            ),
            key=lambda other: (other.segment_id, other.entity_id),
        )
        cursor = 0.0
        for other in candidates:
            if other.entity_id == entity_id:
                return cursor + 0.5 * float(other.length)
            cursor += float(other.length)
        return None

    def _event_trace_neighborhood(self, entity_id: str) -> set[str]:
        neighborhood: set[str] = set()
        if entity_id.startswith("tj:"):
            tj = self.snapshot.triple_junctions.get(entity_id)
            if tj is None:
                return {trace_entity_key("TJ", entity_id)}
            neighborhood.add(trace_entity_key("TJ", entity_id))
            neighborhood.update(
                trace_entity_key("GB", boundary_id)
                for boundary_id in tj.adjoining_boundaries
                if boundary_id in self.snapshot.boundaries
            )
            return neighborhood

        neighborhood.add(trace_entity_key("GB", entity_id))
        for tj in self._boundary_to_tjs.get(entity_id, ()):
            neighborhood.add(trace_entity_key("TJ", tj.entity_id))
            neighborhood.update(
                trace_entity_key("GB", boundary_id)
                for boundary_id in tj.adjoining_boundaries
                if boundary_id in self.snapshot.boundaries
            )
        return neighborhood

    def _observe_event_record(self, record: dict[str, Any]) -> None:
        if (
            self.event_trace_recorder is None
            or str(record.get("event_type", "")) not in self.EVENT_TRACE_RELEASE_TYPES
        ):
            return
        entity_id = str(record.get("entity_id", ""))
        enriched = dict(record)
        if entity_id.startswith("gb:"):
            enriched["gb_arclength_coordinate"] = self._gb_arclength_coordinate(entity_id)
            enriched["adjoining_tj_ids"] = [
                tj.entity_id for tj in self._boundary_to_tjs.get(entity_id, ())
            ]
        else:
            tj = self.snapshot.triple_junctions.get(entity_id)
            enriched["adjoining_tj_ids"] = [entity_id] if tj is not None else []
        self._event_trace_records_by_entity.setdefault(entity_id, []).append(enriched)
        self.event_trace_recorder.trigger(
            enriched, self._event_trace_neighborhood(entity_id)
        )

    def _event_trace_candidate_keys(self) -> set[str]:
        keys = (
            set(self.event_trace_recorder.requested_keys)
            if self.event_trace_recorder is not None else set()
        )
        for entity_id, domain in self.domains.items():
            if (
                domain.blocked or domain.compatibility_pending or domain.area_loss_pending
                or domain.climb.stage not in {ClimbStage.INACTIVE, ClimbStage.COMPLETE}
            ):
                keys.update(self._event_trace_neighborhood(entity_id))
        for entity_id, domain in self.tj_domains.items():
            if (
                domain.blocked or domain.compatibility_pending or domain.area_loss_pending
                or domain.climb.stage not in {ClimbStage.INACTIVE, ClimbStage.COMPLETE}
            ):
                keys.update(self._event_trace_neighborhood(entity_id))
        return keys

    def _event_trace_inventory(self) -> dict[str, float]:
        if self.area_loss_enabled:
            diagnostics = self.defect_inventory.diagnostics()
            return {
                "N_required": float(diagnostics["N_required"]),
                "N_accommodated_GB": float(diagnostics["N_accommodated_GB"]),
                "N_accommodated_TJ": float(diagnostics["N_accommodated_TJ"]),
                "N_stored": float(self.defect_inventory.stored_total),
                "conservation_residual": float(diagnostics["conservation_residual"]),
            }
        return {
            "N_required": float(sum(d.free_volume.required_total for d in self.domains.values())),
            "N_accommodated_GB": float(sum(d.free_volume.accommodated_total for d in self.domains.values())),
            "N_accommodated_TJ": 0.0,
            "N_stored": float(sum(d.free_volume.deficit for d in self.domains.values())),
            "conservation_residual": 0.0,
        }

    def _event_trace_rows(self, requested: set[str]) -> dict[str, dict[str, Any]]:
        rows: dict[str, dict[str, Any]] = {}
        inventory = self._event_trace_inventory()
        stiffness = float(self.config.parameters.get("free_volume_stiffness", 0.05))
        global_mu = (
            stiffness * self.defect_inventory.stored_signed
            if self.area_loss_enabled else 0.0
        )
        elapsed = max(self.solver.time - self.previous_entity_time, 0.0)
        for entity_id, segment in self.snapshot.boundaries.items():
            key = trace_entity_key("GB", entity_id)
            if key not in requested:
                continue
            domain = self.domains.get(entity_id)
            if domain is None:
                continue
            local_signed = (
                self.defect_inventory.active_vacancy.get(entity_id, 0.0)
                - self.defect_inventory.active_interstitial.get(entity_id, 0.0)
                if self.area_loss_enabled else float(domain.free_volume.deficit)
            )
            velocity = float(segment.velocity)
            capillary_work_rate = float(
                self.config.pf.gb_energy * segment.curvature * velocity
            )
            resolved_shear = float(self._boundary_resolved_shear(domain, segment))
            force_balance = self._boundary_force_balance(domain, segment)
            shear_work_rate = float(resolved_shear * abs(velocity))
            chemical_work = float(global_mu * local_signed)
            latest = self._event_trace_records_by_entity.get(entity_id, [{}])[-1]
            adjoining = [tj.entity_id for tj in self._boundary_to_tjs.get(entity_id, ())]
            rows[key] = {
                "time": float(self.solver.time),
                "temperature": float(self.config.pf.temperature), "seed": int(self.config.seed),
                "entity_type": "GB", "entity_id": entity_id,
                "grain_ids": f"{segment.grain_i};{segment.grain_j}",
                "position": json.dumps(segment.points.mean(axis=0).tolist()) if len(segment.points) else "",
                "gb_arclength_coordinate": self._gb_arclength_coordinate(entity_id),
                "adjoining_tj_ids": json.dumps(adjoining),
                "local_normal_velocity": velocity,
                "local_signed_normal_displacement": self._event_trace_local_displacement.get(
                    entity_id, velocity * elapsed
                ),
                "gb_length": float(segment.length), "domain_length": float(segment.length),
                "curvature": float(segment.curvature),
                "G_pending": int(domain.compatibility_pending), "T_pending": 0,
                "C_pending": int(domain.area_loss_pending),
                "climb_stage": domain.climb.stage.value,
                "shear_state_s": float(domain.shear.state),
                "tau_int": float(domain.shear.internal_shear_stress),
                "free_volume_signed_inventory": float(local_signed),
                "p_cap": force_balance.p_cap,
                "p_chem": force_balance.p_chem,
                "p_shear": force_balance.p_shear,
                "p_event": force_balance.p_event,
                "p_net": force_balance.p_net,
                "chi_s": force_balance.chi_s,
                "local_shear_energy": float(domain.shear.energy),
                "p_cap_V_n": capillary_work_rate, "tau_V_tau": shear_work_rate,
                "Delta_mu_v_N_v": chemical_work,
                "W_total": capillary_work_rate + shear_work_rate + chemical_work,
                "DeltaG0": latest.get("DeltaG0", ""),
                "DeltaG_eff": latest.get("effective_DeltaG", ""),
                "instantaneous_rate": latest.get("instantaneous_rate", ""),
                "cumulative_hazard": float(domain.activation.clock.cumulative_hazard),
                "hazard_threshold": float(domain.activation.clock.threshold),
                "tj_velocity": "", "tj_compatibility_residual": "",
                "residual_burgers": "", "tj_sink_admissible": "",
                "signed_point_defect_flux": float(local_signed),
                **inventory,
            }
        for entity_id, tj in self.snapshot.triple_junctions.items():
            key = trace_entity_key("TJ", entity_id)
            if key not in requested:
                continue
            domain = self.tj_domains.get(entity_id)
            if domain is None:
                continue
            velocity = np.zeros(2, dtype=float)
            residual = np.asarray([], dtype=float)
            try:
                _, velocity, residual, _, _ = self._tj_kinematics(tj)
            except ValueError:
                pass
            latest = self._event_trace_records_by_entity.get(entity_id, [{}])[-1]
            local_velocity = float(np.linalg.norm(velocity))
            rows[key] = {
                "time": float(self.solver.time),
                "temperature": float(self.config.pf.temperature), "seed": int(self.config.seed),
                "entity_type": "TJ", "entity_id": entity_id,
                "grain_ids": ";".join(map(str, tj.grain_ids)),
                "position": json.dumps(np.asarray(tj.position, dtype=float).tolist()),
                "gb_arclength_coordinate": "", "adjoining_tj_ids": json.dumps([entity_id]),
                "local_normal_velocity": local_velocity,
                "local_signed_normal_displacement": local_velocity * elapsed,
                "gb_length": "", "domain_length": "", "curvature": "",
                "G_pending": 0,
                "T_pending": int(domain.blocked and not domain.area_loss_pending),
                "C_pending": int(domain.area_loss_pending),
                "climb_stage": domain.climb.stage.value,
                "shear_state_s": float(domain.shear.state),
                "tau_int": float(domain.shear.internal_shear_stress),
                "free_volume_signed_inventory": (
                    float(self.defect_inventory.stored_signed) if self.area_loss_enabled else 0.0
                ),
                "p_cap": "", "p_chem": "", "p_shear": "", "p_event": "",
                "p_net": "", "chi_s": "", "local_shear_energy": float(domain.shear.energy),
                "p_cap_V_n": "", "tau_V_tau": "",
                "Delta_mu_v_N_v": float(global_mu * self.defect_inventory.stored_signed)
                if self.area_loss_enabled else 0.0,
                "W_total": latest.get("work_total", ""),
                "DeltaG0": latest.get("DeltaG0", ""),
                "DeltaG_eff": latest.get("effective_DeltaG", ""),
                "instantaneous_rate": latest.get("instantaneous_rate", ""),
                "cumulative_hazard": float(domain.activation.clock.cumulative_hazard),
                "hazard_threshold": float(domain.activation.clock.threshold),
                "tj_velocity": json.dumps(velocity.tolist()),
                "tj_compatibility_residual": json.dumps(residual.tolist()),
                "residual_burgers": json.dumps(tj.residual_burgers.tolist()),
                "tj_sink_admissible": latest.get("candidate_allowed", ""),
                "signed_point_defect_flux": (
                    float(self.defect_inventory.stored_signed) if self.area_loss_enabled else 0.0
                ),
                **inventory,
            }
        return rows

    def _capture_event_trace(self) -> None:
        if self.event_trace_recorder is None:
            return
        requested = self._event_trace_candidate_keys()
        self.event_trace_recorder.capture(
            int(self.solver.step_number), self._event_trace_rows(requested)
        )

    def _write_defect_history(self) -> None:
        if self._defect_history_writer is None:
            return
        self.defect_inventory.assert_conserved()
        residual = abs(self.defect_inventory.conservation_residual)
        self.max_abs_conservation_residual = max(
            self.max_abs_conservation_residual, residual
        )
        row = {
            "time": self.solver.time,
            "step": self.solver.step_number,
            "gb_measure": self._total_gb_measure(),
            "released_excess_volume": self.last_released_excess_volume,
            "released_defect_quota": self.last_released_defect_quota,
            "max_abs_conservation_residual": self.max_abs_conservation_residual,
            **self.defect_inventory.diagnostics(),
        }
        self._defect_history_writer.writerow(row)
        self._defect_history_handle.flush()

    def _write_tracks(self) -> None:
        super()._write_tracks()
        if self.area_loss_enabled:
            self._write_defect_history()
        self._write_mechanism_state()

    def _write_mechanism_state(self) -> None:
        if self._mechanism_state_writer is None:
            return
        trigger = float(self.config.parameters.get("climb_trigger_quota", 0.25))
        for key, segment in sorted(self.snapshot.boundaries.items()):
            domain = self.domains.get(key)
            if domain is None:
                continue
            local_inventory = (
                self.defect_inventory.active_vacancy.get(key, 0.0)
                + self.defect_inventory.active_interstitial.get(key, 0.0)
                if self.area_loss_enabled else max(domain.free_volume.deficit, 0.0)
            )
            self._mechanism_state_writer.writerow({
                "run_id": self.run_id, "time": self.solver.time,
                "step": self.solver.step_number, "entity_type": "GB",
                "entity_id": key,
                "grain_ids": f"{segment.grain_i};{segment.grain_j}",
                "G_pending": int(domain.compatibility_pending), "T_pending": 0,
                "C_pending": int(domain.area_loss_pending or local_inventory > trigger),
                "blocked": int(domain.blocked), "climb_stage": domain.climb.stage.value,
                "local_defect_inventory": local_inventory,
                "GB_sink_activity": int(key in self._last_gb_sink_entities),
                "TJ_sink_activity": 0,
            })
        for key, tj in sorted(self.snapshot.triple_junctions.items()):
            domain = self.tj_domains.get(key)
            if domain is None:
                continue
            self._mechanism_state_writer.writerow({
                "run_id": self.run_id, "time": self.solver.time,
                "step": self.solver.step_number, "entity_type": "TJ",
                "entity_id": key, "grain_ids": ";".join(map(str, tj.grain_ids)),
                "G_pending": 0,
                "T_pending": int(domain.blocked and not domain.area_loss_pending),
                "C_pending": int(domain.area_loss_pending),
                "blocked": int(domain.blocked), "climb_stage": domain.climb.stage.value,
                "local_defect_inventory": self.defect_inventory.stored_total if self.area_loss_enabled else 0.0,
                "GB_sink_activity": 0,
                "TJ_sink_activity": int(key in self._last_tj_sink_entities),
            })
        self._mechanism_state_handle.flush()

    def _register_area_loss_source(self, current_ids: set[str]) -> None:
        """Create signed defect demand only from material-wide GB measure loss."""
        current_measure = self._total_gb_measure()
        if self._area_loss_initializing or self.previous_total_gb_measure is None:
            self.previous_total_gb_measure = current_measure
            self.last_released_excess_volume = 0.0
            self.last_released_defect_quota = 0.0
            return

        parameters = self.config.parameters
        density = gb_excess_volume_density(
            excess_volume_per_area=(
                float(parameters["excess_volume_per_area"])
                if "excess_volume_per_area" in parameters else None
            ),
            delta_gb=parameters.get("gb_excess_width"),
            rho_lattice=parameters.get("rho_lattice"),
            rho_gb=parameters.get("rho_gb"),
        )
        formation_volume = float(parameters.get("point_defect_formation_volume", 0.01))
        alpha = float(parameters.get("gb_excess_volume_alpha", 1.0))
        quota = released_point_defect_quota(
            self.previous_total_gb_measure,
            current_measure,
            density,
            formation_volume,
            alpha,
        )
        self.last_released_excess_volume = quota * formation_volume
        self.last_released_defect_quota = quota
        self.defect_inventory.retire_missing(current_ids)

        if quota > 0.0:
            local_losses = {
                key: max(self.domains[key].previous_length - segment.length, 0.0)
                for key, segment in self.snapshot.boundaries.items()
                if key in self.domains and self.domains[key].previous_length > 0.0
            }
            local_total = sum(local_losses.values())
            if local_total > 0.0:
                for key, loss in sorted(local_losses.items()):
                    if loss > 0.0:
                        self.defect_inventory.require(quota * loss / local_total, key)
            else:
                # Topology may remove an entire domain before a local shortening
                # can be attributed. Keep that demand in a material reservoir.
                self.defect_inventory.require(quota, "material-reservoir")
        self.previous_total_gb_measure = current_measure
        self.max_abs_conservation_residual = max(
            self.max_abs_conservation_residual,
            abs(self.defect_inventory.conservation_residual),
        )

    def _available_flux_sign(self) -> int:
        vacancy = sum(self.defect_inventory.active_vacancy.values()) + self.defect_inventory.retired_vacancy
        interstitial = (
            sum(self.defect_inventory.active_interstitial.values())
            + self.defect_inventory.retired_interstitial
        )
        if vacancy <= 0.0 and interstitial <= 0.0:
            return 0
        return 1 if vacancy >= interstitial else -1

    def _area_loss_stage_rates(
        self,
        domain: DomainPhysics,
        length: float,
        *,
        prefix: str,
        barrier_penalty: float = 0.0,
    ) -> tuple[float, float, float]:
        p = self.config.parameters
        temperature = float(self.config.pf.temperature)

        def parameter(name: str, fallback: str, default: float) -> float:
            return float(p.get(prefix + name, p.get(fallback, default)))

        def arrhenius(prefactor: float, barrier: float) -> float:
            return float(prefactor * np.exp(-max(barrier, 0.0) / (K_B_EV * temperature)))

        nucleation = arrhenius(
            parameter("nucleation_prefactor", "nucleation_prefactor", 1e5),
            parameter("nucleation_barrier_ev", "nucleation_barrier_ev", 0.45)
            + max(float(barrier_penalty), 0.0),
        )
        exchange_current = arrhenius(
            parameter("exchange_prefactor", "exchange_prefactor", 1e5),
            parameter("exchange_barrier_ev", "exchange_barrier_ev", 0.55),
        )
        stiffness = float(p.get("free_volume_stiffness", 0.05))
        chemical_potential = stiffness * self.defect_inventory.stored_signed
        exchange_flux = abs(butler_volmer_flux(
            chemical_potential,
            temperature,
            exchange_current,
            float(p.get(prefix + "exchange_transfer_coefficient", p.get("exchange_transfer_coefficient", 0.5))),
        ))
        exchange = exchange_flux / max(domain.climb.required_quota, np.finfo(float).tiny)
        diffusion = diffusivity(
            temperature,
            parameter("transport_prefactor", "transport_prefactor", 1e5),
            parameter("transport_barrier_ev", "transport_barrier_ev", 0.65),
        )
        path_length = float(p.get(prefix + "transport_length", max(length, self.config.pf.grid_spacing)))
        transport = 1.0 / transport_time(
            path_length,
            diffusion,
            float(p.get(prefix + "transport_geometry_factor", p.get("transport_geometry_factor", 1.0))),
        )
        return nucleation, exchange, transport

    def _record_area_loss_row(
        self,
        *,
        domain: DomainPhysics,
        entity_type: str,
        entity_id: str,
        grain_ids: str,
        position: Any,
        event_type: str,
        event_time: float,
        sink_path: str,
        rate: float = 0.0,
        threshold: float = 0.0,
        stage_residence_time: float = 0.0,
        signed_quota: float = 0.0,
        inventory_before: float | None = None,
        inventory_after: float | None = None,
        compatibility_norm: float = 0.0,
        compatibility_residual: Any = "",
        requested_velocity: Any = "",
        compatible_velocity: Any = "",
        burgers_before: Any = "",
        burgers_after: Any = "",
        residual_energy_change: float = 0.0,
        candidate_allowed: bool = True,
        candidate_reason: str = "",
        work_interfacial: float = 0.0,
        work_tj_residual: float = 0.0,
        work_chemical: float = 0.0,
    ) -> None:
        domain.event_counter += 1
        before = self.defect_inventory.stored_total if inventory_before is None else inventory_before
        after = self.defect_inventory.stored_total if inventory_after is None else inventory_after
        self.ledger.write({
            "run_id": self.run_id,
            "time": event_time,
            "step": self.solver.step_number,
            "temperature": self.config.pf.temperature,
            "seed": self.config.seed,
            "event_id": f"{entity_id}:{domain.event_counter}",
            "event_type": event_type,
            "grain_ids": grain_ids,
            "entity_id": entity_id,
            "position": position,
            "geometry_measure_Q": self._total_gb_measure(),
            "barrier_type": "area_loss_climb",
            "instantaneous_rate": rate,
            "cumulative_hazard": threshold,
            "random_hazard_threshold": threshold,
            "hit_count": 1,
            "required_hits_K": 1,
            "point_defect_quota": abs(signed_quota),
            "Nv": signed_quota,
            "Git_SHA": self.sha,
            "sink_path": sink_path,
            "signed_defect_quota": signed_quota,
            "inventory_before": before,
            "inventory_after": after,
            "conservation_residual": self.defect_inventory.conservation_residual,
            "compatibility_norm": compatibility_norm,
            "compatibility_residual": compatibility_residual,
            "tj_velocity_requested": requested_velocity,
            "tj_velocity_compatible": compatible_velocity,
            "burgers_before": burgers_before,
            "burgers_after": burgers_after,
            "residual_energy_change": residual_energy_change,
            "candidate_allowed": int(candidate_allowed),
            "candidate_reason": candidate_reason,
            "stage_residence_time": stage_residence_time,
            "work_applied": 0.0,
            "work_gb_internal": 0.0,
            "work_tj_residual": work_tj_residual,
            "work_chemical": work_chemical,
            "work_interfacial": work_interfacial,
            "work_total": work_interfacial + work_tj_residual + work_chemical,
        })

    def _record_sink_transitions(
        self,
        domain: DomainPhysics,
        entity_type: str,
        entity_id: str,
        grain_ids: str,
        position: Any,
        sink_path: str,
        rates: tuple[float, float, float],
    ) -> None:
        destination = {
            ClimbStage.EXCHANGE: ("nucleation", rates[0]),
            ClimbStage.TRANSPORT: ("exchange", rates[1]),
            ClimbStage.COMPLETE: ("transport", rates[2]),
        }
        history = domain.climb.history
        search_start = max(0, len(history) - len(domain.climb.last_transitions) - 1)
        for transition_time, transition_stage, threshold in domain.climb.last_transitions:
            stage, rate = destination[transition_stage]
            transition_index = next(
                index for index in range(search_start + 1, len(history))
                if history[index][1] == transition_stage
                and np.isclose(history[index][0], transition_time)
            )
            stage_start_time = float(history[transition_index - 1][0])
            actual_residence = max(0.0, float(transition_time) - stage_start_time)
            search_start = transition_index
            self._record_area_loss_row(
                domain=domain,
                entity_type=entity_type,
                entity_id=entity_id,
                grain_ids=grain_ids,
                position=position,
                event_type=f"{sink_path}_{stage}",
                event_time=transition_time,
                sink_path=sink_path,
                rate=rate,
                threshold=threshold,
                stage_residence_time=actual_residence,
            )

    def _select_signed_disconnection(self, flux_sign: int) -> DisconnectionMode:
        candidates = [
            mode for mode in self.modes
            if int(np.sign(mode.point_defect_quota)) == int(np.sign(flux_sign))
        ]
        if not candidates:
            raise RuntimeError(f"no disconnection mode carries defect sign {flux_sign}")
        return min(candidates, key=lambda mode: (mode.barrier_ev, mode.mode_id))

    def _advance_area_loss_gb_sink(
        self, domain: DomainPhysics, segment: GBSegment
    ) -> None:
        if (
            self._area_loss_initializing
            or not self.area_loss_enabled
            or "gb_defect_sink" not in self.config.active_modules
        ):
            return
        trigger = float(self.config.parameters.get("climb_trigger_quota", 0.25))
        local = (
            self.defect_inventory.active_vacancy.get(domain.entity_id, 0.0)
            + self.defect_inventory.active_interstitial.get(domain.entity_id, 0.0)
        )
        reservoir = (
            self.defect_inventory.retired_inventory
            + self.defect_inventory.active_vacancy.get("material-reservoir", 0.0)
            + self.defect_inventory.active_interstitial.get("material-reservoir", 0.0)
        )
        first_domain = min(self.snapshot.boundaries) if self.snapshot.boundaries else ""
        eligible = local > trigger or (reservoir > trigger and domain.entity_id == first_domain)
        if not eligible:
            domain.area_loss_pending = False
            domain.blocked = domain.compatibility_pending
            return

        flux_sign = self._available_flux_sign()
        if flux_sign == 0:
            return
        domain.area_loss_pending = True
        domain.blocked = True
        if domain.climb.stage in {ClimbStage.INACTIVE, ClimbStage.COMPLETE}:
            domain.climb.activate(self.solver.time - self.config.pf.time_step)
        rates = self._area_loss_stage_rates(
            domain, segment.length, prefix="gb_sink_"
        )
        complete = domain.climb.advance(
            self.config.pf.time_step,
            self.solver.time - self.config.pf.time_step,
            *rates,
        )
        position = segment.points.mean(axis=0).tolist() if len(segment.points) else ""
        self._record_sink_transitions(
            domain,
            "GB",
            domain.entity_id,
            f"{segment.grain_i};{segment.grain_j}",
            position,
            "gb_disconnection_sink",
            rates,
        )
        if not complete:
            return
        self._pending_sink_completions.append({
            "event_time": domain.climb.last_completion_time or self.solver.time,
            "sink": "gb", "flux_sign": flux_sign, "domain": domain,
            "segment": segment, "position": position,
        })

    def _tj_kinematics(
        self, tj: TripleJunction
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
        segments = [
            self.snapshot.boundaries[key]
            for key in sorted(tj.adjoining_boundaries)
            if key in self.snapshot.boundaries
        ]
        if len(segments) < 2:
            raise ValueError("TJ has fewer than two resolved adjoining GB branches")
        normals = np.asarray([segment.normal for segment in segments], dtype=float)
        desired = np.asarray([segment.velocity for segment in segments], dtype=float)
        weights = np.asarray([max(segment.length, self.config.pf.grid_spacing) for segment in segments])
        compatible, residual, norm = best_compatible_tj_velocity(normals, desired, weights)
        requested = np.sum(weights[:, None] * desired[:, None] * normals, axis=0) / weights.sum()

        force = np.zeros(2, dtype=float)
        position = np.asarray(tj.position, dtype=float)
        shape = np.asarray(self.config.pf.shape, dtype=float)
        for segment in segments:
            delta = np.asarray(segment.points, dtype=float) - position
            if self.config.pf.boundary_conditions == "periodic":
                delta -= np.round(delta / shape) * shape
            endpoint = delta[np.argmax(np.linalg.norm(delta, axis=1))]
            endpoint /= max(np.linalg.norm(endpoint), np.finfo(float).tiny)
            force += float(self.config.pf.gb_energy) * endpoint
        return requested, compatible, residual, norm, force

    def _advance_area_loss_tj_sinks(self) -> None:
        if (
            self._area_loss_initializing
            or not self.area_loss_enabled
            or "tj_defect_sink" not in self.config.active_modules
        ):
            return
        current = set(self.snapshot.triple_junctions)
        self.tj_domains = {
            key: value for key, value in self.tj_domains.items() if key in current
        }
        trigger = float(self.config.parameters.get("climb_trigger_quota", 0.25))
        if self.defect_inventory.stored_total <= trigger:
            for domain in self.tj_domains.values():
                domain.area_loss_pending = False
            return

        flux_sign = self._available_flux_sign()
        for key, tj in sorted(self.snapshot.triple_junctions.items()):
            if key not in self.tj_domains:
                fake = GBSegment(tj.grain_ids[0], tj.grain_ids[1], 0)
                fake.points = np.asarray([tj.position])
                fake.length = self.config.pf.grid_spacing
                domain = self._new_domain(fake)
                domain.entity_id = key
                self.tj_domains[key] = domain
            domain = self.tj_domains[key]
            candidate_cadence = max(
                1, int(self.config.parameters.get("candidate_diagnostic_cadence", 1))
            )
            candidate_diagnostic_due = self.solver.step_number % candidate_cadence == 0
            adjacent_inventory = {
                boundary_id: (
                    self.defect_inventory.active_vacancy.get(boundary_id, 0.0)
                    + self.defect_inventory.active_interstitial.get(boundary_id, 0.0)
                )
                for boundary_id in tj.adjoining_boundaries
            }
            preferred_entity = max(
                adjacent_inventory,
                key=lambda boundary_id: (adjacent_inventory[boundary_id], boundary_id),
                default=None,
            )
            preferred_amount = (
                adjacent_inventory.get(preferred_entity, 0.0)
                if preferred_entity is not None else 0.0
            )
            reservoir = (
                self.defect_inventory.retired_inventory
                + self.defect_inventory.active_vacancy.get("material-reservoir", 0.0)
                + self.defect_inventory.active_interstitial.get("material-reservoir", 0.0)
            )
            fallback_key = min(self.snapshot.triple_junctions) if self.snapshot.triple_junctions else ""
            if preferred_amount <= trigger and not (reservoir > trigger and key == fallback_key):
                domain.area_loss_pending = False
                domain.blocked = domain.compatibility_pending
                continue
            if preferred_amount <= trigger:
                preferred_entity = "material-reservoir"
                preferred_amount = reservoir
            try:
                requested, compatible, residual, compatibility_norm, capillary_force = self._tj_kinematics(tj)
            except ValueError as exc:
                if candidate_diagnostic_due:
                    self._record_area_loss_row(
                        domain=domain, entity_type="TJ", entity_id=key,
                        grain_ids=";".join(map(str, tj.grain_ids)), position=tj.position,
                        event_type="tj_sink_candidate_rejected", event_time=self.solver.time,
                        sink_path="tj_sink", candidate_allowed=False,
                        candidate_reason=str(exc),
                    )
                continue
            direction = np.asarray(compatible, dtype=float)
            if np.linalg.norm(direction) <= np.finfo(float).tiny:
                direction = -capillary_force
            if np.linalg.norm(direction) <= np.finfo(float).tiny:
                direction = np.asarray([1.0, 0.0])
            direction /= np.linalg.norm(direction)
            step_length = float(self.config.parameters.get("tj_sink_step_length", self.config.pf.grid_spacing))
            displacement = flux_sign * step_length * direction
            burgers_increment = (
                flux_sign * float(self.config.parameters.get("tj_sink_burgers", 0.25))
                * direction
            )
            tolerance = float(self.config.parameters.get("tj_sink_compatibility_tolerance", 0.25))
            stiffness = float(self.config.parameters.get("tj_residual_stiffness_ev", 1.0))
            strict_tolerance = (
                float(self.config.parameters["tj_sink_burgers_tolerance"])
                if "tj_sink_burgers_tolerance" in self.config.parameters else None
            )
            decision = validate_tj_sink_candidate(
                available_signed_quota=self.defect_inventory.stored_signed,
                requested_sign=flux_sign,
                compatibility_norm=compatibility_norm,
                compatibility_tolerance=tolerance,
                residual_burgers=tj.residual_burgers,
                burgers_increment=burgers_increment,
                residual_stiffness=stiffness,
                strict_burgers_tolerance=strict_tolerance,
            )
            interfacial_change = float(-capillary_force @ displacement)
            common = dict(
                domain=domain, entity_type="TJ", entity_id=key,
                grain_ids=";".join(map(str, tj.grain_ids)), position=tj.position,
                sink_path="tj_sink", compatibility_norm=compatibility_norm,
                compatibility_residual=residual.tolist(),
                requested_velocity=requested.tolist(), compatible_velocity=compatible.tolist(),
                burgers_before=tj.residual_burgers.tolist(),
                burgers_after=(tj.residual_burgers + burgers_increment).tolist(),
                residual_energy_change=decision.residual_energy_change,
                work_interfacial=interfacial_change,
                work_tj_residual=decision.residual_energy_change,
            )
            if not decision.allowed:
                domain.area_loss_pending = False
                if candidate_diagnostic_due:
                    self._record_area_loss_row(
                        **common, event_type="tj_sink_candidate_rejected",
                        event_time=self.solver.time, candidate_allowed=False,
                        candidate_reason=decision.reason,
                    )
                continue
            domain.area_loss_pending = True
            domain.blocked = True
            if domain.climb.stage in {ClimbStage.INACTIVE, ClimbStage.COMPLETE}:
                domain.climb.activate(self.solver.time - self.config.pf.time_step)
            penalty = max(interfacial_change, 0.0) + max(decision.residual_energy_change, 0.0)
            rates = self._area_loss_stage_rates(
                domain,
                float(self.config.parameters.get("tj_correlation_length", self.config.pf.grid_spacing)),
                prefix="tj_sink_",
                barrier_penalty=penalty,
            )
            if candidate_diagnostic_due:
                self._record_area_loss_row(
                    **common, event_type="tj_sink_candidate", event_time=self.solver.time,
                    rate=rates[0], candidate_allowed=True, candidate_reason=decision.reason,
                )
            complete = domain.climb.advance(
                self.config.pf.time_step,
                self.solver.time - self.config.pf.time_step,
                *rates,
            )
            self._record_sink_transitions(
                domain, "TJ", key, ";".join(map(str, tj.grain_ids)),
                tj.position, "tj_sink", rates,
            )
            if not complete:
                continue
            self._pending_sink_completions.append({
                "event_time": domain.climb.last_completion_time or self.solver.time,
                "sink": "tj", "flux_sign": flux_sign, "domain": domain,
                "tj": tj, "key": key, "burgers_increment": burgers_increment,
                "common": common, "preferred_entity": preferred_entity,
                "preferred_amount": preferred_amount,
            })

    def _resolve_area_loss_sink_completions(self) -> None:
        """Resolve parallel GB/TJ renewal completions in continuous-time order."""
        trigger = float(self.config.parameters.get("climb_trigger_quota", 0.25))
        ordered = sorted(
            self._pending_sink_completions,
            key=lambda item: (
                float(item["event_time"]), str(item["sink"]),
                item["domain"].entity_id,
            ),
        )
        for candidate in ordered:
            domain = candidate["domain"]
            flux_sign = int(candidate["flux_sign"])
            available_sign = self._available_flux_sign()
            before = self.defect_inventory.stored_total
            if available_sign != flux_sign or before <= 0.0:
                if candidate["sink"] == "tj":
                    common = candidate["common"]
                    self._record_area_loss_row(
                        **common,
                        event_type="tj_sink_candidate_lost_competition",
                        event_time=candidate["event_time"],
                        candidate_allowed=False,
                        candidate_reason="parallel_sink_consumed_available_flux",
                    )
                domain.area_loss_pending = False
                domain.blocked = domain.compatibility_pending
                continue

            if candidate["sink"] == "gb":
                segment = candidate["segment"]
                requested = flux_sign * float(
                    self.config.parameters.get("climb_release_quota", 1.0)
                )
                accepted = self.defect_inventory.accommodate(
                    requested, "gb", preferred_entity=domain.entity_id
                )
                after = self.defect_inventory.stored_total
                mode = self._select_signed_disconnection(flux_sign)
                chemical_work = (
                    float(self.config.parameters.get("free_volume_stiffness", 0.05))
                    * self.defect_inventory.stored_signed * accepted
                )
                self._record_area_loss_row(
                    domain=domain, entity_type="GB", entity_id=domain.entity_id,
                    grain_ids=f"{segment.grain_i};{segment.grain_j}",
                    position=candidate["position"], event_type="gb_sink_completion",
                    event_time=candidate["event_time"], sink_path="gb_disconnection_sink",
                    signed_quota=accepted, inventory_before=before, inventory_after=after,
                    burgers_before=[0.0, 0.0],
                    burgers_after=np.asarray(mode.burgers).tolist(),
                    candidate_reason="selected_by_parallel_completion_time",
                    work_chemical=chemical_work,
                )
                if accepted != 0.0:
                    self._last_gb_sink_entities.add(domain.entity_id)
                local_after = (
                    self.defect_inventory.active_vacancy.get(domain.entity_id, 0.0)
                    + self.defect_inventory.active_interstitial.get(domain.entity_id, 0.0)
                )
                domain.area_loss_pending = local_after > trigger
            else:
                tj = candidate["tj"]
                preferred_entity = candidate["preferred_entity"]
                local_available = (
                    self.defect_inventory.active_vacancy.get(preferred_entity, 0.0)
                    + self.defect_inventory.active_interstitial.get(preferred_entity, 0.0)
                )
                if preferred_entity == "material-reservoir":
                    local_available += self.defect_inventory.retired_inventory
                if local_available <= 0.0:
                    common = candidate["common"]
                    self._record_area_loss_row(
                        **common, event_type="tj_sink_candidate_lost_competition",
                        event_time=candidate["event_time"], candidate_allowed=False,
                        candidate_reason="adjacent_gb_flux_consumed_by_competing_sink",
                    )
                    domain.area_loss_pending = False
                    domain.blocked = domain.compatibility_pending
                    continue
                requested_amount = min(local_available, float(self.config.parameters.get(
                    "tj_sink_release_quota",
                    self.config.parameters.get("climb_release_quota", 1.0),
                )))
                requested = flux_sign * requested_amount
                accepted = self.defect_inventory.accommodate(
                    requested, "tj", preferred_entity=preferred_entity
                )
                after = self.defect_inventory.stored_total
                if accepted != 0.0:
                    tj.add_burgers(candidate["burgers_increment"])
                chemical_work = (
                    float(self.config.parameters.get("free_volume_stiffness", 0.05))
                    * self.defect_inventory.stored_signed * accepted
                )
                common = candidate["common"]
                self._record_area_loss_row(
                    **{**common, "burgers_after": tj.residual_burgers.tolist()},
                    event_type="tj_sink_completion", event_time=candidate["event_time"],
                    signed_quota=accepted, inventory_before=before, inventory_after=after,
                    candidate_allowed=True,
                    candidate_reason="selected_by_parallel_completion_time",
                    work_chemical=chemical_work,
                )
                if accepted != 0.0:
                    self._last_tj_sink_entities.add(candidate["key"])
                domain.area_loss_pending = self.defect_inventory.stored_total > trigger
            domain.blocked = domain.compatibility_pending or domain.area_loss_pending
        self._pending_sink_completions.clear()

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
        domain.blocked = (
            domain.compatibility_pending or domain.area_loss_pending or climb_blocked
        )

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
            "tj_burgers_residual", "tj_geometric_surrogate", "tj_defect_sink",
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
        self._event_trace_local_displacement.clear()
        self._event_trace_records_by_entity.clear()
        self._last_gb_sink_entities.clear()
        self._last_tj_sink_entities.clear()
        self._pending_sink_completions.clear()
        self._boundary_to_tjs = self._index_boundary_tjs()
        mobility = np.ones(cfg.pf.shape)
        self.driving_field.fill(0.0)
        entity_elapsed = self.solver.time - self.previous_entity_time
        current_ids = set(self.snapshot.boundaries)
        if self.area_loss_enabled:
            self._register_area_loss_source(current_ids)
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
            self._event_trace_local_displacement[key] = local_normal_displacement
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

            if not self.area_loss_enabled and modules.intersection({
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

            if self.area_loss_enabled:
                self._advance_area_loss_gb_sink(domain, segment)
            else:
                self._advance_climb(domain, segment, swept_measure)

            force_balance = self._boundary_force_balance(domain, segment)
            # Capillarity is already in the PF functional and corrected climb
            # has no continuous chemical pressure.  Only the algebraically
            # identical external shear and event-release terms are applied.
            pair_force = force_balance.p_shear + force_balance.p_event
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
        self._advance_area_loss_tj_sinks()
        self._resolve_area_loss_sink_completions()
        self._apply_physical_tj_gate(mobility)
        self.solver.set_mobility_scale(mobility)
        if self.full_field is not None:
            self.full_field.solve()
        self._capture_event_trace()
        self.previous_entity_eta = self.solver.eta.copy()
        self.previous_entity_time = self.solver.time
        self._apply_diffuse_blocked_gate()
        if self.area_loss_enabled:
            self.defect_inventory.assert_conserved()
