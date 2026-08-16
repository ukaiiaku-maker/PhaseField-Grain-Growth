from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from grain_growth_pf.climb.free_volume import FreeVolumeState
from grain_growth_pf.climb.exchange import butler_volmer_flux
from grain_growth_pf.climb.serial_cycle import SerialClimbCycle
from grain_growth_pf.climb.transport import diffusivity, transport_time
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.disconnections.barriers import assign_barriers
from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving
from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.disconnections.shear_coupling import event_shear_increment, event_volumetric_increment
from grain_growth_pf.disconnections.spectrum import isotropic_surrogate_library
from grain_growth_pf.encounters.geometric_hazard import GeometricEncounterClock
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.entities.tracker import EntityTracker
from grain_growth_pf.io.event_ledger import EventLedger
from grain_growth_pf.io.checkpoints import atomic_savez_compressed, atomic_write_text
from grain_growth_pf.io.provenance import file_sha256, git_sha, write_manifest
from grain_growth_pf.mechanics.local_shear_memory import LocalShearMemory
from grain_growth_pf.mechanics.qiu_full_field import QiuFullField
from grain_growth_pf.obstacles.particles import ParticleField
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.kinematics import interface_kinematics
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from grain_growth_pf.stochastic.multihit import CompletionEvent, MultiHitProcess
from grain_growth_pf.stochastic.hazard import HazardEvent


@dataclass
class DomainPhysics:
    entity_id: str
    rng: np.random.Generator
    encounter_density: float
    hits: int
    hit_interpretation: str
    shear_stiffness: float
    shear_relaxation_time: float | None
    excess_volume: float
    formation_volume: float
    free_volume_stiffness: float
    encounter: GeometricEncounterClock = field(init=False)
    activation: MultiHitProcess = field(init=False)
    shear: LocalShearMemory = field(init=False)
    free_volume: FreeVolumeState = field(init=False)
    climb: SerialClimbCycle = field(init=False)
    blocked: bool = False
    previous_length: float = 0.0
    previous_area_i: float = 0.0
    previous_area_j: float = 0.0
    previous_time: float = 0.0
    normal_displacement_ledger: float = 0.0
    normal_release_remaining: float = 0.0
    packet_window_elapsed: float = 0.0
    event_counter: int = 0

    def __post_init__(self) -> None:
        self.encounter = GeometricEncounterClock(self.encounter_density, self.rng)
        self.activation = MultiHitProcess(self.hits, self.rng, self.hit_interpretation)
        self.shear = LocalShearMemory(self.shear_stiffness, relaxation_time=self.shear_relaxation_time)
        self.free_volume = FreeVolumeState(self.excess_volume, self.formation_volume, self.free_volume_stiffness)
        self.climb = SerialClimbCycle(self.rng)

    def state_dict(self) -> dict[str, Any]:
        return {
            "rng": self.rng.bit_generator.state,
            "encounter": {"cumulative_hazard": self.encounter.cumulative_hazard,
                          "threshold": self.encounter.threshold, "total_measure": self.encounter.total_measure},
            "activation": {"cumulative_hazard": self.activation.clock.cumulative_hazard,
                           "threshold": self.activation.clock.threshold, "last_rate": self.activation.clock.last_rate,
                           "hit_count": self.activation.hit_count,
                           "packet_window_elapsed": self.packet_window_elapsed},
            "shear": {"state": self.shear.state, "dissipated_energy": self.shear.dissipated_energy},
            "free_volume": {"required_total": self.free_volume.required_total,
                            "accommodated_total": self.free_volume.accommodated_total,
                            "dissipated_energy": self.free_volume.dissipated_energy},
            "climb": {"stage": self.climb.stage.value, "required_quota": self.climb.required_quota,
                      "completed_quota": self.climb.completed_quota,
                      "clock_hazard": self.climb.clock.cumulative_hazard,
                      "clock_threshold": self.climb.clock.threshold,
                      "clock_last_rate": self.climb.clock.last_rate,
                      "last_completion_time": self.climb.last_completion_time},
            "blocked": self.blocked, "previous_length": self.previous_length,
            "previous_area_i": self.previous_area_i, "previous_area_j": self.previous_area_j,
            "previous_time": self.previous_time,
            "normal_displacement_ledger": self.normal_displacement_ledger,
            "normal_release_remaining": self.normal_release_remaining,
            "event_counter": self.event_counter,
        }

    def load_state_dict(self, state: dict[str, Any]) -> None:
        from grain_growth_pf.climb.serial_cycle import ClimbStage
        self.rng.bit_generator.state = state["rng"]
        for name, value in state["encounter"].items(): setattr(self.encounter, name, value)
        self.activation.clock.cumulative_hazard = state["activation"]["cumulative_hazard"]
        self.activation.clock.threshold = state["activation"]["threshold"]
        self.activation.clock.last_rate = state["activation"]["last_rate"]
        self.activation.hit_count = state["activation"]["hit_count"]
        self.packet_window_elapsed = state["activation"].get("packet_window_elapsed", 0.0)
        self.shear.state = state["shear"]["state"]
        self.shear.dissipated_energy = state["shear"]["dissipated_energy"]
        self.free_volume.required_total = state["free_volume"]["required_total"]
        self.free_volume.accommodated_total = state["free_volume"]["accommodated_total"]
        self.free_volume.dissipated_energy = state["free_volume"].get("dissipated_energy", 0.0)
        self.climb.stage = ClimbStage(state["climb"]["stage"])
        self.climb.required_quota = state["climb"]["required_quota"]
        self.climb.completed_quota = state["climb"]["completed_quota"]
        self.climb.clock.cumulative_hazard = state["climb"]["clock_hazard"]
        self.climb.clock.threshold = state["climb"]["clock_threshold"]
        self.climb.clock.last_rate = state["climb"]["clock_last_rate"]
        self.climb.last_completion_time = state["climb"].get("last_completion_time")
        self.blocked = state["blocked"]
        self.previous_length = state["previous_length"]
        self.previous_area_i = state.get("previous_area_i", 0.0)
        self.previous_area_j = state.get("previous_area_j", 0.0)
        self.previous_time = state.get("previous_time", 0.0)
        self.normal_displacement_ledger = state["normal_displacement_ledger"]
        self.normal_release_remaining = state.get("normal_release_remaining", 0.0)
        self.event_counter = state["event_counter"]


def _child_rng(seed: int, entity_id: str) -> np.random.Generator:
    # Stable hierarchy independent of dictionary/Python hash ordering.
    words = np.frombuffer(entity_id.encode(), dtype=np.uint8).astype(np.uint32)
    sequence = np.random.SeedSequence([seed, *words.tolist()])
    return np.random.default_rng(sequence)


class EventResolvedSimulation:
    """Couples persistent entity clocks and compatibility state to the PF solver."""

    def __init__(self, config: ModelConfig, output_dir: str | Path, resume: bool = False,
                 code_sha: str | None = None):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=resume)
        self.sha = code_sha or git_sha()
        initial_state_path = config.parameters.get("initial_state_file")
        used_cached_initial_condition = bool(initial_state_path)
        if initial_state_path:
            with np.load(initial_state_path) as state:
                eta = state["eta"].copy()
                seeds = state["seed_positions"].copy()
                orientations = state["orientations"].copy()
                active_original_ids = state["active_original_ids"].astype(int).copy()
                equilibration_steps = int(state["equilibration_steps"])
        else:
            eta, seeds, orientations = voronoi_polycrystal(
                config.pf.shape, int(config.parameters.get("initial_grains", 50)), config.seed,
                width=config.pf.interface_width / 2,
                periodic=config.pf.boundary_conditions == "periodic",
            )
            active_original_ids = np.arange(len(orientations))
            equilibration_steps = int(config.parameters.get("equilibration_steps", 0))
        self.orientations = orientations
        effective_pf = config.pf
        if "arrhenius_intrinsic" in config.active_modules:
            from dataclasses import replace
            barrier = float(config.parameters.get("intrinsic_barrier_ev", 0.45))
            effective_pf = replace(config.pf, intrinsic_mobility=(
                config.pf.intrinsic_mobility * np.exp(-barrier / (K_B_EV * config.pf.temperature))
            ))
        self.driving_field = np.zeros_like(eta)
        self.solver = MultiphaseFieldSolver(eta, effective_pf, driving=None)
        if not resume:
            write_manifest(self.output_dir / "manifest.json", config.to_dict(), "equilibrating", {
                "initial_seed_positions": seeds.tolist(),
                "original_orientations": orientations.tolist(),
            }, code_sha=self.sha)
            if not used_cached_initial_condition:
                for _ in range(equilibration_steps):
                    self.solver.step()
            target_grains = config.parameters.get("equilibrate_to_grains")
            if target_grains is not None and not used_cached_initial_condition:
                target = int(target_grains)
                maximum = int(config.parameters.get("equilibration_max_steps", 5000))
                while np.count_nonzero(self.solver.active_phases) > target and equilibration_steps < maximum:
                    self.solver.step()
                    equilibration_steps += 1
                    if equilibration_steps % 100 == 0:
                        (self.output_dir / "equilibration_progress.json").write_text(json.dumps({
                            "steps": equilibration_steps,
                            "active_grains": int(np.count_nonzero(self.solver.active_phases)),
                            "target_grains": target,
                            "git_sha": self.sha,
                        }, indent=2) + "\n")
                if np.count_nonzero(self.solver.active_phases) > target:
                    raise RuntimeError(
                        f"pre-equilibration retained {np.count_nonzero(self.solver.active_phases)} "
                        f"grains after the configured maximum of {maximum} steps"
                    )
            if bool(config.parameters.get("compact_after_equilibration", True)) and not used_cached_initial_condition:
                active_original_ids = np.flatnonzero(self.solver.active_phases)
                self.solver.eta = self.solver.eta[active_original_ids].copy()
                self.solver.active_phases = np.ones(len(active_original_ids), dtype=bool)
                self.orientations = orientations[active_original_ids].copy()
            self.solver.time = 0.0
            self.solver.step_number = 0
            self.driving_field = np.zeros_like(self.solver.eta)
        self.solver.driving = self._driving
        self.tracker = EntityTracker(
            self.orientations, config.pf.grid_spacing,
            float(config.parameters.get("event_domain_length", 8.0)),
            periodic=config.pf.boundary_conditions == "periodic",
        )
        self.snapshot = self.tracker.update(self.solver.labels)
        self._boundary_to_tjs = self._index_boundary_tjs()
        self.domains: dict[str, DomainPhysics] = {}
        self.tj_domains: dict[str, DomainPhysics] = {}
        self.modes = isotropic_surrogate_library(
            b_shells=tuple(config.parameters.get("b_shells", (0.25, 0.5, 1.0))),
            directions=int(config.parameters.get("mode_directions", 8)),
            step_heights=tuple(config.parameters.get("step_heights", (0.25,))),
            barrier_core_ev=float(config.parameters.get("barrier_core_ev", 0.25)),
            b_coefficient_ev=float(config.parameters.get("b_coefficient_ev", 0.25)),
            h_coefficient_ev=float(config.parameters.get("h_coefficient_ev", 0.10)),
            b_power=float(config.parameters.get("b_power", 2.0)),
            attempt_frequency=float(config.parameters.get("attempt_frequency", 1e2)),
            seed=config.seed,
            disorder_std_ev=float(config.parameters.get("mode_disorder_std_ev", 0.0)),
        )
        barrier_distribution = config.parameters.get("barrier_distribution")
        if barrier_distribution and barrier_distribution != "gb_character":
            raw_bounds = config.parameters.get("barrier_bounds_ev")
            bounds = tuple(map(float, raw_bounds)) if raw_bounds is not None else None
            self.modes = assign_barriers(
                self.modes,
                str(barrier_distribution),
                config.seed + 1771,
                float(config.parameters.get("barrier_mean_ev", 0.5)),
                float(config.parameters.get("barrier_std_ev", 0.1)),
                bounds,
            )
        self.full_field = QiuFullField(config.pf.shape) if config.mechanics_backend == "qiu_full_field" else None
        particle_modules = {"random_spatial_pinning", "particle_zener"}.intersection(config.active_modules)
        self.particles = ParticleField.random(
            int(config.parameters.get("particle_count", 20)),
            float(config.parameters.get("particle_radius", 1.5)), config.pf.shape,
            config.seed + 991,
        ) if particle_modules else None
        event_format = str(config.parameters.get("event_ledger_format", "csv")).lower()
        if event_format not in {"csv", "parquet"}:
            raise ValueError(f"unsupported event_ledger_format {event_format!r}")
        event_name = (
            "events.parquet"
            if event_format == "parquet"
            else "events.csv.gz"
            if bool(config.parameters.get("compress_event_ledger", False))
            else "events.csv"
        )
        self.ledger = EventLedger(self.output_dir / event_name)
        track_path = self.output_dir / "grain_tracks.csv"
        self.track_handle = track_path.open("a" if resume else "w", newline="", encoding="utf-8")
        self.track_writer = csv.DictWriter(self.track_handle, fieldnames=(
            "run_id", "time", "step", "grain_id", "area", "radius", "perimeter", "neighbors"
        ))
        if not resume or track_path.stat().st_size == 0:
            self.track_writer.writeheader()
        boundary_path = self.output_dir / "boundary_tracks.csv"
        self.boundary_handle = boundary_path.open("a" if resume else "w", newline="", encoding="utf-8")
        self.boundary_writer = csv.DictWriter(self.boundary_handle, fieldnames=(
            "run_id", "time", "step", "entity_id", "grain_i", "grain_j", "length",
            "curvature", "normal_velocity", "blocked", "resolved_shear", "free_volume_deficit"
        ))
        if not resume or boundary_path.stat().st_size == 0:
            self.boundary_writer.writeheader()
        self.run_id = self.output_dir.name
        self.energy_records: list[dict[str, float]] = []
        self.accumulated_shear_strain = 0.0
        self.accumulated_volumetric_strain = 0.0
        self.previous_entity_eta = self.solver.eta.copy()
        self.previous_entity_time = self.solver.time
        if resume:
            self._load_checkpoint()
        else:
            self._update_physics()
            write_manifest(self.output_dir / "manifest.json", config.to_dict(), "running", {
                "initial_seed_positions": seeds.tolist(), "orientations": self.orientations.tolist(),
                "original_orientations": orientations.tolist(),
                "equilibration_steps_completed_before_time_zero": equilibration_steps,
                "grains_after_equilibration": len(self.snapshot.grains),
                "active_original_grain_ids": active_original_ids.tolist(),
                "initial_condition_source": str(initial_state_path) if initial_state_path else "generated_in_run",
            }, code_sha=self.sha)
            self._write_tracks()

    def _driving(self, _eta: np.ndarray, _time: float) -> np.ndarray:
        return self.driving_field

    def _new_domain(self, segment: GBSegment) -> DomainPhysics:
        p = self.config.parameters
        modules = set(self.config.active_modules)
        hits = int(p.get("required_hits", 3 if any("multihit" in m for m in modules) else 1))
        interpretation = "packet_reset" if "multihit_packet_reset" in modules else "persistent_hits"
        domain = DomainPhysics(
            segment.entity_id, _child_rng(self.config.seed, segment.entity_id),
            float(p.get("encounter_density", 0.08)), hits, interpretation,
            float(p.get("shear_stiffness", 0.15)),
            (float(p.get("shear_relaxation_time", 1.0)) if "deterministic_relaxation" in modules else None),
            float(p.get("excess_volume_per_area", 0.02)),
            float(p.get("point_defect_formation_volume", 0.01)), float(p.get("free_volume_stiffness", 0.05)),
        )
        domain.climb.required_quota = float(p.get("climb_release_quota", 1.0))
        return domain

    def _begin_activation_window(self, domain: DomainPhysics) -> None:
        domain.activation.begin_window()
        domain.packet_window_elapsed = 0.0

    def _advance_activation(self, domain: DomainPhysics, rate: float, dt: float,
                            start_time: float, stop_after_completion: bool = False,
                            ) -> tuple[list[CompletionEvent], list[tuple[HazardEvent, int, bool]]]:
        """Advance an activation clock with explicit finite packet renewal windows."""
        if domain.activation.interpretation != "packet_reset":
            completions = domain.activation.advance(
                rate, dt, start_time, stop_after_completion
            )
            hits = list(zip(domain.activation.last_hit_events,
                            domain.activation.last_hit_counts,
                            domain.activation.last_hit_completions))
            return completions, hits
        window_time = float(self.config.parameters.get("packet_window_time", 1.0))
        if not np.isfinite(window_time) or window_time <= 0:
            raise ValueError("packet_window_time must be finite and positive")
        completions: list[CompletionEvent] = []
        hits: list[tuple[HazardEvent, int, bool]] = []
        remaining = float(dt)
        current_time = float(start_time)
        tolerance = 16 * np.finfo(float).eps * max(1.0, window_time)
        while remaining > tolerance:
            available = max(window_time - domain.packet_window_elapsed, 0.0)
            if available <= tolerance:
                domain.activation.begin_window()
                domain.packet_window_elapsed = 0.0
                available = window_time
            span = min(remaining, available)
            found = domain.activation.advance(
                rate, span, current_time, stop_after_completion
            )
            completions.extend(found)
            hits.extend(zip(domain.activation.last_hit_events,
                            domain.activation.last_hit_counts,
                            domain.activation.last_hit_completions))
            consumed = (
                max(0.0, found[-1].time - current_time)
                if found and stop_after_completion else span
            )
            domain.packet_window_elapsed += consumed
            current_time += span
            remaining -= span
            if found and stop_after_completion:
                return completions, hits
            if found:
                domain.packet_window_elapsed = max(0.0, current_time - found[-1].time)
            elif domain.packet_window_elapsed >= window_time - tolerance:
                domain.activation.begin_window()
                domain.packet_window_elapsed = 0.0
        return completions, hits

    def _record_activation_hits(
        self, domain: DomainPhysics, rate: float,
        hits: list[tuple[HazardEvent, int, bool]], *,
        segment: GBSegment | None = None, grain_ids: str = "",
        position: Any = "", geometry_measure: float = 0.0,
        tj_travel: float = 0.0,
    ) -> None:
        """Write every stochastic first passage, including sub-completion hits."""
        if segment is not None:
            grain_ids = f"{segment.grain_i};{segment.grain_j}"
            if position == "":
                position = (
                    json.dumps(segment.points.mean(axis=0).tolist())
                    if len(segment.points) else ""
                )
            geometry_measure = domain.encounter.total_measure
        for event, hit_count, completed in hits:
            domain.event_counter += 1
            self.ledger.write({
                "run_id": self.run_id, "time": event.event_time,
                "step": self.solver.step_number,
                "temperature": self.config.pf.temperature, "seed": self.config.seed,
                "event_id": f"{domain.entity_id}:{domain.event_counter}",
                "event_type": "activation_hit" if segment is not None else "tj_activation_hit",
                "grain_ids": grain_ids, "entity_id": domain.entity_id,
                "position": position, "geometry_measure_Q": geometry_measure,
                "curvature": segment.curvature if segment is not None else "",
                "local_velocity": segment.velocity if segment is not None else "",
                "local_normal_free_volume_stress": domain.free_volume.chemical_potential,
                "shear_state_s": domain.shear.state,
                "free_volume_state_q": domain.free_volume.deficit,
                "instantaneous_rate": rate,
                "cumulative_hazard": event.threshold,
                "random_hazard_threshold": event.threshold,
                "hit_count": hit_count, "required_hits_K": domain.hits,
                "release_Delta_s": 0.0, "release_Delta_q": 0.0,
                "GB_area_change": 0.0,
                "TJ_travel": tj_travel,
                "point_defect_quota": domain.free_volume.deficit,
                "normal_step_h": 0.0, "burgers_vector_b": "", "Nv": 0.0,
                "shear_strain_increment": 0.0,
                "volumetric_strain_increment": 0.0,
                "packet_size": self.config.parameters.get("packet_size", 1.0),
                "Git_SHA": self.sha,
            })

    def _activation_rates(self, domain: DomainPhysics, segment: GBSegment) -> tuple[
        list[DisconnectionMode], np.ndarray, float, np.ndarray, float
    ]:
        capillary = self.config.pf.gb_energy * segment.curvature
        candidates = [m for m in self.modes if (m.family != "easy" if domain.blocked else True)]
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
            candidates = [m for m in candidates if m.delta_s > 0 and m.delta_q > 0]
        normal = np.asarray(segment.normal, dtype=float)
        tangent = np.asarray((-normal[1], normal[0]))
        stress = None
        if self.full_field is not None and len(segment.points):
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            y, x = position[0] % self.config.pf.shape[0], position[1] % self.config.pf.shape[1]
            stress = self.full_field.stress[:, :, y, x]
        burgers = np.asarray([mode.burgers for mode in candidates], dtype=float)
        magnitudes = np.linalg.norm(burgers, axis=1)
        directions = burgers / np.maximum(magnitudes[:, None], np.finfo(float).tiny)
        resolved_shear = domain.shear.internal_shear_stress * (directions @ tangent)
        if stress is not None:
            resolved_shear += np.einsum("mi,ij,j->m", directions, stress, normal)
        normal_pressure = capillary + domain.free_volume.chemical_potential
        vacancy_mu = domain.free_volume.chemical_potential
        barriers = np.fromiter((mode.barrier_ev for mode in candidates), dtype=float)
        prefactors = np.fromiter(
            (mode.site_multiplicity * mode.attempt_frequency for mode in candidates),
            dtype=float,
        )
        work = (
            normal_pressure * np.fromiter(
                (mode.activation_volume_normal for mode in candidates), dtype=float
            )
            + resolved_shear * np.fromiter(
                (mode.activation_volume_shear for mode in candidates), dtype=float
            )
            + vacancy_mu * np.fromiter(
                (mode.activation_vacancies for mode in candidates), dtype=float
            )
        )
        effective_barriers = np.maximum(0.0, barriers - work)
        rates = prefactors * np.exp(
            -effective_barriers / (K_B_EV * self.config.pf.temperature)
        )
        return candidates, rates, float(normal_pressure), resolved_shear, float(vacancy_mu)

    @staticmethod
    def _mode_driving(normal_pressure: float, resolved_shear: np.ndarray,
                      vacancy_mu: float, selected: int) -> ModeDriving:
        return ModeDriving(normal_pressure, float(resolved_shear[selected]), vacancy_mu)

    def _activation_mode(self, domain: DomainPhysics, segment: GBSegment) -> tuple[DisconnectionMode, float, ModeDriving]:
        candidates, rates, normal_pressure, resolved_shear, vacancy_mu = self._activation_rates(
            domain, segment
        )
        total = float(rates.sum())
        if total == 0:
            return candidates[0], 0.0, self._mode_driving(
                normal_pressure, resolved_shear, vacancy_mu, 0
            )
        mode_rng = domain.rng
        selected = int(mode_rng.choice(len(candidates), p=rates / total))
        return candidates[selected], total, self._mode_driving(
            normal_pressure, resolved_shear, vacancy_mu, selected
        )

    def _boundary_resolved_shear(self, domain: DomainPhysics, segment: GBSegment) -> float:
        value = domain.shear.internal_shear_stress
        if self.full_field is not None and len(segment.points):
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            normal = np.asarray(segment.normal, dtype=float)
            tangent = np.asarray((-normal[1], normal[0]))
            value += self.full_field.resolved_shear(position, tangent, normal)
        return float(value)

    def _index_boundary_tjs(self) -> dict[str, tuple[Any, ...]]:
        """Index the (normally zero to two) TJs adjoining each GB domain."""
        indexed: dict[str, list[Any]] = {}
        for tj in self.snapshot.triple_junctions.values():
            for boundary_id in tj.adjoining_boundaries:
                indexed.setdefault(boundary_id, []).append(tj)
        return {boundary_id: tuple(tjs) for boundary_id, tjs in indexed.items()}

    def _record_event(self, domain: DomainPhysics, segment: GBSegment, mode: DisconnectionMode,
                      rate: float, driving: ModeDriving, event_type: str, delta_length: float,
                      event_time: float | None = None, *, ledger_position: str | None = None,
                      field_position: tuple[int, int] | None = None) -> None:
        domain.event_counter += 1
        packet = float(self.config.parameters.get("packet_size", 1.0))
        b = np.asarray(mode.burgers) * packet
        h = mode.step_height * packet
        dq = mode.point_defect_quota * packet
        rve_measure = float(np.prod(self.config.pf.shape) * self.config.pf.grid_spacing**2)
        shear_strain = event_shear_increment(float(np.linalg.norm(b)), segment.length, rve_measure)
        volumetric_strain = event_volumetric_increment(dq, domain.formation_volume, rve_measure)
        self.accumulated_shear_strain += shear_strain
        self.accumulated_volumetric_strain += volumetric_strain
        domain.normal_displacement_ledger += h
        if self.full_field is not None and len(segment.points):
            position = field_position or tuple(segment.points[len(segment.points) // 2].astype(int))
            normal = np.asarray(segment.normal, dtype=float)
            plastic_distortion = np.outer(b, normal)
            strain_increment = 0.5 * (plastic_distortion + plastic_distortion.T)
            strain_increment += dq * domain.formation_volume * np.eye(2)
            self.full_field.add_event(position, strain_increment)
        for tj in self._boundary_to_tjs.get(segment.entity_id, ()):
            tj.add_burgers(b)
        ds_release = mode.delta_s * packet
        released = domain.shear.release(ds_release)
        q_release = domain.free_volume.accommodate(abs(dq) if dq else mode.delta_q * packet)
        self.ledger.write({
            "run_id": self.run_id, "time": self.solver.time if event_time is None else event_time,
            "step": self.solver.step_number,
            "temperature": self.config.pf.temperature, "seed": self.config.seed,
            "event_id": f"{domain.entity_id}:{domain.event_counter}", "event_type": event_type,
            "grain_ids": f"{segment.grain_i};{segment.grain_j}", "entity_id": domain.entity_id,
            "position": ledger_position if ledger_position is not None else (
                json.dumps(segment.points.mean(axis=0).tolist()) if len(segment.points) else ""
            ),
            "geometry_measure_Q": domain.encounter.total_measure,
            "grain_size": 0.5 * (self.snapshot.grains[segment.grain_i].equivalent_radius + self.snapshot.grains[segment.grain_j].equivalent_radius),
            "curvature": segment.curvature, "local_velocity": segment.velocity,
            "barrier_type": mode.family, "DeltaG0": mode.barrier_ev,
            "effective_DeltaG": mode.effective_barrier_ev(driving),
            "activation_volume": f"{mode.activation_volume_normal};{mode.activation_volume_shear}",
            "local_shear_stress": driving.resolved_shear,
            "local_normal_free_volume_stress": driving.normal_pressure,
            "shear_state_s": domain.shear.state, "free_volume_state_q": domain.free_volume.deficit,
            "Ns": mode.site_multiplicity, "nu0": mode.attempt_frequency, "instantaneous_rate": rate,
            "cumulative_hazard": domain.activation.clock.cumulative_hazard,
            "random_hazard_threshold": domain.activation.clock.threshold,
            "hit_count": domain.activation.hit_count, "required_hits_K": domain.hits,
            "release_Delta_s": ds_release, "release_Delta_q": q_release,
            "GB_area_change": delta_length, "TJ_travel": 0,
            "point_defect_quota": domain.free_volume.deficit,
            "normal_step_h": h, "burgers_vector_b": json.dumps(b.tolist()), "Nv": dq,
            "shear_strain_increment": shear_strain,
            "volumetric_strain_increment": volumetric_strain,
            "packet_size": packet, "Git_SHA": self.sha,
        })
        if released >= 0:
            domain.blocked = False

    def _stage_rates(self, domain: DomainPhysics,
                     segment: GBSegment) -> tuple[float, float, float]:
        p, temperature = self.config.parameters, self.config.pf.temperature
        def arrhenius(prefactor: float, barrier: float) -> float:
            return prefactor * np.exp(-barrier / (K_B_EV * temperature))
        nucleation = arrhenius(
            float(p.get("nucleation_prefactor", 1e4)),
            float(p.get("nucleation_barrier_ev", 0.45)),
        )
        exchange_current = arrhenius(
            float(p.get("exchange_prefactor", 1e4)),
            float(p.get("exchange_barrier_ev", 0.55)),
        )
        exchange_flux = abs(butler_volmer_flux(
            domain.free_volume.chemical_potential, temperature, exchange_current,
            float(p.get("exchange_transfer_coefficient", 0.5)),
        ))
        exchange = exchange_flux / max(
            domain.climb.required_quota, np.finfo(float).tiny
        )
        diffusion = diffusivity(
            temperature, float(p.get("transport_prefactor", 1e4)),
            float(p.get("transport_barrier_ev", 0.70)),
        )
        transport_length = float(p.get(
            "transport_length", max(segment.length, self.config.pf.grid_spacing)
        ))
        transport = 1.0 / transport_time(
            transport_length, diffusion,
            float(p.get("transport_geometry_factor", 1.0)),
        )
        return nucleation, exchange, transport

    def _record_climb_transition(self, domain: DomainPhysics, segment: GBSegment,
                                 event_time: float, stage: str, rate: float,
                                 threshold: float) -> None:
        domain.event_counter += 1
        self.ledger.write({
            "run_id": self.run_id, "time": event_time, "step": self.solver.step_number,
            "temperature": self.config.pf.temperature, "seed": self.config.seed,
            "event_id": f"{domain.entity_id}:{domain.event_counter}",
            "event_type": stage, "grain_ids": f"{segment.grain_i};{segment.grain_j}",
            "entity_id": domain.entity_id,
            "position": segment.points.mean(axis=0).tolist() if len(segment.points) else "",
            "geometry_measure_Q": domain.encounter.total_measure,
            "curvature": segment.curvature, "local_velocity": segment.velocity,
            "local_normal_free_volume_stress": domain.free_volume.chemical_potential,
            "free_volume_state_q": domain.free_volume.deficit,
            "instantaneous_rate": rate,
            "cumulative_hazard": threshold,
            "random_hazard_threshold": threshold,
            "hit_count": 1, "required_hits_K": 1,
            "release_Delta_s": 0.0, "release_Delta_q": 0.0,
            "GB_area_change": 0.0, "TJ_travel": 0.0,
            "point_defect_quota": domain.free_volume.deficit,
            "normal_step_h": 0.0, "burgers_vector_b": "", "Nv": 0.0,
            "shear_strain_increment": 0.0,
            "volumetric_strain_increment": 0.0,
            "packet_size": self.config.parameters.get("packet_size", 1.0),
            "Git_SHA": self.sha,
        })

    def _advance_climb(self, domain: DomainPhysics, segment: GBSegment, delta_length: float) -> None:
        modules = set(self.config.active_modules)
        if not modules.intersection({"nucleation_limited", "multihit_nucleation", "exchange_limited", "transport_limited", "serial_climb", "independent_and"}):
            return
        if domain.free_volume.deficit <= float(self.config.parameters.get("climb_trigger_quota", 0.25)):
            return
        domain.blocked = True
        rn, re, rt = self._stage_rates(domain, segment)
        complete = False
        event_time = None
        if modules.intersection({"serial_climb", "independent_and"}):
            if domain.climb.stage.value in {"inactive", "quota_completion"}:
                domain.climb.activate(self.solver.time)
            complete = domain.climb.advance(self.config.pf.time_step,
                self.solver.time - self.config.pf.time_step, rn, re, rt)
            event_time = domain.climb.last_completion_time
            transition_rates = {
                "exchange": ("climb_nucleation", rn),
                "transport": ("climb_exchange", re),
                "quota_completion": ("climb_transport", rt),
            }
            for transition_time, transition_stage, threshold in domain.climb.last_transitions:
                name, rate = transition_rates[transition_stage.value]
                self._record_climb_transition(
                    domain, segment, transition_time, name, rate, threshold
                )
        else:
            rate = rn if modules.intersection({"nucleation_limited", "multihit_nucleation"}) else (re if "exchange_limited" in modules else rt)
            completions, hits = self._advance_activation(
                domain, rate, self.config.pf.time_step,
                self.solver.time - self.config.pf.time_step,
                stop_after_completion=True,
            )
            self._record_activation_hits(domain, rate, hits, segment=segment)
            complete = bool(completions)
            event_time = completions[0].time if completions else None
        if complete:
            release = float(self.config.parameters.get("climb_release_quota", 1.0))
            domain.free_volume.accommodate(release)
            domain.blocked = domain.free_volume.deficit > float(self.config.parameters.get("climb_trigger_quota", 0.25))
            mode, total, driving = self._activation_mode(domain, segment)
            self._record_event(domain, segment, mode, total, driving,
                               "climb_quota_completion", delta_length, event_time)
            domain.blocked = domain.free_volume.deficit > float(self.config.parameters.get("climb_trigger_quota", 0.25))

    def _update_tj_physics(self, mobility: np.ndarray) -> None:
        modules = set(self.config.active_modules)
        enabled = bool(modules.intersection({"tj_compatibility", "tj_pinning", "tj_burgers_strict", "tj_burgers_residual", "tj_geometric_surrogate"}))
        if not enabled:
            self.tj_domains.clear()
            return
        current = set(self.snapshot.triple_junctions)
        self.tj_domains = {key: value for key, value in self.tj_domains.items() if key in current}
        for key, tj in self.snapshot.triple_junctions.items():
            if key not in self.tj_domains:
                fake = GBSegment(tj.grain_ids[0], tj.grain_ids[1], 0)
                fake.points = np.asarray([tj.position]); fake.length = self.config.pf.grid_spacing
                fake_key = fake.entity_id
                fake.segment_id = abs(sum(tj.grain_ids))
                domain = self._new_domain(fake)
                domain.entity_id = key
                self.tj_domains[key] = domain
            domain = self.tj_domains[key]
            delta_path = max(0.0, tj.travel_distance - domain.previous_length)
            domain.previous_length = tj.travel_distance
            explicit_residual = bool(
                modules.intersection({"tj_burgers_strict", "tj_burgers_residual"})
                and np.linalg.norm(tj.residual_burgers) > 1e-10
            )
            if self.solver.step_number > 0 and not domain.blocked:
                encounter = (
                    [] if explicit_residual else
                    domain.encounter.advance(delta_path, maximum_events=1)
                )
                if explicit_residual or encounter:
                    domain.blocked = True
                    self._begin_activation_window(domain)
            if domain.blocked:
                rate = float(self.config.parameters.get("tj_attempt_frequency", 1e3)) * np.exp(
                    -float(self.config.parameters.get("tj_barrier_ev", 0.6)) / (K_B_EV * self.config.pf.temperature))
                completions, hits = self._advance_activation(
                    domain, rate, self.config.pf.time_step,
                    self.solver.time - self.config.pf.time_step,
                    stop_after_completion=True,
                )
                self._record_activation_hits(
                    domain, rate, hits,
                    grain_ids=";".join(map(str, tj.grain_ids)),
                    position=tj.position, geometry_measure=tj.travel_distance,
                    tj_travel=delta_path,
                )
                if completions:
                    domain.event_counter += 1
                    if explicit_residual:
                        target = -tj.residual_burgers
                        mode = min(self.modes, key=lambda m: np.linalg.norm(np.asarray(m.burgers) - target))
                        tj.add_burgers(np.asarray(mode.burgers))
                    domain.blocked = bool("tj_burgers_strict" in modules and np.linalg.norm(tj.residual_burgers) > 1e-10)
                    self.ledger.write({
                        "run_id": self.run_id, "time": completions[0].time,
                        "step": self.solver.step_number,
                        "temperature": self.config.pf.temperature, "seed": self.config.seed,
                        "event_id": f"{key}:{domain.event_counter}", "event_type": "tj_compatibility_release",
                        "grain_ids": ";".join(map(str, tj.grain_ids)), "entity_id": key,
                        "position": tj.position, "geometry_measure_Q": tj.travel_distance,
                        "TJ_travel": delta_path, "instantaneous_rate": rate,
                        "cumulative_hazard": domain.activation.clock.cumulative_hazard,
                        "random_hazard_threshold": domain.activation.clock.threshold,
                        "hit_count": domain.activation.hit_count, "required_hits_K": domain.hits,
                        "burgers_vector_b": tj.residual_burgers.tolist(), "Git_SHA": self.sha,
                    })
            if domain.blocked:
                y, x = np.rint(tj.position).astype(int) % np.asarray(self.config.pf.shape)
                radius = int(self.config.parameters.get("tj_correlation_radius", 2))
                for oy in range(-radius, radius + 1):
                    for ox in range(-radius, radius + 1):
                        mobility[(y + oy) % mobility.shape[0], (x + ox) % mobility.shape[1]] = 0.0

    def _update_physics(self) -> None:
        cfg, modules = self.config, set(self.config.active_modules)
        self._boundary_to_tjs = self._index_boundary_tjs()
        mobility = np.ones(cfg.pf.shape)
        self.driving_field.fill(0.0)
        entity_elapsed = self.solver.time - self.previous_entity_time
        current_ids = set(self.snapshot.boundaries)
        self.domains = {key: state for key, state in self.domains.items() if key in current_ids}
        for key, segment in self.snapshot.boundaries.items():
            if key not in self.domains:
                self.domains[key] = self._new_domain(segment)
            domain = self.domains[key]
            delta_length = 0.0 if not domain.previous_length else abs(segment.length - domain.previous_length)
            domain.previous_length = segment.length
            grain_i = self.snapshot.grains[segment.grain_i]
            grain_j = self.snapshot.grains[segment.grain_j]
            elapsed = self.solver.time - domain.previous_time
            if domain.previous_area_i and domain.previous_area_j and elapsed > 0:
                displacement_i = (grain_i.area - domain.previous_area_i) / max(grain_i.perimeter, cfg.pf.grid_spacing)
                displacement_j = (grain_j.area - domain.previous_area_j) / max(grain_j.perimeter, cfg.pf.grid_spacing)
                normal_displacement = 0.5 * (displacement_i - displacement_j)
                segment.velocity = normal_displacement / elapsed
            else:
                normal_displacement = 0.0
                segment.velocity = 0.0
            if domain.normal_release_remaining * normal_displacement > 0.0:
                consumed = min(
                    abs(domain.normal_release_remaining), abs(normal_displacement)
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
            if (domain.normal_release_remaining == 0.0
                    and abs(domain.normal_displacement_ledger) >= release_distance):
                direction = np.sign(domain.normal_displacement_ledger)
                domain.normal_release_remaining = direction * release_distance
                domain.normal_displacement_ledger -= direction * release_distance
            domain.previous_area_i = grain_i.area
            domain.previous_area_j = grain_j.area
            domain.previous_time = self.solver.time
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
            ci = np.asarray(self.snapshot.grains[segment.grain_i].centroid)
            cj = np.asarray(self.snapshot.grains[segment.grain_j].centroid)
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
                domain.shear.migrate(beta, normal_displacement, cfg.pf.time_step)
            if "qiu_reference_shear" in modules and self.full_field is not None and normal_displacement:
                beta = float(cfg.parameters.get("easy_beta", 0.35))
                tangent = np.asarray((-normal[1], normal[0]))
                displacement = beta * normal_displacement * tangent
                strain = 0.5 * (np.outer(displacement, normal) + np.outer(normal, displacement))
                position = tuple(segment.points[len(segment.points) // 2].astype(int))
                self.full_field.add_event(position, strain)
            if modules.intersection({"free_volume", "serial_climb", "nucleation_limited", "multihit_nucleation", "exchange_limited", "transport_limited", "mixed_shear_climb_event", "independent_and"}):
                domain.free_volume.require_for_area_change(delta_length)
            compatibility_trigger = (
                abs(domain.shear.state) >= float(cfg.parameters.get("shear_trigger", 0.25))
                and domain.free_volume.deficit >= float(cfg.parameters.get("climb_trigger_quota", 0.25))
            )
            if (modules.intersection({"mixed_shear_climb_event", "independent_and"})
                    and compatibility_trigger and not domain.blocked):
                domain.blocked = True
                self._begin_activation_window(domain)
                domain.climb.activate(self.solver.time)

            encounter_enabled = cfg.compatibility_model == "geometric_surrogate" or bool(
                modules.intersection({
                    "gb_compatibility", "gb_area_point_defect_pinning", "gb_pinning"
                })
            )
            if (self.solver.step_number > 0 and encounter_enabled and not domain.blocked
                    and domain.encounter.advance(delta_length, maximum_events=1)):
                domain.blocked = True
                self._begin_activation_window(domain)
                domain.climb.activate(self.solver.time)

            if domain.blocked and self.solver.step_number > 0:
                candidates, rates, normal_pressure, resolved_shear, vacancy_mu = self._activation_rates(
                    domain, segment
                )
                total_rate = float(rates.sum())
                completions, hits = self._advance_activation(
                    domain, total_rate, cfg.pf.time_step,
                    self.solver.time - cfg.pf.time_step,
                    stop_after_completion=True,
                )
                self._record_activation_hits(
                    domain, total_rate, hits, segment=segment, position=ledger_position
                )
                if completions:
                    if total_rate:
                        selected = int(domain.rng.choice(len(candidates), p=rates / total_rate))
                    else:
                        selected = 0
                    mode = candidates[selected]
                    driving = self._mode_driving(
                        normal_pressure, resolved_shear, vacancy_mu, selected
                    )
                    self._record_event(domain, segment, mode, total_rate, driving,
                                       "compatibility_release", delta_length, completions[0].time,
                                       ledger_position=ledger_position,
                                       field_position=field_position)
            elif (cfg.compatibility_model == "explicit_modes" and not encounter_enabled
                  and self.solver.step_number > 0):
                # Event-resolved easy-mode flux; completion changes finite hidden
                # kinematics even when it does not gate mobility.
                candidates, rates, normal_pressure, resolved_shear, vacancy_mu = self._activation_rates(
                    domain, segment
                )
                total_rate = float(rates.sum())
                completions, hits = self._advance_activation(
                    domain, total_rate, cfg.pf.time_step,
                    self.solver.time - cfg.pf.time_step,
                )
                self._record_activation_hits(
                    domain, total_rate, hits, segment=segment, position=ledger_position
                )
                for completion in completions:
                    if total_rate:
                        selected = int(domain.rng.choice(len(candidates), p=rates / total_rate))
                    else:
                        selected = 0
                    mode = candidates[selected]
                    driving = self._mode_driving(
                        normal_pressure, resolved_shear, vacancy_mu, selected
                    )
                    self._record_event(domain, segment, mode, total_rate, driving,
                                       "disconnection_mode", delta_length, completion.time,
                                       ledger_position=ledger_position,
                                       field_position=field_position)

            self._advance_climb(domain, segment, delta_length)

            pair_force = float(cfg.parameters.get("easy_beta", 0.35)) * self._boundary_resolved_shear(domain, segment)
            if domain.normal_release_remaining:
                pair_force += np.sign(domain.normal_release_remaining) * float(
                    cfg.parameters.get("event_normal_pressure", 1.0)
                )
            for yx in segment.points.astype(int):
                y, x = yx % np.asarray(cfg.pf.shape)
                if domain.blocked:
                    mobility[y, x] = min(mobility[y, x], float(cfg.parameters.get("pinned_mobility_fraction", 0.0)))
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
        self.solver.set_mobility_scale(mobility)
        if self.full_field is not None:
            self.full_field.solve()
        self.previous_entity_eta = self.solver.eta.copy()
        self.previous_entity_time = self.solver.time

    def _write_tracks(self) -> None:
        for grain in self.snapshot.grains.values():
            self.track_writer.writerow({
                "run_id": self.run_id, "time": self.solver.time, "step": self.solver.step_number,
                "grain_id": grain.grain_id, "area": grain.area, "radius": grain.equivalent_radius,
                "perimeter": grain.perimeter, "neighbors": len(grain.neighbors),
            })
        self.track_handle.flush()
        for segment in self.snapshot.boundaries.values():
            domain = self.domains.get(segment.entity_id)
            if domain is None:
                continue
            self.boundary_writer.writerow({
                "run_id": self.run_id, "time": self.solver.time, "step": self.solver.step_number,
                "entity_id": segment.entity_id, "grain_i": segment.grain_i, "grain_j": segment.grain_j,
                "length": segment.length, "curvature": segment.curvature,
                "normal_velocity": segment.velocity, "blocked": int(domain.blocked),
                "resolved_shear": self._boundary_resolved_shear(domain, segment),
                "free_volume_deficit": domain.free_volume.deficit,
            })
        self.boundary_handle.flush()

    def _save_checkpoint(self) -> None:
        event_ledger_offset = self.ledger.checkpoint()
        stream_offsets = {}
        for name, handle in (
            ("grain_tracks_offset", self.track_handle),
            ("boundary_tracks_offset", self.boundary_handle),
        ):
            handle.flush()
            os.fsync(handle.fileno())
            stream_offsets[name] = handle.tell()
        state = {
            "time": self.solver.time, "step_number": self.solver.step_number,
            "domains": {key: domain.state_dict() for key, domain in self.domains.items()},
            "tj_domains": {key: domain.state_dict() for key, domain in self.tj_domains.items()},
            "energy_records": self.energy_records,
            "accumulated_shear_strain": self.accumulated_shear_strain,
            "accumulated_volumetric_strain": self.accumulated_volumetric_strain,
            "previous_entity_time": self.previous_entity_time,
            "event_ledger_offset": event_ledger_offset,
            **stream_offsets,
        }
        serialized_state = json.dumps(state, indent=2) + "\n"
        arrays: dict[str, np.ndarray] = {
            "eta": self.solver.eta, "mobility_scale": self.solver.mobility_scale,
            "driving_field": self.driving_field, "active_phases": self.solver.active_phases,
            "orientations": self.orientations,
            "previous_entity_eta": self.previous_entity_eta,
            # The archive is the authoritative checkpoint generation.  Keeping
            # its metadata here prevents an interruption between the two atomic
            # replacements from pairing new arrays with stale JSON.
            "checkpoint_state_json": np.asarray(serialized_state),
        }
        if self.full_field is not None:
            arrays["eigenstrain"] = self.full_field.eigenstrain
        atomic_savez_compressed(self.output_dir / "checkpoint.npz", **arrays)
        atomic_write_text(self.output_dir / "checkpoint.json", serialized_state)

    def _load_checkpoint(self) -> None:
        with np.load(self.output_dir / "checkpoint.npz") as arrays:
            if "checkpoint_state_json" in arrays:
                state = json.loads(str(arrays["checkpoint_state_json"]))
            else:
                # Backward compatibility for production runs created before
                # checkpoint metadata was embedded in the archive.
                state = json.loads((self.output_dir / "checkpoint.json").read_text())
            if "event_ledger_offset" in state:
                self.ledger.truncate(int(state["event_ledger_offset"]))
            for name, handle in (
                ("grain_tracks_offset", self.track_handle),
                ("boundary_tracks_offset", self.boundary_handle),
            ):
                if name in state:
                    offset = int(state[name])
                    handle.flush()
                    size = Path(handle.name).stat().st_size
                    if offset < 0 or offset > size:
                        raise ValueError(
                            f"invalid {name} checkpoint offset {offset} for {size} bytes"
                        )
                    handle.seek(offset)
                    handle.truncate()
                    handle.flush()
                    os.fsync(handle.fileno())
                    handle.seek(0, 2)
            self.solver.eta = arrays["eta"].copy()
            self.solver.mobility_scale = arrays["mobility_scale"].copy()
            self.solver.active_phases = arrays["active_phases"].astype(bool).copy()
            if "orientations" in arrays:
                self.orientations = arrays["orientations"].copy()
            self.driving_field = arrays["driving_field"].copy()
            self.previous_entity_eta = (
                arrays["previous_entity_eta"].copy()
                if "previous_entity_eta" in arrays else arrays["eta"].copy()
            )
            if self.full_field is not None and "eigenstrain" in arrays:
                self.full_field.eigenstrain = arrays["eigenstrain"].copy()
        self.solver.time = float(state["time"])
        self.solver.step_number = int(state["step_number"])
        self.previous_entity_time = float(state.get("previous_entity_time", self.solver.time))
        self.tracker = EntityTracker(
            self.orientations, self.config.pf.grid_spacing,
            float(self.config.parameters.get("event_domain_length", 8.0)),
            periodic=self.config.pf.boundary_conditions == "periodic",
        )
        self.snapshot = self.tracker.update(self.solver.labels)
        for key, domain_state in state["domains"].items():
            if key in self.snapshot.boundaries:
                domain = self._new_domain(self.snapshot.boundaries[key])
                domain.load_state_dict(domain_state)
                self.domains[key] = domain
        for key, domain_state in state.get("tj_domains", {}).items():
            if key in self.snapshot.triple_junctions:
                tj = self.snapshot.triple_junctions[key]
                fake = GBSegment(tj.grain_ids[0], tj.grain_ids[1], 0)
                fake.points = np.asarray([tj.position]); fake.length = self.config.pf.grid_spacing
                domain = self._new_domain(fake); domain.entity_id = key
                domain.load_state_dict(domain_state); self.tj_domains[key] = domain
        self.energy_records = state["energy_records"]
        self.accumulated_shear_strain = float(state.get("accumulated_shear_strain", 0.0))
        self.accumulated_volumetric_strain = float(state.get("accumulated_volumetric_strain", 0.0))

    def run(self) -> Path:
        failure: str | None = None
        try:
            entity_every_step = bool(self.config.active_modules) or self.config.compatibility_model != "off"
            for _ in range(max(0, self.config.max_steps - self.solver.step_number)):
                diag = self.solver.step()
                update_entities = entity_every_step or self.solver.step_number % self.config.output_cadence == 0
                if update_entities:
                    self.snapshot = self.tracker.update(self.solver.labels)
                    self._update_physics()
                stored_shear = sum(d.shear.energy for d in self.domains.values())
                stored_free_volume = sum(d.free_volume.energy for d in self.domains.values())
                dissipated_shear = sum(d.shear.dissipated_energy for d in self.domains.values())
                dissipated_free_volume = sum(
                    d.free_volume.dissipated_energy for d in self.domains.values()
                )
                self.energy_records.append({
                    "time": diag.time,
                    "interfacial": diag.interfacial_energy,
                    "stored": stored_shear + stored_free_volume,
                    "stored_shear": stored_shear,
                    "stored_free_volume": stored_free_volume,
                    "dissipated_shear": dissipated_shear,
                    "dissipated_free_volume": dissipated_free_volume,
                })
                if self.solver.step_number % self.config.output_cadence == 0:
                    self._write_tracks()
                    self._save_checkpoint()
                if update_entities and len(self.snapshot.grains) <= self.config.termination_grains:
                    break
        except BaseException as exc:
            failure = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            if (failure is None and self.solver.step_number > 0
                    and self.solver.step_number % self.config.output_cadence != 0):
                self._write_tracks()
                self._save_checkpoint()
            self.ledger.close()
            self.track_handle.close()
            self.boundary_handle.close()
            (self.output_dir / "energy.json").write_text(json.dumps(self.energy_records, indent=2) + "\n")
            restart_artifacts = []
            for name in ("checkpoint.npz", "checkpoint.json"):
                artifact = self.output_dir / name
                if artifact.exists():
                    restart_artifacts.append({
                        "path": str(artifact),
                        "sha256": file_sha256(artifact),
                        "size_bytes": artifact.stat().st_size,
                    })
            write_manifest(self.output_dir / "manifest.json", self.config.to_dict(),
                           "failed" if failure else "completed", {
                               "failure": failure, "steps_completed": self.solver.step_number,
                               "final_grains": len(self.snapshot.grains),
                               "accumulated_shear_strain": self.accumulated_shear_strain,
                               "accumulated_volumetric_strain": self.accumulated_volumetric_strain,
                               "event_ledger": str(self.ledger.path),
                               "restart_artifacts": restart_artifacts,
                           }, code_sha=self.sha)
        return self.output_dir
