from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from grain_growth_pf.climb.free_volume import FreeVolumeState
from grain_growth_pf.climb.serial_cycle import SerialClimbCycle
from grain_growth_pf.config import ModelConfig
from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving
from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.disconnections.shear_coupling import event_shear_increment, event_volumetric_increment
from grain_growth_pf.disconnections.spectrum import isotropic_surrogate_library
from grain_growth_pf.encounters.geometric_hazard import GeometricEncounterClock
from grain_growth_pf.entities.gb_segment import GBSegment
from grain_growth_pf.entities.tracker import EntityTracker
from grain_growth_pf.io.event_ledger import EventLedger
from grain_growth_pf.io.provenance import git_sha, write_manifest
from grain_growth_pf.mechanics.local_shear_memory import LocalShearMemory
from grain_growth_pf.mechanics.qiu_full_field import QiuFullField
from grain_growth_pf.obstacles.particles import ParticleField
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from grain_growth_pf.stochastic.multihit import MultiHitProcess


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
                           "hit_count": self.activation.hit_count},
            "shear": {"state": self.shear.state, "dissipated_energy": self.shear.dissipated_energy},
            "free_volume": {"required_total": self.free_volume.required_total,
                            "accommodated_total": self.free_volume.accommodated_total},
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
        self.shear.state = state["shear"]["state"]
        self.shear.dissipated_energy = state["shear"]["dissipated_energy"]
        self.free_volume.required_total = state["free_volume"]["required_total"]
        self.free_volume.accommodated_total = state["free_volume"]["accommodated_total"]
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
        self.event_counter = state["event_counter"]


def _child_rng(seed: int, entity_id: str) -> np.random.Generator:
    # Stable hierarchy independent of dictionary/Python hash ordering.
    words = np.frombuffer(entity_id.encode(), dtype=np.uint8).astype(np.uint32)
    sequence = np.random.SeedSequence([seed, *words.tolist()])
    return np.random.default_rng(sequence)


class EventResolvedSimulation:
    """Couples persistent entity clocks and compatibility state to the PF solver."""

    def __init__(self, config: ModelConfig, output_dir: str | Path, resume: bool = False):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=resume)
        self.sha = git_sha()
        eta, seeds, orientations = voronoi_polycrystal(
            config.pf.shape, int(config.parameters.get("initial_grains", 50)), config.seed,
            width=config.pf.interface_width / 2,
            periodic=config.pf.boundary_conditions == "periodic",
        )
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
        equilibration_steps = int(config.parameters.get("equilibration_steps", 0))
        active_original_ids = np.arange(len(orientations))
        if not resume:
            write_manifest(self.output_dir / "manifest.json", config.to_dict(), "equilibrating", {
                "initial_seed_positions": seeds.tolist(),
                "original_orientations": orientations.tolist(),
            }, code_sha=self.sha)
            for _ in range(equilibration_steps):
                self.solver.step()
            target_grains = config.parameters.get("equilibrate_to_grains")
            if target_grains is not None:
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
            if bool(config.parameters.get("compact_after_equilibration", True)):
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
        )
        self.full_field = QiuFullField(config.pf.shape) if config.mechanics_backend == "qiu_full_field" else None
        particle_modules = {"random_spatial_pinning", "particle_zener"}.intersection(config.active_modules)
        self.particles = ParticleField.random(
            int(config.parameters.get("particle_count", 20)),
            float(config.parameters.get("particle_radius", 1.5)), config.pf.shape,
            config.seed + 991,
        ) if particle_modules else None
        self.ledger = EventLedger(self.output_dir / "events.csv")
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
            }, code_sha=self.sha)
            self._write_tracks()

    def _driving(self, _eta: np.ndarray, _time: float) -> np.ndarray:
        return self.driving_field

    def _new_domain(self, segment: GBSegment) -> DomainPhysics:
        p = self.config.parameters
        modules = set(self.config.active_modules)
        hits = int(p.get("required_hits", 3 if any("multihit" in m for m in modules) else 1))
        interpretation = "packet_reset" if "multihit_packet_reset" in modules else "persistent_hits"
        return DomainPhysics(
            segment.entity_id, _child_rng(self.config.seed, segment.entity_id),
            float(p.get("encounter_density", 0.08)), hits, interpretation,
            float(p.get("shear_stiffness", 0.15)),
            (float(p.get("shear_relaxation_time", 1.0)) if "deterministic_relaxation" in modules else None),
            float(p.get("excess_volume_per_area", 0.02)),
            float(p.get("point_defect_formation_volume", 0.01)), float(p.get("free_volume_stiffness", 0.05)),
        )

    def _activation_rates(self, domain: DomainPhysics, segment: GBSegment) -> tuple[list[DisconnectionMode], np.ndarray, list[ModeDriving]]:
        capillary = self.config.pf.gb_energy * segment.curvature
        candidates = [m for m in self.modes if (m.family != "easy" if domain.blocked else True)]
        if "mixed_shear_climb_event" in self.config.active_modules:
            candidates = [m for m in candidates if m.delta_s > 0 and m.delta_q > 0]
        normal = np.asarray(segment.normal, dtype=float)
        tangent = np.asarray((-normal[1], normal[0]))
        stress = None
        if self.full_field is not None and len(segment.points):
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            y, x = position[0] % self.config.pf.shape[0], position[1] % self.config.pf.shape[1]
            stress = self.full_field.stress[:, :, y, x]
        drivings = []
        for mode in candidates:
            b = np.asarray(mode.burgers, dtype=float)
            b_direction = b / max(np.linalg.norm(b), np.finfo(float).tiny)
            resolved_shear = domain.shear.internal_shear_stress * float(b_direction @ tangent)
            if stress is not None:
                resolved_shear += mode.resolved_shear(stress, normal)
            drivings.append(ModeDriving(
                normal_pressure=capillary + domain.free_volume.chemical_potential,
                resolved_shear=resolved_shear,
                vacancy_chemical_potential=domain.free_volume.chemical_potential,
            ))
        rates = np.asarray([
            mode.rate(self.config.pf.temperature, driving)
            for mode, driving in zip(candidates, drivings)
        ])
        return candidates, rates, drivings

    def _activation_mode(self, domain: DomainPhysics, segment: GBSegment) -> tuple[DisconnectionMode, float, ModeDriving]:
        candidates, rates, drivings = self._activation_rates(domain, segment)
        total = float(rates.sum())
        if total == 0:
            return candidates[0], 0.0, drivings[0]
        mode_rng = domain.rng
        selected = int(mode_rng.choice(len(candidates), p=rates / total))
        return candidates[selected], total, drivings[selected]

    def _boundary_resolved_shear(self, domain: DomainPhysics, segment: GBSegment) -> float:
        value = domain.shear.internal_shear_stress
        if self.full_field is not None and len(segment.points):
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            normal = np.asarray(segment.normal, dtype=float)
            tangent = np.asarray((-normal[1], normal[0]))
            value += self.full_field.resolved_shear(position, tangent, normal)
        return float(value)

    def _record_event(self, domain: DomainPhysics, segment: GBSegment, mode: DisconnectionMode,
                      rate: float, driving: ModeDriving, event_type: str, delta_length: float,
                      event_time: float | None = None) -> None:
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
            position = tuple(segment.points[len(segment.points) // 2].astype(int))
            normal = np.asarray(segment.normal, dtype=float)
            plastic_distortion = np.outer(b, normal)
            strain_increment = 0.5 * (plastic_distortion + plastic_distortion.T)
            strain_increment += dq * domain.formation_volume * np.eye(2)
            self.full_field.add_event(position, strain_increment)
        for tj in self.snapshot.triple_junctions.values():
            if segment.entity_id in tj.adjoining_boundaries:
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
            "position": segment.points.mean(axis=0).tolist() if len(segment.points) else "",
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
            "normal_step_h": h, "burgers_vector_b": b.tolist(), "Nv": dq,
            "shear_strain_increment": shear_strain,
            "volumetric_strain_increment": volumetric_strain,
            "packet_size": packet, "Git_SHA": self.sha,
        })
        if released >= 0:
            domain.blocked = False

    def _stage_rates(self) -> tuple[float, float, float]:
        p, temperature = self.config.parameters, self.config.pf.temperature
        def arrhenius(prefactor: float, barrier: float) -> float:
            return prefactor * np.exp(-barrier / (K_B_EV * temperature))
        return (
            arrhenius(float(p.get("nucleation_prefactor", 1e4)), float(p.get("nucleation_barrier_ev", 0.45))),
            arrhenius(float(p.get("exchange_prefactor", 1e4)), float(p.get("exchange_barrier_ev", 0.55))),
            arrhenius(float(p.get("transport_prefactor", 1e4)), float(p.get("transport_barrier_ev", 0.70))),
        )

    def _advance_climb(self, domain: DomainPhysics, segment: GBSegment, delta_length: float) -> None:
        modules = set(self.config.active_modules)
        if not modules.intersection({"nucleation_limited", "multihit_nucleation", "exchange_limited", "transport_limited", "serial_climb", "independent_and"}):
            return
        if domain.free_volume.deficit <= float(self.config.parameters.get("climb_trigger_quota", 0.25)):
            return
        domain.blocked = True
        rn, re, rt = self._stage_rates()
        complete = False
        event_time = None
        if modules.intersection({"serial_climb", "independent_and"}):
            if domain.climb.stage.value in {"inactive", "quota_completion"}:
                domain.climb.activate(self.solver.time)
            complete = domain.climb.advance(self.config.pf.time_step,
                self.solver.time - self.config.pf.time_step, rn, re, rt)
            event_time = domain.climb.last_completion_time
        else:
            rate = rn if modules.intersection({"nucleation_limited", "multihit_nucleation"}) else (re if "exchange_limited" in modules else rt)
            completions = domain.activation.advance(rate, self.config.pf.time_step,
                                                     self.solver.time - self.config.pf.time_step)
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
            explicit_residual = modules.intersection({"tj_burgers_strict", "tj_burgers_residual"}) and np.linalg.norm(tj.residual_burgers) > 1e-10
            if self.solver.step_number > 0 and not domain.blocked:
                if domain.encounter.advance(delta_path) or explicit_residual:
                    domain.blocked = True
                    domain.activation.begin_window()
            if domain.blocked:
                rate = float(self.config.parameters.get("tj_attempt_frequency", 1e3)) * np.exp(
                    -float(self.config.parameters.get("tj_barrier_ev", 0.6)) / (K_B_EV * self.config.pf.temperature))
                completions = domain.activation.advance(
                    rate, self.config.pf.time_step, self.solver.time - self.config.pf.time_step
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
        mobility = np.ones(cfg.pf.shape)
        self.driving_field.fill(0.0)
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
            domain.previous_area_i = grain_i.area
            domain.previous_area_j = grain_j.area
            domain.previous_time = self.solver.time
            ri, rj = grain_i.equivalent_radius, grain_j.equivalent_radius
            segment.curvature = 0.5 * (1.0 / max(rj, cfg.pf.grid_spacing) - 1.0 / max(ri, cfg.pf.grid_spacing))
            ci = np.asarray(self.snapshot.grains[segment.grain_i].centroid)
            cj = np.asarray(self.snapshot.grains[segment.grain_j].centroid)
            normal = cj - ci
            if cfg.pf.boundary_conditions == "periodic":
                box = np.asarray(cfg.pf.shape, dtype=float)
                normal -= np.round(normal / box) * box
            normal /= max(np.linalg.norm(normal), np.finfo(float).tiny)
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
                domain.activation.begin_window()
                domain.climb.activate(self.solver.time)

            encounter_enabled = cfg.compatibility_model == "geometric_surrogate" or bool(
                modules.intersection({"gb_area_point_defect_pinning", "gb_pinning"})
            )
            if (self.solver.step_number > 0 and encounter_enabled and not domain.blocked
                    and domain.encounter.advance(delta_length)):
                domain.blocked = True
                domain.activation.begin_window()
                domain.climb.activate(self.solver.time)

            candidates, rates, drivings = self._activation_rates(domain, segment)
            total_rate = float(rates.sum())
            if domain.blocked and self.solver.step_number > 0:
                completions = domain.activation.advance(total_rate, cfg.pf.time_step, self.solver.time - cfg.pf.time_step)
                if completions:
                    if total_rate:
                        selected = int(domain.rng.choice(len(candidates), p=rates / total_rate))
                    else:
                        selected = 0
                    mode, driving = candidates[selected], drivings[selected]
                    self._record_event(domain, segment, mode, total_rate, driving,
                                       "compatibility_release", delta_length, completions[0].time)
            elif cfg.compatibility_model == "explicit_modes" and self.solver.step_number > 0:
                # Event-resolved easy-mode flux; completion changes finite hidden
                # kinematics even when it does not gate mobility.
                completions = domain.activation.advance(total_rate, cfg.pf.time_step, self.solver.time - cfg.pf.time_step)
                for completion in completions:
                    if total_rate:
                        selected = int(domain.rng.choice(len(candidates), p=rates / total_rate))
                    else:
                        selected = 0
                    mode, driving = candidates[selected], drivings[selected]
                    self._record_event(domain, segment, mode, total_rate, driving,
                                       "disconnection_mode", delta_length, completion.time)

            self._advance_climb(domain, segment, delta_length)

            pair_force = float(cfg.parameters.get("easy_beta", 0.35)) * self._boundary_resolved_shear(domain, segment)
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
        arrays: dict[str, np.ndarray] = {
            "eta": self.solver.eta, "mobility_scale": self.solver.mobility_scale,
            "driving_field": self.driving_field, "active_phases": self.solver.active_phases,
            "orientations": self.orientations,
        }
        if self.full_field is not None:
            arrays["eigenstrain"] = self.full_field.eigenstrain
        np.savez_compressed(self.output_dir / "checkpoint.npz", **arrays)
        state = {
            "time": self.solver.time, "step_number": self.solver.step_number,
            "domains": {key: domain.state_dict() for key, domain in self.domains.items()},
            "tj_domains": {key: domain.state_dict() for key, domain in self.tj_domains.items()},
            "energy_records": self.energy_records,
            "accumulated_shear_strain": self.accumulated_shear_strain,
            "accumulated_volumetric_strain": self.accumulated_volumetric_strain,
        }
        (self.output_dir / "checkpoint.json").write_text(json.dumps(state, indent=2) + "\n")

    def _load_checkpoint(self) -> None:
        state = json.loads((self.output_dir / "checkpoint.json").read_text())
        with np.load(self.output_dir / "checkpoint.npz") as arrays:
            self.solver.eta = arrays["eta"].copy()
            self.solver.mobility_scale = arrays["mobility_scale"].copy()
            self.solver.active_phases = arrays["active_phases"].astype(bool).copy()
            if "orientations" in arrays:
                self.orientations = arrays["orientations"].copy()
            self.driving_field = arrays["driving_field"].copy()
            if self.full_field is not None and "eigenstrain" in arrays:
                self.full_field.eigenstrain = arrays["eigenstrain"].copy()
        self.solver.time = float(state["time"])
        self.solver.step_number = int(state["step_number"])
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
                stored = sum(d.shear.energy + d.free_volume.energy for d in self.domains.values())
                self.energy_records.append({"time": diag.time, "interfacial": diag.interfacial_energy, "stored": stored})
                if self.solver.step_number % self.config.output_cadence == 0:
                    self._write_tracks()
                    self._save_checkpoint()
                if update_entities and len(self.snapshot.grains) <= self.config.termination_grains:
                    break
        except Exception as exc:
            failure = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            self.ledger.close()
            self.track_handle.close()
            self.boundary_handle.close()
            (self.output_dir / "energy.json").write_text(json.dumps(self.energy_records, indent=2) + "\n")
            write_manifest(self.output_dir / "manifest.json", self.config.to_dict(),
                           "failed" if failure else "completed", {
                               "failure": failure, "steps_completed": self.solver.step_number,
                               "final_grains": len(self.snapshot.grains),
                               "accumulated_shear_strain": self.accumulated_shear_strain,
                               "accumulated_volumetric_strain": self.accumulated_volumetric_strain,
                           }, code_sha=self.sha)
        return self.output_dir
