import csv
import json

import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.simulation import EventResolvedSimulation


def test_event_resolved_smoke_writes_reproducible_schema(tmp_path):
    config = ModelConfig(
        regime="G1", seed=17,
        pf=PFConfig(shape=(24, 24), interface_width=3, time_step=0.02,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "single_hit_poisson", "shear_memory", "free_volume"),
        output_cadence=2, max_steps=6, termination_grains=1,
        parameters={"initial_grains": 6, "equilibration_steps": 0, "encounter_density": 5.0,
                    "attempt_frequency": 100.0, "event_domain_length": 100.0},
    )
    output = tmp_path / "run"
    EventResolvedSimulation(config, output).run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["config"]["seed"] == 17
    with (output / "grain_tracks.csv").open() as handle:
        tracks = list(csv.DictReader(handle))
    with (output / "events.csv").open() as handle:
        fields = csv.DictReader(handle).fieldnames
    assert tracks
    assert "random_hazard_threshold" in fields
    assert "burgers_vector_b" in fields


def test_qiu_full_field_backend_smoke(tmp_path):
    config = ModelConfig(
        regime="E1", seed=9,
        pf=PFConfig(shape=(20, 20), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="qiu_full_field", compatibility_model="explicit_modes",
        active_modules=("event_modes", "shear_feedback"), output_cadence=1,
        max_steps=2, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 0,
                    "attempt_frequency": 100.0, "barrier_core_ev": 0.0,
                    "b_coefficient_ev": 0.0, "h_coefficient_ev": 0.0,
                    "event_domain_length": 100.0},
    )
    output = tmp_path / "qiu"
    simulation = EventResolvedSimulation(config, output)
    simulation.run()
    assert simulation.full_field is not None
    assert simulation.full_field.stress.shape == (2, 2, 20, 20)
    assert np.any(simulation.full_field.eigenstrain != 0)
    assert np.any(simulation.full_field.stress != 0)
    assert np.any(simulation.driving_field != 0)
    with (tmp_path / "qiu" / "events.csv").open() as handle:
        strain_sum = sum(float(row["shear_strain_increment"]) for row in csv.DictReader(handle))
    assert np.isclose(strain_sum, simulation.accumulated_shear_strain)


def test_continuous_qiu_reference_converts_boundary_motion_to_eigenstrain(tmp_path):
    config = ModelConfig(
        regime="Q1", seed=19,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="qiu_full_field", active_modules=("qiu_reference_shear",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "easy_beta": 0.5},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "qiu-continuous")
    for segment in simulation.snapshot.boundaries.values():
        domain = simulation.domains[segment.entity_id]
        domain.previous_area_i -= 1.0
        domain.previous_area_j += 1.0
        domain.previous_time = -config.pf.time_step
    simulation.solver.step_number = 1
    simulation._update_physics()
    assert simulation.full_field is not None
    assert np.any(simulation.full_field.eigenstrain != 0)
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_mixed_event_releases_shear_and_free_volume_together(tmp_path):
    config = ModelConfig(
        regime="SC2", seed=23,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="local_memory", compatibility_model="explicit_modes",
        active_modules=("mixed_shear_climb_event",), output_cadence=1,
        max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "attempt_frequency": 100.0,
                    "barrier_core_ev": 0.0, "b_coefficient_ev": 0.0,
                    "h_coefficient_ev": 0.0, "shear_trigger": 0.1,
                    "climb_trigger_quota": 0.1},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "mixed")
    for domain in simulation.domains.values():
        domain.shear.state = 1.0
        domain.free_volume.required_total = 1.0
    simulation.solver.step_number = 1
    simulation.solver.time = config.pf.time_step
    simulation._update_physics()
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()
    with (tmp_path / "mixed" / "events.csv").open() as handle:
        events = list(csv.DictReader(handle))
    assert events
    assert all(row["barrier_type"] != "easy" for row in events)
    assert any(float(row["release_Delta_s"]) > 0 and float(row["release_Delta_q"]) > 0
               for row in events)


def test_event_simulation_checkpoint_restart_is_exact(tmp_path):
    config = ModelConfig(
        regime="G2", seed=44,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "multihit_persistent", "shear_memory", "free_volume"),
        output_cadence=1, max_steps=4, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 0,
                    "encounter_density": 2.0, "attempt_frequency": 2.0,
                    "event_domain_length": 100.0},
    )
    continuous = EventResolvedSimulation(config, tmp_path / "continuous")
    continuous.run()

    interrupted = EventResolvedSimulation(config, tmp_path / "interrupted")
    for _ in range(2):
        interrupted.solver.step()
        interrupted.snapshot = interrupted.tracker.update(interrupted.solver.labels)
        interrupted._update_physics()
    interrupted._save_checkpoint()
    interrupted.ledger.close(); interrupted.track_handle.close(); interrupted.boundary_handle.close()

    resumed = EventResolvedSimulation(config, tmp_path / "interrupted", resume=True)
    resumed.run()
    assert resumed.solver.step_number == continuous.solver.step_number
    assert (resumed.solver.eta == continuous.solver.eta).all()
    assert {k: v.state_dict() for k, v in resumed.domains.items()} == {
        k: v.state_dict() for k, v in continuous.domains.items()
    }


def test_named_temperature_and_tj_particle_regimes_are_distinct(tmp_path):
    common = dict(shape=(18, 18), interface_width=3, time_step=0.01,
                  intrinsic_mobility=1.0, adaptive_stepping=True)
    low = ModelConfig(regime="B1", seed=3, pf=PFConfig(**common, temperature=700),
                      active_modules=("arrhenius_intrinsic",), max_steps=1,
                      parameters={"initial_grains": 5, "intrinsic_barrier_ev": 0.4})
    high = ModelConfig(regime="B1", seed=3, pf=PFConfig(**common, temperature=1050),
                       active_modules=("arrhenius_intrinsic",), max_steps=1,
                       parameters={"initial_grains": 5, "intrinsic_barrier_ev": 0.4})
    low_sim = EventResolvedSimulation(low, tmp_path / "low")
    high_sim = EventResolvedSimulation(high, tmp_path / "high")
    assert high_sim.solver.config.intrinsic_mobility > low_sim.solver.config.intrinsic_mobility
    low_sim.ledger.close(); low_sim.track_handle.close()
    low_sim.boundary_handle.close()
    high_sim.ledger.close(); high_sim.track_handle.close(); high_sim.boundary_handle.close()

    tj_config = ModelConfig(
        regime="P5", seed=12, pf=PFConfig(**common, temperature=900),
        compatibility_model="off", active_modules=("gb_pinning", "tj_pinning", "random_spatial_pinning"),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 7, "equilibration_steps": 0, "particle_count": 100,
                    "particle_radius": 3.0, "encounter_density": 10.0},
    )
    simulation = EventResolvedSimulation(tj_config, tmp_path / "pins")
    simulation.run()
    assert simulation.particles is not None
    assert simulation.tj_domains
    assert (simulation.solver.mobility_scale == 0).any()


def test_equilibration_precedes_physical_time_and_hazard(tmp_path):
    config = ModelConfig(
        regime="G1", seed=5,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01, intrinsic_mobility=0.2),
        compatibility_model="geometric_surrogate", active_modules=("gb_compatibility",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 3, "encounter_density": 100.0},
    )
    output = tmp_path / "equilibrated"
    simulation = EventResolvedSimulation(config, output)
    assert simulation.solver.time == 0
    assert simulation.solver.step_number == 0
    simulation.run()
    with (output / "grain_tracks.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert any(float(row["time"]) == 0 for row in rows)
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["equilibration_steps_completed_before_time_zero"] == 3


def test_target_equilibration_compacts_phases_before_time_zero(tmp_path):
    config = ModelConfig(
        regime="B0", seed=1,
        pf=PFConfig(shape=(24, 24), interface_width=3, time_step=0.04,
                    intrinsic_mobility=2.0, adaptive_stepping=True),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 8, "equilibration_steps": 0,
                    "equilibrate_to_grains": 7, "equilibration_max_steps": 500},
    )
    output = tmp_path / "target-equilibrated"
    simulation = EventResolvedSimulation(config, output)
    assert simulation.solver.time == 0
    assert simulation.solver.eta.shape[0] == 7
    assert len(simulation.orientations) == 7
    simulation.run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["grains_after_equilibration"] == 7
    assert manifest["equilibration_steps_completed_before_time_zero"] > 0
    resumed = EventResolvedSimulation(config, output, resume=True)
    resumed.run()
    assert resumed.solver.eta.shape[0] == 7
    assert len(resumed.orientations) == 7


def test_output_cadence_does_not_change_stochastic_trajectory(tmp_path):
    common = dict(
        regime="E0", seed=88,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="explicit_modes", active_modules=("event_modes",),
        max_steps=5, termination_grains=1,
        parameters={"initial_grains": 5, "barrier_core_ev": 0.0,
                    "b_coefficient_ev": 0.0, "h_coefficient_ev": 0.0,
                    "attempt_frequency": 10.0},
    )
    frequent = EventResolvedSimulation(ModelConfig(**common, output_cadence=1), tmp_path / "frequent")
    sparse = EventResolvedSimulation(ModelConfig(**common, output_cadence=5), tmp_path / "sparse")
    frequent.run(); sparse.run()
    assert np.array_equal(frequent.solver.eta, sparse.solver.eta)
    with (tmp_path / "frequent" / "events.csv").open() as handle:
        frequent_events = [{k: v for k, v in row.items() if k != "run_id"} for row in csv.DictReader(handle)]
    with (tmp_path / "sparse" / "events.csv").open() as handle:
        sparse_events = [{k: v for k, v in row.items() if k != "run_id"} for row in csv.DictReader(handle)]
    assert frequent_events == sparse_events
