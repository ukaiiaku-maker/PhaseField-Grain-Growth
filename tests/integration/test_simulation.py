import csv
import json

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
                    "attempt_frequency": 1.0, "event_domain_length": 100.0},
    )
    output = tmp_path / "qiu"
    simulation = EventResolvedSimulation(config, output)
    simulation.run()
    assert simulation.full_field is not None
    assert simulation.full_field.stress.shape == (2, 2, 20, 20)


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
    interrupted.ledger.close(); interrupted.track_handle.close()

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
    high_sim.ledger.close(); high_sim.track_handle.close()

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
