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
