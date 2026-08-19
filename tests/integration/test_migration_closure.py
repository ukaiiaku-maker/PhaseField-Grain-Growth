import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.disconnections.mode import ModeDriving
from grain_growth_pf.migration_closure import MigrationClosureSimulation


def _config(tmp_path, *, closure="gate_only", profile="line"):
    return ModelConfig(
        regime=f"closure-{closure}-{profile}",
        seed=731,
        pf=PFConfig(
            shape=(24, 24), interface_width=4.0, time_step=0.01,
            intrinsic_mobility=0.2, adaptive_stepping=True,
        ),
        output_cadence=1,
        max_steps=1,
        termination_grains=1,
        parameters={
            "initial_grains": 6,
            "equilibration_steps": 0,
            "event_domain_length": 100.0,
            "migration_closure": closure,
            "blocked_gate_profile": profile,
            "blocked_gate_halfwidth": 2.0,
            "pinned_mobility_fraction": 0.0,
        },
    )


def _close(simulation):
    simulation.ledger.close()
    simulation.track_handle.close()
    simulation.boundary_handle.close()


def test_gate_only_event_does_not_add_normal_pf_displacement(tmp_path):
    simulation = MigrationClosureSimulation(
        _config(tmp_path, closure="gate_only"), tmp_path / "gate-only"
    )
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    mode = next(mode for mode in simulation.modes if mode.step_height != 0)
    before = domain.normal_displacement_ledger
    simulation._record_event(domain, segment, mode, 1.0, ModeDriving(), "test", 0.0)
    assert domain.normal_displacement_ledger == before
    _close(simulation)


def test_legacy_hybrid_event_retains_normal_pf_displacement(tmp_path):
    simulation = MigrationClosureSimulation(
        _config(tmp_path, closure="hybrid"), tmp_path / "hybrid"
    )
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    mode = next(mode for mode in simulation.modes if mode.step_height != 0)
    simulation._record_event(domain, segment, mode, 1.0, ModeDriving(), "test", 0.0)
    expected = mode.step_height * float(simulation.config.parameters.get("packet_size", 1.0))
    assert np.isclose(domain.normal_displacement_ledger, expected)
    _close(simulation)


def test_diffuse_gate_suppresses_full_interface_band(tmp_path):
    simulation = MigrationClosureSimulation(
        _config(tmp_path, closure="gate_only", profile="diffuse"),
        tmp_path / "diffuse",
    )
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    domain.blocked = True
    simulation.solver.set_mobility_scale(1.0)
    simulation._apply_diffuse_blocked_gate()
    mobility = simulation.solver.mobility_scale
    centerline_pixels = len(np.unique(segment.points.astype(int), axis=0))
    assert np.min(mobility) == 0.0
    assert np.count_nonzero(mobility < 0.999999) > centerline_pixels
    assert np.any((mobility > 0.0) & (mobility < 1.0))
    _close(simulation)
