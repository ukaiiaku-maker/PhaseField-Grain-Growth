from pathlib import Path
import csv

import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.migration_closure import MigrationClosureSimulation


def _close(simulation):
    simulation.ledger.close()
    simulation.track_handle.close()
    simulation.boundary_handle.close()
    simulation._activation_work_handle.close()


def _base(tmp_path: Path, regime: str, modules=(), compatibility="off", mechanics="none"):
    return ModelConfig(
        regime=regime,
        seed=919,
        pf=PFConfig(
            shape=(32, 32), interface_width=4.0, time_step=0.02,
            intrinsic_mobility=0.5, adaptive_stepping=True,
        ),
        mechanics_backend=mechanics,
        compatibility_model=compatibility,
        active_modules=tuple(modules),
        output_cadence=2,
        max_steps=12,
        termination_grains=1,
        parameters={
            "initial_grains": 10,
            "equilibration_steps": 0,
            "event_domain_length": 6.0,
            "encounter_density": 0.5,
            "migration_closure": "gate_only",
            "blocked_gate_profile": "line",
            "arclength_domains": True,
            "tj_correlation_radius": 0,
            "tj_correlation_length": 2.0,
            "shear_stiffness": 0.2,
            "easy_beta": 0.5,
            "attempt_frequency": 1e5,
            "barrier_core_ev": 0.15,
            "b_coefficient_ev": 0.01,
            "h_coefficient_ev": 0.01,
        },
    )


def test_pure_tj_case_does_not_advance_hidden_gb_encounter_clock(tmp_path):
    config = _base(
        tmp_path, "T-only", modules=("tj_compatibility", "multihit_persistent"),
        compatibility="geometric_surrogate",
    )
    simulation = MigrationClosureSimulation(config, tmp_path / "tj-only")
    simulation.run()
    assert simulation.domains
    assert all(np.isclose(domain.encounter.total_measure, 0.0) for domain in simulation.domains.values())


def test_local_shear_memory_accumulates_without_geometric_barrier(tmp_path):
    config = _base(
        tmp_path, "S-only", modules=("shear_memory",), mechanics="local_memory"
    )
    simulation = MigrationClosureSimulation(config, tmp_path / "shear-only")
    simulation.run()
    states = np.asarray([domain.shear.state for domain in simulation.domains.values()])
    assert states.size
    assert np.max(np.abs(states)) > 0.0


def test_activation_work_diagnostic_is_written_for_gb_release(tmp_path):
    config = _base(
        tmp_path, "G-only", modules=("gb_compatibility", "single_hit_poisson"),
        compatibility="geometric_surrogate",
    )
    config = ModelConfig.from_dict({
        **config.to_dict(),
        "max_steps": 40,
        "parameters": {
            **config.parameters,
            "encounter_density": 5.0,
            "attempt_frequency": 1e8,
            "barrier_core_ev": 0.01,
        },
    })
    output = tmp_path / "gb-only"
    simulation = MigrationClosureSimulation(config, output)
    simulation.run()
    work = output / "activation_work.csv"
    assert work.exists()
    assert work.stat().st_size > len(",".join(simulation.ACTIVATION_WORK_FIELDS))


def test_tj_release_uses_local_shear_activation_work(tmp_path):
    config = _base(
        tmp_path,
        "TS",
        modules=("tj_compatibility", "multihit_persistent", "shear_memory"),
        compatibility="geometric_surrogate",
        mechanics="local_memory",
    )
    output = tmp_path / "tj-shear-work"
    simulation = MigrationClosureSimulation(config, output)
    tj = next(iter(simulation.snapshot.triple_junctions.values()))
    tj_domain = simulation.tj_domains[tj.entity_id]
    packet = float(config.parameters.get("packet_size", 1.0))
    mode = min(
        simulation.modes,
        key=lambda candidate: np.linalg.norm(
            packet * np.asarray(candidate.burgers) + tj.residual_burgers
        ),
    )
    b_direction = np.asarray(mode.burgers, dtype=float)
    b_direction /= np.linalg.norm(b_direction)
    # tangent=(-normal_y, normal_x)=b_direction for every adjoining domain.
    normal = (float(b_direction[1]), float(-b_direction[0]))
    for boundary_id in tj.adjoining_boundaries:
        simulation.snapshot.boundaries[boundary_id].normal = normal
        simulation.domains[boundary_id].shear.state = 1.0

    driving = simulation._tj_mode_driving(tj, mode)
    assert not np.isclose(driving.resolved_shear, 0.0)

    tj_domain.blocked = True
    tj_domain.activation.hit_count = tj_domain.hits - 1
    tj_domain.activation.clock.threshold = 0.0
    simulation._update_tj_physics(np.ones(config.pf.shape))
    simulation._activation_work_handle.flush()

    with (output / "activation_work.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    tj_rows = [row for row in rows if row["event_type"] == "tj_compatibility_release"]
    assert tj_rows
    assert not np.isclose(float(tj_rows[-1]["work_shear"]), 0.0)
    assert np.isclose(
        float(tj_rows[-1]["effective_DeltaG"]),
        max(
            0.0,
            float(config.parameters.get("tj_barrier_ev", 0.6))
            - mode.activation_work_ev(driving),
        ),
    )
    _close(simulation)


def test_corrected_tj_gate_ignores_legacy_pixel_radius(tmp_path):
    config = _base(
        tmp_path,
        "T",
        modules=("tj_compatibility", "multihit_persistent"),
        compatibility="geometric_surrogate",
    )
    config = ModelConfig.from_dict({
        **config.to_dict(),
        "parameters": {
            **config.parameters,
            "tj_correlation_radius": 12,
            "tj_correlation_length": 0.0,
        },
    })
    simulation = MigrationClosureSimulation(config, tmp_path / "physical-tj-gate")
    for domain in simulation.tj_domains.values():
        domain.blocked = False
    tj = next(iter(simulation.snapshot.triple_junctions.values()))
    simulation.tj_domains[tj.entity_id].blocked = True
    mobility = np.ones(config.pf.shape)
    simulation._apply_physical_tj_gate(mobility)
    assert simulation._tj_gate_radius_pixels() == 0
    assert np.count_nonzero(mobility == 0.0) == 1
    _close(simulation)
