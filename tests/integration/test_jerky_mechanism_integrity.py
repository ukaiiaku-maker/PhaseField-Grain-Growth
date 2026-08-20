from pathlib import Path

import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.migration_closure import MigrationClosureSimulation


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
