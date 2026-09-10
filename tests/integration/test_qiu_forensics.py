import json
from dataclasses import replace

import numpy as np
import pyarrow.dataset as ds

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.simulation import EventResolvedSimulation


def _config(*, diagnostics: bool) -> ModelConfig:
    parameters = {
        "initial_grains": 7,
        "equilibration_steps": 0,
        "event_domain_length": 100.0,
        "easy_beta": 0.35,
        "checkpoint_cadence": 2,
        "energy_diagnostic_cadence": 1,
    }
    if diagnostics:
        parameters.update({
            "qiu_diagnostics_enabled": True,
            "qiu_diagnostic_flush_rows": 2,
            "qiu_diagnostic_field_start_step": 2,
            "qiu_diagnostic_field_cadence": 2,
            "qiu_guard_clip_fraction": 1.0,
            "qiu_guard_extinction_count": 1000,
        })
    return ModelConfig(
        regime="QIU-FORENSIC-TEST", seed=20260910,
        pf=PFConfig(
            shape=(24, 24), interface_width=3.0, time_step=0.01,
            intrinsic_mobility=0.1, adaptive_stepping=True,
        ),
        mechanics_backend="qiu_full_field", compatibility_model="off",
        active_modules=("qiu_reference_shear",), output_cadence=2,
        max_steps=4, termination_grains=1, parameters=parameters,
    )


def _checkpoint(path):
    with np.load(path / "checkpoint.npz") as state:
        return {
            name: state[name].copy()
            for name in ("eta", "eigenstrain", "active_phases")
        }, json.loads(str(state["checkpoint_state_json"]))


def test_qiu_forensics_is_bitwise_trajectory_invariant(tmp_path):
    plain_path = tmp_path / "plain"
    forensic_path = tmp_path / "forensic"
    EventResolvedSimulation(_config(diagnostics=False), plain_path, code_sha="same").run()
    EventResolvedSimulation(_config(diagnostics=True), forensic_path, code_sha="same").run()

    plain, plain_state = _checkpoint(plain_path)
    forensic, forensic_state = _checkpoint(forensic_path)
    for name in plain:
        assert np.array_equal(plain[name], forensic[name])
    assert plain_state["time"] == forensic_state["time"]
    assert plain_state["step_number"] == forensic_state["step_number"]
    assert np.count_nonzero(plain["active_phases"]) == np.count_nonzero(
        forensic["active_phases"]
    )

    scalar = ds.dataset(
        forensic_path / "per_step_diagnostics.parquet", format="parquet"
    ).to_table().to_pandas()
    assert scalar["step"].tolist() == [1, 2, 3, 4]
    assert {
        "max_raw_order_parameter_increment", "clipped_fraction",
        "mechanical_equilibrium_residual", "source_work_error",
        "compactness_p95", "extinction_cluster_size",
    }.issubset(scalar.columns)
    assert (forensic_path / "diagnostic_fields" / "step-0000002-cadence.npz").exists()
    assert (forensic_path / "diagnostic_fields" / "step-0000004-cadence.npz").exists()


def test_qiu_forensics_checkpoint_resume_has_no_duplicate_steps(tmp_path):
    output = tmp_path / "resume"
    first = _config(diagnostics=True)
    first = replace(first, max_steps=2)
    EventResolvedSimulation(first, output, code_sha="same").run()
    EventResolvedSimulation(_config(diagnostics=True), output, resume=True, code_sha="same").run()
    table = ds.dataset(output / "per_step_diagnostics.parquet", format="parquet").to_table()
    assert table.column("step").to_pylist() == [1, 2, 3, 4]


def test_qiu_forensics_can_start_from_an_uninstrumented_legacy_checkpoint(tmp_path):
    output = tmp_path / "late-enable"
    first = replace(_config(diagnostics=False), max_steps=2)
    EventResolvedSimulation(first, output, code_sha="same").run()
    EventResolvedSimulation(_config(diagnostics=True), output, resume=True, code_sha="same").run()
    table = ds.dataset(output / "per_step_diagnostics.parquet", format="parquet").to_table()
    assert table.column("step").to_pylist() == [3, 4]


def test_qiu_guard_saves_full_state_and_terminates_as_diagnostic_capture(tmp_path):
    config = _config(diagnostics=True)
    config = replace(config, max_steps=4, parameters={
        **config.parameters, "qiu_guard_clip_fraction": -1.0,
    })
    output = tmp_path / "guard"
    simulation = EventResolvedSimulation(config, output, code_sha="same")
    simulation.run()
    manifest = json.loads((output / "manifest.json").read_text())
    capture = json.loads((output / "diagnostic_capture.json").read_text())
    assert manifest["status"] == "diagnostic_capture"
    assert manifest["steps_completed"] == 1
    assert capture["scientific_status"] == "diagnostic_capture_not_completed"
    assert (output / "checkpoint.npz").exists()
    assert (output / "diagnostic_fields" / "step-0000001-guard.npz").exists()
