import json
from dataclasses import replace

import numpy as np
import pyarrow.dataset as ds

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.simulation import EventResolvedSimulation
from run_qiu_fft_v2_qualification import FFTEigenstrainFrameSimulation


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
        "qiu_guard_clip_warmup_steps": 0,
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


def test_qiu_guard_continue_saves_one_guard_field_then_returns_to_cadence(tmp_path):
    config = _config(diagnostics=True)
    config = replace(config, max_steps=4, parameters={
        **config.parameters, "qiu_guard_clip_fraction": -1.0,
        "qiu_guard_clip_warmup_steps": 0, "qiu_guard_terminate": False,
    })
    output = tmp_path / "guard-continue"
    EventResolvedSimulation(config, output, code_sha="same").run()

    manifest = json.loads((output / "manifest.json").read_text())
    fields = sorted(path.name for path in (output / "diagnostic_fields").glob("*.npz"))
    assert manifest["status"] == "completed"
    assert manifest["steps_completed"] == 4
    assert fields == [
        "step-0000001-guard.npz",
        "step-0000002-cadence.npz",
        "step-0000004-cadence.npz",
    ]


def _fft_v2_config(max_steps=4):
    return ModelConfig(
        regime="FFT_EIGENSTRAIN_V2", seed=97,
        pf=PFConfig(
            shape=(24, 24), interface_width=3.0, time_step=0.01,
            intrinsic_mobility=0.1, adaptive_stepping=True,
        ),
        mechanics_backend="fft_eigenstrain_v2", compatibility_model="off",
        active_modules=("fft_eigenstrain_shear",), output_cadence=2,
        max_steps=max_steps, termination_grains=1,
        parameters={
            "initial_grains": 7, "equilibration_steps": 0,
            "checkpoint_cadence": 2, "elastic_constitutive_state": "plane_strain",
        },
    )


def test_fft_v2_simulation_uses_distributed_source_and_work_conjugate_force(tmp_path):
    simulation = EventResolvedSimulation(_fft_v2_config(), tmp_path / "fft-v2")
    simulation.run()
    assert simulation.fft_coupling is not None
    assert np.count_nonzero(simulation.full_field.eigenstrain) > 4
    assert np.any(simulation.driving_field != 0.0)
    assert simulation.full_field.last_equilibrium_residual < 1e-10


def test_fft_v2_checkpoint_restart_is_exact(tmp_path):
    parameters = {
        **_fft_v2_config().parameters, "external_delta_eta_target": 1e-12,
    }
    continuous_config = replace(_fft_v2_config(4), parameters=parameters)
    first_config = replace(_fft_v2_config(2), parameters=parameters)
    continuous = tmp_path / "continuous"
    resumed = tmp_path / "resumed"
    EventResolvedSimulation(continuous_config, continuous, code_sha="same").run()
    EventResolvedSimulation(first_config, resumed, code_sha="same").run()
    EventResolvedSimulation(continuous_config, resumed, resume=True, code_sha="same").run()
    left, left_state = _checkpoint(continuous)
    right, right_state = _checkpoint(resumed)
    for name in left:
        assert np.array_equal(left[name], right[name])
    with np.load(continuous / "checkpoint.npz") as a, np.load(resumed / "checkpoint.npz") as b:
        assert np.array_equal(a["driving_field"], b["driving_field"])
    assert left_state["time"] == right_state["time"]
    assert left_state["time"] < 4 * continuous_config.pf.time_step


def test_fft_v2_diagnostics_and_output_cadence_are_trajectory_invariant(tmp_path):
    base = _fft_v2_config(5)
    diagnostic = replace(base, parameters={
        **base.parameters,
        "qiu_diagnostics_enabled": True,
        "qiu_diagnostic_field_start_step": 100,
        "qiu_guard_clip_fraction": 1.0,
        "qiu_guard_extinction_count": 1000,
    })
    sparse_output = replace(base, output_cadence=3)
    paths = [tmp_path / name for name in ("base", "diagnostic", "sparse-output")]
    for config, path in zip((base, diagnostic, sparse_output), paths):
        EventResolvedSimulation(config, path, code_sha="same").run()
    checkpoints = [_checkpoint(path)[0] for path in paths]
    for candidate in checkpoints[1:]:
        for name in checkpoints[0]:
            assert np.array_equal(checkpoints[0][name], candidate[name])
    with np.load(paths[0] / "checkpoint.npz") as left:
        for path in paths[1:]:
            with np.load(path / "checkpoint.npz") as right:
                assert np.array_equal(left["driving_field"], right["driving_field"])


def test_fft_v2_external_limit_reduces_dt_and_complete_energy_never_increases(tmp_path):
    base = _fft_v2_config(6)
    config = replace(base, parameters={
        **base.parameters,
        "external_delta_eta_target": 1e-12,
        "qiu_diagnostics_enabled": True,
        "qiu_diagnostic_field_start_step": 100,
        "qiu_guard_clip_fraction": 1.0,
        "qiu_guard_extinction_count": 1000,
    })
    output = tmp_path / "limited"
    EventResolvedSimulation(config, output, code_sha="same").run()
    scalar = ds.dataset(output / "per_step_diagnostics.parquet", format="parquet").to_table().to_pandas()
    assert np.any(scalar["used_dt"].to_numpy()[1:] < config.pf.time_step)
    assert np.all(scalar["max_raw_order_parameter_increment"] < 0.03)
    energy = json.loads((output / "energy.json").read_text())
    total = np.asarray([row["total_complete"] for row in energy])
    assert np.all(np.diff(total) <= 2e-10 * np.maximum(total[:-1], 1.0))


def test_fft_v2_converges_as_external_increment_target_is_tightened(tmp_path):
    def evolve(target, name):
        base = _fft_v2_config(10000)
        config = replace(
            base,
            pf=replace(base.pf, shape=(16, 16), intrinsic_mobility=0.2),
            parameters={
                **base.parameters, "initial_grains": 4,
                "external_delta_eta_target": target,
                "elastic_shear_modulus": 1e4,
            },
        )
        simulation = EventResolvedSimulation(config, tmp_path / name, code_sha="same")
        steps = 0
        while simulation.solver.time < 0.05 - 1e-14:
            simulation._advance_fft_v2_step(
                maximum_dt=0.05 - simulation.solver.time
            )
            simulation.snapshot = simulation.tracker.update(simulation.solver.labels)
            simulation._update_physics(fft_source_applied=True)
            steps += 1
        eta = simulation.solver.eta.copy()
        simulation.ledger.close()
        simulation.track_handle.close()
        simulation.boundary_handle.close()
        return eta, steps

    targets = (1e-5, 5e-6, 2.5e-6, 1.25e-6)
    results = [evolve(target, f"target-{index}") for index, target in enumerate(targets)]
    reference = results[-1][0]
    errors = [np.linalg.norm(result[0] - reference) for result in results[:-1]]
    assert errors[0] > errors[1] > errors[2]
    assert results[0][1] < results[1][1] < results[2][1] < results[3][1]


def test_fft_v2_movie_path_uses_base_coupling_and_writes_restart_frames(tmp_path):
    config = replace(_fft_v2_config(3), parameters={
        **_fft_v2_config().parameters, "video_frame_cadence": 2,
    })
    output = tmp_path / "movie"
    FFTEigenstrainFrameSimulation(config, output, code_sha="same").run()
    frames = sorted((output / "frames").glob("frame-*.npz"))
    assert [int(path.stem.rsplit("-", 1)[1]) for path in frames] == [0, 2, 3]
    with np.load(frames[-1]) as frame:
        assert set(("labels", "stress", "eigenstrain", "step", "time")).issubset(frame.files)
