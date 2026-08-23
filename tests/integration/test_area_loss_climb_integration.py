from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.migration_closure import MigrationClosureSimulation


def _config(regime: str, modules: tuple[str, ...], max_steps: int = 2) -> ModelConfig:
    return ModelConfig(
        regime=regime,
        seed=5101,
        pf=PFConfig(
            shape=(28, 28), interface_width=4.0, time_step=0.02,
            intrinsic_mobility=0.5, adaptive_stepping=True,
        ),
        active_modules=modules,
        output_cadence=1,
        max_steps=max_steps,
        termination_grains=1,
        parameters={
            "initial_grains": 10,
            "equilibration_steps": 0,
            "event_domain_length": 6.0,
            "arclength_domains": True,
            "migration_closure": "gate_only",
            "blocked_gate_profile": "line",
            "excess_volume_per_area": 0.01,
            "point_defect_formation_volume": 0.02,
            "free_volume_stiffness": 0.05,
            "climb_trigger_quota": 0.01,
            "climb_release_quota": 1.0,
            "nucleation_prefactor": 1e12,
            "exchange_prefactor": 1e12,
            "transport_prefactor": 1e12,
            "nucleation_barrier_ev": 0.0,
            "exchange_barrier_ev": 0.0,
            "transport_barrier_ev": 0.0,
            "tj_sink_nucleation_prefactor": 1e12,
            "tj_sink_exchange_prefactor": 1e12,
            "tj_sink_transport_prefactor": 1e12,
            "tj_sink_nucleation_barrier_ev": 0.0,
            "tj_sink_exchange_barrier_ev": 0.0,
            "tj_sink_transport_barrier_ev": 0.0,
            "tj_sink_compatibility_tolerance": 1e9,
            "event_ledger_format": "csv",
        },
    )


def _close(simulation: MigrationClosureSimulation) -> None:
    simulation.ledger.close()
    simulation.track_handle.close()
    simulation.boundary_handle.close()
    simulation._activation_work_handle.close()
    if simulation._defect_history_handle is not None:
        simulation._defect_history_handle.close()
    if simulation._mechanism_state_handle is not None:
        simulation._mechanism_state_handle.close()
    if simulation.event_trace_recorder is not None:
        simulation.event_trace_recorder.close()


def test_area_loss_inventory_checkpoint_restart_is_exact(tmp_path):
    config = _config("C_GB", ("area_loss_climb", "gb_defect_sink"))
    output = tmp_path / "restart"
    simulation = MigrationClosureSimulation(config, output)
    entity = min(simulation.snapshot.boundaries)
    simulation.defect_inventory.require(1.25, entity)
    simulation.previous_total_gb_measure = simulation._total_gb_measure()
    simulation._save_checkpoint()
    expected = simulation.defect_inventory.state_dict()
    expected_measure = simulation.previous_total_gb_measure
    _close(simulation)

    resumed = MigrationClosureSimulation(config, output, resume=True)
    assert resumed.defect_inventory.state_dict() == expected
    assert resumed.previous_total_gb_measure == expected_measure
    resumed.defect_inventory.assert_conserved()
    _close(resumed)


def test_gb_sink_completion_consumes_global_quota_once(tmp_path):
    config = _config("C_GB", ("area_loss_climb", "gb_defect_sink"))
    simulation = MigrationClosureSimulation(config, tmp_path / "gb-sink")
    key = min(simulation.snapshot.boundaries)
    segment = simulation.snapshot.boundaries[key]
    domain = simulation.domains[key]
    simulation.defect_inventory.require(0.5, key)
    simulation._advance_area_loss_gb_sink(domain, segment)
    simulation._resolve_area_loss_sink_completions()
    assert np.isclose(simulation.defect_inventory.accommodated_gb, 0.5)
    assert simulation.defect_inventory.stored_total == 0.0
    simulation._advance_area_loss_gb_sink(domain, segment)
    simulation._resolve_area_loss_sink_completions()
    assert np.isclose(simulation.defect_inventory.accommodated_gb, 0.5)
    _close(simulation)


def test_tj_sink_requires_geometry_and_consumes_signed_flux_once(tmp_path):
    config = _config(
        "C_GBTJ", ("area_loss_climb", "gb_defect_sink", "tj_defect_sink")
    )
    simulation = MigrationClosureSimulation(config, tmp_path / "tj-sink")
    simulation.defect_inventory.require(0.5, "material-reservoir")
    simulation._advance_area_loss_tj_sinks()
    simulation._resolve_area_loss_sink_completions()
    assert np.isclose(simulation.defect_inventory.accommodated_tj, 0.5)
    assert simulation.defect_inventory.stored_total == 0.0
    residuals = [np.linalg.norm(tj.residual_burgers) for tj in simulation.snapshot.triple_junctions.values()]
    assert max(residuals) > 0.0
    _close(simulation)


def test_manifest_names_real_area_loss_variant(tmp_path):
    config = _config("C_GB", ("area_loss_climb", "gb_defect_sink"), max_steps=1)
    output = tmp_path / "manifest"
    MigrationClosureSimulation(config, output).run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["config"]["regime"] == "C_GB"
    assert "area_loss_climb" in manifest["config"]["active_modules"]
    assert (output / "defect_inventory.csv").exists()


def _checkpoint_without_trace(path: Path) -> tuple[dict[str, np.ndarray], dict]:
    with np.load(path / "checkpoint.npz") as archive:
        arrays = {
            key: archive[key].copy()
            for key in archive.files
            if key != "checkpoint_state_json"
        }
        state = json.loads(str(archive["checkpoint_state_json"]))
    state.get("extension_state", {}).pop("event_trace", None)
    for key in ("event_ledger_offset", "grain_tracks_offset", "boundary_tracks_offset"):
        state.pop(key, None)
    return arrays, state


def _events_without_run_id(path: Path) -> list[dict[str, str]]:
    with (path / "events.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row.pop("run_id", None)
    return rows


def test_event_trace_on_off_and_window_size_do_not_change_trajectory(tmp_path):
    modules = ("area_loss_climb", "gb_defect_sink", "tj_defect_sink")
    baseline_config = _config("C_GBTJ", modules, max_steps=8)
    baseline = tmp_path / "off"
    MigrationClosureSimulation(baseline_config, baseline, code_sha="test-sha").run()
    baseline_arrays, baseline_state = _checkpoint_without_trace(baseline)

    for name, pre_steps, post_steps in (("short", 1, 2), ("long", 5, 7)):
        traced_config = _config("C_GBTJ", modules, max_steps=8)
        traced_config.parameters.update({
            "event_trace_enabled": True,
            "event_trace_pre_steps": pre_steps,
            "event_trace_post_steps": post_steps,
            "event_trace_stride": 1,
        })
        traced = tmp_path / name
        MigrationClosureSimulation(traced_config, traced, code_sha="test-sha").run()
        traced_arrays, traced_state = _checkpoint_without_trace(traced)

        assert traced_arrays.keys() == baseline_arrays.keys()
        for key in baseline_arrays:
            assert np.array_equal(traced_arrays[key], baseline_arrays[key])
        assert traced_state == baseline_state
        assert _events_without_run_id(traced) == _events_without_run_id(baseline)
        assert len((traced / "event_trace_events.csv").read_text().splitlines()) > 1
        trace_lines = (traced / "event_traces.csv").read_text().splitlines()
        assert len(trace_lines) > 1
