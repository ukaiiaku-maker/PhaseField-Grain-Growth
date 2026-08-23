#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.migration_closure import MigrationClosureSimulation


def _config(regime: str, modules: tuple[str, ...], traced: bool) -> ModelConfig:
    parameters: dict[str, Any] = {
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
        "event_trace_enabled": traced,
        "event_trace_pre_steps": 25,
        "event_trace_post_steps": 50,
        "event_trace_stride": 1,
        "event_trace_format": "parquet",
    }
    return ModelConfig(
        regime=regime,
        seed=5101,
        pf=PFConfig(
            shape=(28, 28),
            interface_width=4.0,
            time_step=0.02,
            intrinsic_mobility=0.5,
            adaptive_stepping=True,
        ),
        active_modules=modules,
        output_cadence=1,
        max_steps=80,
        termination_grains=1,
        parameters=parameters,
    )


def _checkpoint(path: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
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


def _events(path: Path) -> list[dict[str, str]]:
    with (path / "events.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row.pop("run_id", None)
    return rows


def _audit(path: Path) -> dict[str, Any]:
    inventory = pd.read_csv(path / "defect_inventory.csv").iloc[-1]
    events = pd.DataFrame(_events(path))
    state = pd.read_csv(path / "mechanism_state.csv")
    completions = {}
    for event_type, column in (
        ("gb_sink_completion", "N_accommodated_GB"),
        ("tj_sink_completion", "N_accommodated_TJ"),
    ):
        selected = events[events["event_type"] == event_type]
        consumed = pd.to_numeric(
            selected.get("signed_defect_quota", pd.Series(dtype=float)), errors="coerce"
        ).abs().sum()
        completions[event_type] = {
            "events": int(len(selected)),
            "ledger_consumed": float(consumed),
            "inventory_accommodated": float(inventory[column]),
            "matches_inventory": bool(np.isclose(consumed, float(inventory[column]))),
        }
    return {
        "N_required": float(inventory["N_required"]),
        "N_accommodated_GB": float(inventory["N_accommodated_GB"]),
        "N_accommodated_TJ": float(inventory["N_accommodated_TJ"]),
        "N_stored": float(inventory["N_active_deficit"] + inventory["N_retired"]),
        "conservation_residual": float(inventory["conservation_residual"]),
        "max_abs_conservation_residual": float(inventory["max_abs_conservation_residual"]),
        "G_pending_entity_steps": int(pd.to_numeric(state["G_pending"]).sum()),
        "completion_accounting": completions,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    args = parser.parse_args()
    root = Path(args.output)
    root.mkdir(parents=True, exist_ok=False)
    sha = git_sha()

    gbtj_modules = ("area_loss_climb", "gb_defect_sink", "tj_defect_sink")
    gb_modules = ("area_loss_climb", "gb_defect_sink")
    off = root / "C_GBTJ-recorder-OFF"
    on = root / "C_GBTJ-recorder-ON"
    gb = root / "C_GB-recorder-ON"
    MigrationClosureSimulation(_config("C_GBTJ", gbtj_modules, False), off, code_sha=sha).run()
    MigrationClosureSimulation(_config("C_GBTJ", gbtj_modules, True), on, code_sha=sha).run()
    MigrationClosureSimulation(_config("C_GB", gb_modules, True), gb, code_sha=sha).run()

    off_arrays, off_state = _checkpoint(off)
    on_arrays, on_state = _checkpoint(on)
    arrays_equal = off_arrays.keys() == on_arrays.keys() and all(
        np.array_equal(off_arrays[key], on_arrays[key]) for key in off_arrays
    )
    event_trajectory_equal = _events(off) == _events(on)
    state_equal = off_state == on_state
    audits = {"C_GB": _audit(gb), "C_GBTJ": _audit(on)}
    gbtj = audits["C_GBTJ"]
    gb_audit = audits["C_GB"]
    checks = {
        "recorder_on_off_checkpoint_arrays_identical": arrays_equal,
        "recorder_on_off_physical_state_identical": state_equal,
        "recorder_on_off_event_trajectory_identical": event_trajectory_equal,
        "event_trace_nonempty": len(pd.read_parquet(on / "event_traces.parquet")) > 0,
        "exact_defect_conservation": all(
            abs(audit["max_abs_conservation_residual"]) <= 1e-10
            for audit in audits.values()
        ),
        "C_GB_zero_TJ_accommodation": gb_audit["N_accommodated_TJ"] == 0.0,
        "C_GBTJ_competing_GB_TJ_sinks": (
            gbtj["N_accommodated_GB"] > 0.0 and gbtj["N_accommodated_TJ"] > 0.0
        ),
        "pure_C_has_no_hidden_G_route": all(
            audit["G_pending_entity_steps"] == 0 for audit in audits.values()
        ),
        "no_duplicate_quota_consumption": all(
            item["matches_inventory"]
            for audit in audits.values()
            for item in audit["completion_accounting"].values()
        ),
    }
    report = {
        "status": "passed" if all(checks.values()) else "failed",
        "source_sha": sha,
        "checks": checks,
        "audits": audits,
        "runs": {"off": str(off), "on": str(on), "C_GB": str(gb)},
    }
    (root / "replay_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(root / "replay_report.json")
    if report["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
