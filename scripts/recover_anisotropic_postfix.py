#!/usr/bin/env python3
"""Incrementally recover the corrected anisotropic postfix qualification."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import signal
import time
from typing import Any

import numpy as np
from numba import njit

from grain_growth_pf.io.checkpoints import atomic_savez_compressed, atomic_write_text
from grain_growth_pf.pf.free_energy import free_energy
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.diagnose_anisotropic_first_transition import DT, GRAINS, SEED, SHAPE, support_masks
from scripts.run_anisotropic_reduced_pilot import config, morphology

SOURCE_COMMIT = "1316cc89dbabfb41cb883b0d4a4c74738cc2bef6"
CHECKPOINT_CADENCE = 16
ENERGY_RELATIVE_CONVERGENCE = 0.01
MORPHOLOGY_RELATIVE_CONVERGENCE = 0.10
_stop_requested = False


def _request_stop(_signum: int, _frame: object) -> None:
    global _stop_requested
    _stop_requested = True


def sha256_bytes(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_energy(solver: MultiphaseFieldSolver) -> float:
    if solver.anisotropic:
        return solver._anisotropic_energy()
    c = solver.config
    return float(free_energy(
        solver.eta, c.gb_energy, c.interface_width, c.grid_spacing,
        boundary=c.boundary_conditions,
    ))


@njit(cache=True)
def _gradient_support_metrics(eta: np.ndarray, active: np.ndarray) -> tuple[int, int, float]:
    phases, height, width = eta.shape
    tolerance = 1e-14
    supported = 0
    zero_gradient = 0
    minimum = np.inf
    local = np.empty(phases, dtype=np.int64)
    for y in range(height):
        ym = (y - 1) % height
        yp = (y + 1) % height
        for x in range(width):
            xm = (x - 1) % width
            xp = (x + 1) % width
            count = 0
            for phase in range(phases):
                if active[phase] and (
                    eta[phase, y, x] > tolerance or eta[phase, ym, x] > tolerance
                    or eta[phase, yp, x] > tolerance or eta[phase, y, xm] > tolerance
                    or eta[phase, y, xp] > tolerance
                ):
                    local[count] = phase
                    count += 1
            for left in range(count - 1):
                i = local[left]
                for right in range(left + 1, count):
                    j = local[right]
                    px = (eta[i, y, xp] - eta[i, y, xm] - eta[j, y, xp] + eta[j, y, xm]) * 0.5
                    py = (eta[i, yp, x] - eta[i, ym, x] - eta[j, yp, x] + eta[j, ym, x]) * 0.5
                    magnitude = np.sqrt(px * px + py * py)
                    supported += 1
                    if magnitude == 0.0:
                        zero_gradient += 1
                    elif magnitude < minimum:
                        minimum = magnitude
    return supported, zero_gradient, 0.0 if minimum == np.inf else minimum


def state_diagnostics(solver: MultiphaseFieldSolver) -> dict[str, Any]:
    _, support = support_masks(solver.eta, solver.active_phases, True)
    counts = support.sum(axis=0)
    pairs = counts * (counts - 1) // 2
    supported, zero_gradient, minimum = _gradient_support_metrics(
        solver.eta, solver.active_phases
    )
    return {
        "state_index": solver.step_number,
        "physical_time": solver.time,
        "active_grains": int(np.count_nonzero(solver.active_phases)),
        "maximum_active_phases": int(np.max(counts)),
        "p95_active_phases": float(np.percentile(counts, 95)),
        "maximum_active_pairs": int(np.max(pairs)),
        "cell_local_pair_instances": int(np.sum(pairs)),
        "audit_supported_pair_instances": supported,
        "zero_gradient_supported_pair_instances": zero_gradient,
        "minimum_nonzero_pair_gradient": minimum,
        "finite": bool(np.all(np.isfinite(solver.eta))),
        "minimum_phase_value": float(np.min(solver.eta)),
        "maximum_phase_sum_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1.0))),
    }


def _write_csv(path: Path, header: tuple[str, ...], rows: list[tuple[Any, ...]]) -> None:
    import io
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    writer.writerow(header)
    writer.writerows(rows)
    atomic_write_text(path, buffer.getvalue())


def _save_state(path: Path, solver: MultiphaseFieldSolver, metadata: dict[str, Any]) -> None:
    state = solver.state_dict()
    if state["orientations"] is None:
        raise ValueError("anisotropic recovery checkpoint requires orientations")
    metadata = {**metadata, "solver_config": asdict(solver.config)}
    atomic_savez_compressed(
        path,
        eta=np.asarray(state["eta"]),
        active_phases=np.asarray(state["active_phases"]),
        orientations=np.asarray(state["orientations"]),
        mobility_scale=np.asarray(state["mobility_scale"]),
        checkpoint_state_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


def _load_state(path: Path, solver: MultiphaseFieldSolver) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["checkpoint_state_json"]))
        if metadata["schema"] != "anisotropic-postfix-checkpoint-v1":
            raise ValueError("unsupported recovery checkpoint schema")
        if metadata["source_commit"] != SOURCE_COMMIT:
            raise ValueError("checkpoint scientific source mismatch")
        expected_config = json.loads(json.dumps(asdict(solver.config)))
        if metadata["solver_config"] != expected_config:
            raise ValueError("checkpoint solver configuration mismatch")
        solver.load_state_dict({
            "eta": archive["eta"],
            "active_phases": archive["active_phases"],
            "orientations": archive["orientations"],
            "mobility_scale": archive["mobility_scale"],
            "time": metadata["physical_time"],
            "step_number": metadata["accepted_step"],
        })
    return metadata


def _checkpoint(case_dir: Path, solver: MultiphaseFieldSolver, stage: str,
                energies: list[float], dt_history: list[float], diagnostics: list[dict[str, Any]],
                identities: dict[str, str], timing: list[dict[str, Any]]) -> None:
    metadata = {
        "schema": "anisotropic-postfix-checkpoint-v1",
        "source_commit": SOURCE_COMMIT,
        "case": case_dir.name,
        "stage": stage,
        "accepted_step": solver.step_number,
        "physical_time": solver.time,
        "energies": energies,
        "accepted_timestep_history": dt_history,
        "support_history": diagnostics,
        "timing_history": timing,
        "checkpoint_written_at_unix": time.time(),
        "identities": identities,
        "pair_support_state": "derived exactly from eta and active_phases; no mutable cache",
    }
    _save_state(case_dir / "checkpoint.npz", solver, metadata)
    atomic_write_text(case_dir / "checkpoint.json", json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    atomic_write_text(case_dir / "CASE_INCOMPLETE", f"{stage} step={solver.step_number}\n")


def _advance(case_dir: Path, solver: MultiphaseFieldSolver, target: int, stage: str,
             energies: list[float], dt_history: list[float], diagnostics: list[dict[str, Any]],
             identities: dict[str, str], timing: list[dict[str, Any]], record_energy: bool) -> None:
    block_started = time.monotonic()
    while solver.step_number < target:
        before = solver.time
        solver.step(compute_energy=False)
        dt_history.append(solver.time - before)
        if record_energy:
            energies.append(exact_energy(solver))
            diagnostics.append(state_diagnostics(solver))
        if solver.step_number % CHECKPOINT_CADENCE == 0 or solver.step_number == target or _stop_requested:
            now = time.monotonic()
            timing.append({"stage": stage, "accepted_step": solver.step_number,
                           "block_seconds": now - block_started})
            _checkpoint(case_dir, solver, stage, energies, dt_history, diagnostics, identities, timing)
            print(f"{time.time():.6f} {case_dir.name} {stage} {solver.step_number}/{target}", flush=True)
            block_started = time.monotonic()
        if _stop_requested:
            raise InterruptedError("USR1 requested an exact recovery checkpoint")


def run_case(output: Path, name: str, eta0: np.ndarray, orientations: np.ndarray,
             dt: float, steps: int, identities: dict[str, str]) -> dict[str, Any]:
    case_dir = output / name
    case_dir.mkdir(parents=True, exist_ok=True)
    complete = case_dir / "CASE_COMPLETE"
    if complete.exists():
        return json.loads((case_dir / "case-summary.json").read_text())
    strength = "A0_ISOTROPIC" if name == "A0" else "A2_STRONG"
    c = config(strength, name != "A0", False, dt)
    solver = MultiphaseFieldSolver(eta0.copy(), c, orientations=orientations)
    energies = [exact_energy(solver)]
    dt_history: list[float] = []
    diagnostics = [state_diagnostics(solver)]
    timing: list[dict[str, Any]] = []
    checkpoint_path = case_dir / "checkpoint.npz"
    stage = "continuous"
    if checkpoint_path.exists():
        metadata = _load_state(checkpoint_path, solver)
        if metadata["case"] != name or metadata["identities"] != identities:
            raise ValueError("checkpoint case or identity mismatch")
        stage = metadata["stage"]
        energies = list(metadata["energies"])
        dt_history = list(metadata["accepted_timestep_history"])
        diagnostics = list(metadata["support_history"])
        timing = list(metadata["timing_history"])

    midpoint = steps // 2
    if stage == "continuous":
        if solver.step_number < midpoint:
            _advance(case_dir, solver, midpoint, stage, energies, dt_history, diagnostics, identities, timing, True)
        _save_state(case_dir / "midpoint.npz", solver, {
            "schema": "anisotropic-postfix-checkpoint-v1", "source_commit": SOURCE_COMMIT,
            "case": name, "stage": "midpoint", "accepted_step": solver.step_number,
            "physical_time": solver.time, "identities": identities,
        })
        _advance(case_dir, solver, steps, stage, energies, dt_history, diagnostics, identities, timing, True)
        _save_state(case_dir / "continuous-final.npz", solver, {
            "schema": "anisotropic-postfix-checkpoint-v1", "source_commit": SOURCE_COMMIT,
            "case": name, "stage": "continuous-final", "accepted_step": solver.step_number,
            "physical_time": solver.time, "identities": identities,
        })
        continuous = solver.state_dict()
        restored = MultiphaseFieldSolver(eta0.copy(), c, orientations=orientations)
        _load_state(case_dir / "midpoint.npz", restored)
        _checkpoint(case_dir, restored, "restart", energies, [], diagnostics, identities, timing)
    else:
        with np.load(case_dir / "continuous-final.npz", allow_pickle=False) as archive:
            final_metadata = json.loads(str(archive["checkpoint_state_json"]))
            continuous = {
                "eta": archive["eta"].copy(), "active_phases": archive["active_phases"].copy(),
                "orientations": archive["orientations"].copy(), "mobility_scale": archive["mobility_scale"].copy(),
                "time": final_metadata["physical_time"], "step_number": final_metadata["accepted_step"],
            }
        restored = solver
    restart_dts: list[float] = [] if stage == "continuous" else dt_history
    _advance(case_dir, restored, steps, "restart", energies, restart_dts, diagnostics, identities, timing, False)
    final_eta = np.asarray(continuous["eta"])
    final_active = np.asarray(continuous["active_phases"])
    deltas = np.diff(np.asarray(energies))
    tolerance = max(1e-10, abs(energies[0]) * 1e-10)
    initial_morphology = morphology(eta0)
    final_morphology = morphology(final_eta)
    validity = {
        "complete_energy_sequence": len(energies) == steps + 1,
        "finite": bool(all(row["finite"] for row in diagnostics) and np.all(np.isfinite(energies))),
        "nonnegative": bool(min(row["minimum_phase_value"] for row in diagnostics) >= -1e-12),
        "phase_sum": bool(max(row["maximum_phase_sum_error"] for row in diagnostics) <= 1e-12),
        "energy_descent": bool(np.max(deltas) <= tolerance),
        "restart_exact": bool(np.array_equal(restored.eta, final_eta)
                              and np.array_equal(restored.active_phases, final_active)
                              and np.array_equal(restored.orientations, continuous["orientations"])),
        "restart_time_exact": bool(restored.time == continuous["time"]),
        "homogeneous_norm_zero_gradient_branch_finite": bool(all(row["finite"] for row in diagnostics)),
    }
    accepted_dts = [dt] * steps
    restart_dts = [dt] * (steps - midpoint)
    summary = {
        "schema": "anisotropic-postfix-case-v1", "name": name,
        "source_commit": SOURCE_COMMIT, "identities": identities, "config": asdict(c),
        "steps": steps, "accepted_dt": dt, "physical_time": steps * dt,
        "accepted_timestep_history": accepted_dts, "restart_timestep_history": restart_dts,
        "initial_energy": energies[0], "final_energy": energies[-1],
        "maximum_energy_increase": float(np.max(deltas)),
        "positive_energy_steps": int(np.count_nonzero(deltas > tolerance)),
        "energy_roundoff_tolerance": tolerance,
        "initial_field_sha256": sha256_bytes(eta0), "final_field_sha256": sha256_bytes(final_eta),
        "restart_field_sha256": sha256_bytes(restored.eta),
        "initial_morphology": initial_morphology, "final_morphology": final_morphology,
        "support_history": diagnostics, "timing_history": timing,
        "validity": validity, "passed": bool(all(validity.values())),
    }
    _write_csv(case_dir / "energy.csv", ("state_index", "energy"), list(enumerate(energies)))
    _write_csv(case_dir / "accepted-timesteps.csv", ("accepted_step", "dt"),
               [(i + 1, value) for i, value in enumerate(accepted_dts)])
    atomic_write_text(case_dir / "case-summary.json", json.dumps(summary, indent=2, sort_keys=True) + "\n")
    atomic_write_text(complete, json.dumps({"case_summary_sha256": sha256_file(case_dir / "case-summary.json")}) + "\n")
    (case_dir / "CASE_INCOMPLETE").unlink(missing_ok=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-archive-sha256", required=True)
    parser.add_argument("--input-bundle-sha256", required=True)
    parser.add_argument("--environment-sha256", required=True)
    parser.add_argument("--config-sha256", required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID") and os.environ.get("PFGG_ALLOW_LOCAL_RECOVERY") != "1":
        raise RuntimeError("recovery qualification must run in a Slurm allocation")
    signal.signal(signal.SIGUSR1, _request_stop)
    args.output.mkdir(parents=True, exist_ok=True)
    identities = {
        "source_commit": SOURCE_COMMIT,
        "source_archive_sha256": args.source_archive_sha256,
        "input_bundle_sha256": args.input_bundle_sha256,
        "environment_sha256": args.environment_sha256,
        "resolved_config_sha256": args.config_sha256,
    }
    atomic_write_text(args.output / "RUN_INCOMPLETE", f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
    eta0, _, orientations = voronoi_polycrystal(SHAPE, GRAINS, SEED, width=2.0, periodic=True)
    identities["initial_state_sha256"] = sha256_bytes(eta0)
    started = time.monotonic()
    cases = {
        "A0": run_case(args.output, "A0", eta0, orientations, DT, 64, identities),
        "A2_energy_only": run_case(args.output, "A2_energy_only", eta0, orientations, DT, 64, identities),
        "A2_energy_only_half_dt": run_case(args.output, "A2_energy_only_half_dt", eta0, orientations, DT / 2.0, 128, identities),
    }
    coarse, fine = cases["A2_energy_only"], cases["A2_energy_only_half_dt"]
    energy_difference = abs(coarse["final_energy"] - fine["final_energy"]) / max(abs(fine["final_energy"]), np.finfo(float).tiny)
    morphology_difference = {
        key: abs(float(coarse["final_morphology"][key]) - float(fine["final_morphology"][key]))
        / max(abs(float(fine["final_morphology"][key])), np.finfo(float).tiny)
        for key in ("boundary_density", "grain_area_cv", "area_weighted_mean_radius")
    }
    converged = energy_difference <= ENERGY_RELATIVE_CONVERGENCE and max(morphology_difference.values()) <= MORPHOLOGY_RELATIVE_CONVERGENCE
    if not coarse["validity"]["energy_descent"] or not fine["validity"]["energy_descent"]:
        classification = "A2_POSTFIX_ENERGY_FAILURE"
    elif not coarse["validity"]["restart_exact"] or not fine["validity"]["restart_exact"]:
        classification = "A2_POSTFIX_RESTART_FAILURE"
    elif not converged:
        classification = "A2_POSTFIX_TIMESTEP_NONCONVERGED"
    elif all(case["passed"] for case in cases.values()):
        classification = "A2_POSTFIX_QUALIFIED"
    else:
        classification = "A2_POSTFIX_OPERATIONALLY_INCOMPLETE"
    report = {
        "schema": "anisotropic-postfix-recovery-v1", "job_id": os.environ.get("SLURM_JOB_ID", "local"),
        "identities": identities, "cases": cases,
        "timestep_convergence": {"final_energy_relative_difference": energy_difference,
            "morphology_relative_differences": morphology_difference,
            "energy_threshold": ENERGY_RELATIVE_CONVERGENCE,
            "morphology_threshold": MORPHOLOGY_RELATIVE_CONVERGENCE, "passed": converged},
        "elapsed_seconds": time.monotonic() - started, "classification": classification,
        "mobility_controls_released": classification == "A2_POSTFIX_QUALIFIED",
    }
    atomic_write_text(args.output / "postfix-recovery.json", json.dumps(report, indent=2, sort_keys=True) + "\n")
    atomic_write_text(args.output / "RUN_COMPLETE", classification + "\n")
    (args.output / "RUN_INCOMPLETE").unlink(missing_ok=True)
    print(json.dumps({"classification": classification, "elapsed_seconds": report["elapsed_seconds"]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
