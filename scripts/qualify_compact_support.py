#!/usr/bin/env python3
"""HPC3 qualification of the compact-support anisotropic operator."""
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

import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.io.checkpoints import atomic_savez_compressed, atomic_write_text
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.diagnose_anisotropic_first_transition import DT, GRAINS, SEED, SHAPE
from scripts.run_anisotropic_reduced_pilot import (
    ENERGY_NORMALIZATION, MOBILITY_NORMALIZATION, morphology,
)

SOURCE_COMMIT = os.environ.get("PFGG_SOURCE_COMMIT", "UNSPECIFIED")
PHYSICAL_TIME_CEILING = 0.25
CHECK_INTERVAL = 64
CHECKPOINT_INTERVAL = 512
KKT_TOLERANCES = (1e-8, 1e-10, 1e-12)
STOP_REQUESTED = False


def request_stop(_signum: int, _frame: object) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True


def hash_array(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    return hashlib.sha256(value.view(np.uint8)).hexdigest()


def config(dt: float, tolerance: float) -> PFConfig:
    return PFConfig(
        shape=SHAPE, grid_spacing=1.0, interface_width=4.0, time_step=dt,
        gb_energy=1.0, intrinsic_mobility=4.0, boundary_conditions="periodic",
        adaptive_stepping=False, grain_extinction_threshold=0.5,
        anisotropy_strength="A2_STRONG", anisotropic_energy=True,
        anisotropic_mobility=True,
        anisotropy_energy_normalization=ENERGY_NORMALIZATION,
        anisotropy_mobility_normalization=MOBILITY_NORMALIZATION,
        anisotropic_support_mode="compact_active_set",
        anisotropic_kkt_tolerance=tolerance,
    )


def grain_areas(eta: np.ndarray) -> list[float]:
    return [float(value) for value in eta.sum(axis=(1, 2)) if value > 1e-10]


def orientation_harmonics(eta: np.ndarray, orientations: np.ndarray) -> dict[str, float]:
    labels = np.argmax(eta, axis=0)
    total = 0
    fourth = 0j
    eighth = 0j
    for axis in (0, 1):
        shifted = np.roll(labels, -1, axis=axis)
        mask = shifted != labels
        ys, xs = np.nonzero(mask)
        for y, x in zip(ys, xs, strict=True):
            i, j = labels[y, x], shifted[y, x]
            normal = 0.0 if axis == 1 else np.pi / 2.0
            angle = normal - 0.5 * (orientations[i] + orientations[j])
            fourth += np.exp(4j * angle)
            eighth += np.exp(8j * angle)
            total += 1
    return {
        "fourth": float(abs(fourth) / max(total, 1)),
        "eighth": float(abs(eighth) / max(total, 1)),
        "sample_count": total,
    }


def state_record(
    solver: MultiphaseFieldSolver, elapsed: float, previous_energy: float | None,
    block_seconds: float = 0.0, block_steps: int = 0,
) -> dict[str, object]:
    energy = solver._anisotropic_energy()
    support = solver._last_compact_support
    if support is None:
        exact = (solver.eta > 0.0).sum(axis=0)
        support_data = {
            "active_mean": float(np.mean(exact)), "active_p95": float(np.percentile(exact, 95)),
            "active_max": int(np.max(exact)), "candidate_mean": None,
            "candidate_p95": None, "candidate_max": None,
            "active_pair_instances": int(np.sum(exact * (exact - 1) // 2)),
            "candidate_pair_instances": None, "kkt_residual": 0.0,
            "entries": 0, "retirements": 0, "line_search_scale": 1.0,
        }
    else:
        support_data = asdict(support)
    metric = morphology(solver.eta)
    return {
        "accepted_step": solver.step_number, "physical_time": solver.time,
        "exact_energy": energy,
        "energy_increment": None if previous_energy is None else energy - previous_energy,
        **support_data,
        "phase_sum_error": float(np.max(np.abs(solver.eta.sum(axis=0) - 1.0))),
        "minimum_phase_value": float(np.min(solver.eta)),
        "grain_count": metric["active_grains"],
        "boundary_density": metric["boundary_density"],
        "grain_area_distribution": grain_areas(solver.eta),
        "crystal_frame_interface_harmonics": orientation_harmonics(solver.eta, solver.orientations),
        "elapsed_seconds": elapsed,
        "runtime_per_accepted_step": block_seconds / max(block_steps, 1),
    }


def stop_reason(initial: dict[str, object], current: dict[str, object]) -> str | None:
    if current["grain_count"] <= 0.9 * initial["grain_count"]:
        return "grain_count_decreased_10_percent"
    if abs(current["boundary_density"] / initial["boundary_density"] - 1.0) >= 0.05:
        return "boundary_density_changed_5_percent"
    initial_radius = np.sqrt(np.asarray(initial["grain_area_distribution"]) / np.pi)
    current_areas = np.asarray(current["grain_area_distribution"])
    current_radius = np.sqrt(current_areas / np.pi)
    r0 = float(np.mean(initial_radius)); r1 = float(np.mean(current_radius))
    if abs(r1 / r0 - 1.0) >= 0.05:
        return "mean_grain_radius_changed_5_percent"
    if current["physical_time"] >= PHYSICAL_TIME_CEILING:
        return "preregistered_A0_based_physical_time_ceiling"
    return None


def save_checkpoint(case_dir: Path, solver: MultiphaseFieldSolver,
                    records: list[dict[str, object]], stage: str) -> None:
    state = solver.state_dict()
    metadata = {
        "schema": "compact-support-checkpoint-v1", "source_commit": SOURCE_COMMIT,
        "case": case_dir.name, "stage": stage, "time": solver.time,
        "step": solver.step_number, "config": asdict(solver.config), "records": records,
        "candidate_graph": "deterministic function of eta and boundary conditions",
    }
    atomic_savez_compressed(
        case_dir / "checkpoint.npz", eta=state["eta"],
        active_phases=state["active_phases"], orientations=state["orientations"],
        mobility_scale=state["mobility_scale"],
        checkpoint_state_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    atomic_write_text(case_dir / "checkpoint.json", json.dumps(metadata, indent=2) + "\n")


def load_checkpoint(path: Path, solver: MultiphaseFieldSolver) -> tuple[list[dict[str, object]], str]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["checkpoint_state_json"]))
        if metadata["schema"] != "compact-support-checkpoint-v1":
            raise ValueError("checkpoint schema mismatch")
        if metadata["config"] != json.loads(json.dumps(asdict(solver.config))):
            raise ValueError("checkpoint configuration mismatch")
        solver.load_state_dict({
            "eta": archive["eta"], "active_phases": archive["active_phases"],
            "orientations": archive["orientations"], "mobility_scale": archive["mobility_scale"],
            "time": metadata["time"], "step_number": metadata["step"],
        })
    return list(metadata["records"]), str(metadata["stage"])


def run_case(output: Path, eta0: np.ndarray, orientations: np.ndarray,
             dt: float, tolerance: float) -> dict[str, object]:
    name = f"dt-{dt:.16g}_kkt-{tolerance:.0e}"
    case_dir = output / name
    case_dir.mkdir(parents=True, exist_ok=True)
    summary_path = case_dir / "case-summary.json"
    if summary_path.exists():
        return json.loads(summary_path.read_text())
    solver = MultiphaseFieldSolver(eta0.copy(), config(dt, tolerance), orientations=orientations)
    records: list[dict[str, object]] = []
    started = time.monotonic()
    stage = "continuous"
    checkpoint = case_dir / "checkpoint.npz"
    if checkpoint.exists():
        records, stage = load_checkpoint(checkpoint, solver)
    elif (case_dir / "midpoint.npz").exists():
        records, stage = load_checkpoint(case_dir / "midpoint.npz", solver)
    if not records:
        records.append(state_record(solver, 0.0, None))
    initial = records[0]
    midpoint_path = case_dir / "midpoint.npz"
    reason = None
    block_started = time.monotonic()
    block_step = solver.step_number
    while reason is None:
        solver.step(compute_energy=False)
        if solver.step_number % CHECK_INTERVAL == 0:
            now = time.monotonic()
            records.append(state_record(
                solver, now - started, records[-1]["exact_energy"],
                now - block_started, solver.step_number - block_step,
            ))
            block_started = now
            block_step = solver.step_number
            reason = stop_reason(initial, records[-1])
        if not midpoint_path.exists() and solver.time >= PHYSICAL_TIME_CEILING / 2.0:
            save_checkpoint(case_dir, solver, records, "midpoint")
            checkpoint.replace(midpoint_path)
        if solver.step_number % CHECKPOINT_INTERVAL == 0 or reason or STOP_REQUESTED:
            save_checkpoint(case_dir, solver, records, "continuous")
            print(name, solver.step_number, solver.time, flush=True)
        if STOP_REQUESTED:
            raise InterruptedError("USR1 requested checkpoint")
    continuous = solver.state_dict()

    restart_exact = None
    restart_time_exact = None
    if tolerance == 1e-10:
        replay = MultiphaseFieldSolver(eta0.copy(), config(dt, tolerance), orientations=orientations)
        _, _ = load_checkpoint(midpoint_path, replay)
        while replay.step_number < solver.step_number:
            replay.step(compute_energy=False)
        restart_exact = bool(np.array_equal(replay.eta, continuous["eta"]))
        restart_time_exact = bool(replay.time == continuous["time"])

    active_max = max(int(row["active_max"]) for row in records)
    active_p95 = max(float(row["active_p95"]) for row in records)
    kkt_max = max(float(row["kkt_residual"]) for row in records)
    increments = [row["energy_increment"] for row in records[1:]]
    early = [row["runtime_per_accepted_step"] for row in records[1:3]]
    late = [row["runtime_per_accepted_step"] for row in records[-2:]]
    cost_ratio = float(np.mean(late) / max(np.mean(early), np.finfo(float).tiny))
    late_support = np.asarray([row["active_mean"] for row in records[-5:]], dtype=float)
    support_plateau = bool(
        late_support.size >= 3
        and np.ptp(late_support) <= 0.10 * max(float(np.mean(late_support)), np.finfo(float).tiny)
        and not np.all(np.diff(late_support) > 0.0)
    )
    finite = all(
        np.isfinite(row["exact_energy"])
        and np.isfinite(row["minimum_phase_value"])
        and np.isfinite(row["phase_sum_error"])
        for row in records
    )
    summary = {
        "schema": "compact-support-case-v1", "name": name, "source_commit": SOURCE_COMMIT,
        "dt": dt, "kkt_tolerance": tolerance, "stop_reason": reason,
        "steps": solver.step_number, "physical_time": solver.time,
        "initial_field_sha256": hash_array(eta0), "final_field_sha256": hash_array(solver.eta),
        "restart_exact": restart_exact, "restart_time_exact": restart_time_exact,
        "maximum_checkpoint_energy_increment": max(increments),
        "maximum_kkt_residual": kkt_max, "maximum_exact_active_phases": active_max,
        "maximum_p95_exact_active_phases": active_p95, "late_early_cost_ratio": cost_ratio,
        "records": records,
        "release_targets": {
            "energy_descent": max(increments) <= 128 * np.finfo(float).eps * abs(records[0]["exact_energy"]),
            "kkt": kkt_max <= tolerance, "p95_active": active_p95 <= 8,
            "max_active": active_max <= 16, "late_early_cost": cost_ratio <= 10,
            "support_plateau": support_plateau, "finite": finite,
            "nonnegative": min(row["minimum_phase_value"] for row in records) >= 0.0,
            "phase_sum": max(row["phase_sum_error"] for row in records) <= 1e-12,
        },
    }
    atomic_write_text(summary_path, json.dumps(summary, indent=2) + "\n")
    with (case_dir / "history.csv").open("w", newline="") as handle:
        scalar_fields = [k for k, v in records[0].items() if not isinstance(v, (list, dict))]
        writer = csv.DictWriter(handle, fieldnames=scalar_fields, lineterminator="\n")
        writer.writeheader()
        for row in records:
            writer.writerow({k: row[k] for k in scalar_fields})
    return summary


def run_a0_reference(output: Path, eta0: np.ndarray, orientations: np.ndarray) -> dict[str, object]:
    case_dir = output / "A0-runtime-reference"
    case_dir.mkdir(parents=True, exist_ok=True)
    summary_path = case_dir / "case-summary.json"
    if summary_path.exists():
        return json.loads(summary_path.read_text())
    cfg = PFConfig(
        shape=SHAPE, grid_spacing=1.0, interface_width=4.0, time_step=DT,
        gb_energy=1.0, intrinsic_mobility=4.0, boundary_conditions="periodic",
        adaptive_stepping=False, grain_extinction_threshold=0.5,
        anisotropy_strength="A0_ISOTROPIC", anisotropic_energy=False,
        anisotropic_mobility=False,
    )
    solver = MultiphaseFieldSolver(eta0.copy(), cfg, orientations=orientations)
    blocks = []
    while solver.time < PHYSICAL_TIME_CEILING:
        started = time.monotonic()
        start_step = solver.step_number
        for _ in range(CHECK_INTERVAL):
            solver.step(compute_energy=False)
        seconds = time.monotonic() - started
        blocks.append({"end_step": solver.step_number, "end_time": solver.time,
                       "seconds": seconds, "seconds_per_step": seconds / (solver.step_number - start_step)})
        if STOP_REQUESTED:
            raise InterruptedError("USR1 requested checkpoint")
    summary = {
        "schema": "compact-support-A0-runtime-v1", "steps": solver.step_number,
        "physical_time": solver.time, "blocks": blocks,
        "early_seconds_per_step": float(np.mean([row["seconds_per_step"] for row in blocks[:2]])),
        "late_seconds_per_step": float(np.mean([row["seconds_per_step"] for row in blocks[-2:]])),
        "final_morphology": morphology(solver.eta),
    }
    atomic_write_text(summary_path, json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("compact-support qualification must run in Slurm")
    signal.signal(signal.SIGUSR1, request_stop)
    args.output.mkdir(parents=True, exist_ok=True)
    atomic_write_text(args.output / "RUN_INCOMPLETE", os.environ["SLURM_JOB_ID"] + "\n")
    eta0, _, orientations = voronoi_polycrystal(SHAPE, GRAINS, SEED, width=2.0, periodic=True)
    a0 = run_a0_reference(args.output, eta0, orientations)
    cases = []
    for tolerance in KKT_TOLERANCES:
        for dt in (DT, DT / 2.0):
            cases.append(run_case(args.output, eta0, orientations, dt, tolerance))
    by_dt = {}
    for case in cases:
        by_dt.setdefault(case["dt"], []).append(case)
    tolerance_exact = all(
        len({case["final_field_sha256"] for case in group}) == 1 for group in by_dt.values()
    )
    principal = [case for case in cases if case["kkt_tolerance"] == 1e-10]
    for case in cases:
        late = float(np.mean([row["runtime_per_accepted_step"] for row in case["records"][-2:]]))
        case["late_A2_A0_cost_ratio"] = late / max(a0["late_seconds_per_step"], np.finfo(float).tiny)
        case["release_targets"]["late_A2_A0_cost"] = case["late_A2_A0_cost_ratio"] <= 20
    restart_pass = all(case["restart_exact"] and case["restart_time_exact"] for case in principal)
    targets = all(all(case["release_targets"].values()) for case in cases)
    energy_pass = all(case["release_targets"]["energy_descent"] for case in cases)
    topology_pass = all(
        case["release_targets"]["finite"]
        and case["release_targets"]["nonnegative"]
        and case["release_targets"]["phase_sum"]
        for case in cases
    )
    classification = (
        "COMPACT_SUPPORT_OPERATOR_QUALIFIED" if targets and tolerance_exact and restart_pass
        else "COMPACT_SUPPORT_OPERATOR_ENERGY_FAILURE" if not energy_pass
        else "COMPACT_SUPPORT_OPERATOR_TOLERANCE_NONCONVERGED" if not tolerance_exact
        else "COMPACT_SUPPORT_OPERATOR_TOPOLOGY_FAILURE" if not topology_pass
        else "COMPACT_SUPPORT_OPERATOR_SCALING_FAILURE"
    )
    report = {
        "schema": "compact-support-qualification-v1", "job_id": os.environ["SLURM_JOB_ID"],
        "source_commit": SOURCE_COMMIT, "physical_time_ceiling": PHYSICAL_TIME_CEILING,
        "ceiling_basis": "A0 64-step boundary-density rate predicts 5% change near t=0.28; ceiling conservatively set to 0.25",
        "A0_runtime_reference": a0, "cases": cases,
        "tolerance_exact": tolerance_exact, "restart_pass": restart_pass,
        "release_targets_pass": targets, "classification": classification,
    }
    atomic_write_text(args.output / "compact_support_qualification.json", json.dumps(report, indent=2) + "\n")
    atomic_write_text(args.output / "RUN_COMPLETE", classification + "\n")
    (args.output / "RUN_INCOMPLETE").unlink(missing_ok=True)
    print(classification, flush=True)
    # A scientific gate failure is a valid completed qualification result.  A
    # nonzero process status would make the runner classify and retrieve the
    # evidence as an incomplete infrastructure failure, obscuring the actual
    # scientific classification written above.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
