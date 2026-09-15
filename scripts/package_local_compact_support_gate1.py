#!/usr/bin/env python3
"""Package a completed local compact-support Gate-1 qualification."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column) for column in columns} for row in rows)


def case_rows(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for case in cases:
        for record in case["records"]:
            area = np.asarray(record["grain_area_distribution"], dtype=float)
            radii = np.sqrt(area / np.pi)
            harmonic = record["crystal_frame_interface_harmonics"]
            rows.append({
                "case": case["name"], "dt": case["dt"],
                "kkt_tolerance": case["kkt_tolerance"], **record,
                "grain_area_mean": float(np.mean(area)),
                "grain_area_std": float(np.std(area)),
                "grain_radius_mean": float(np.mean(radii)),
                "grain_radius_std": float(np.std(radii)),
                "crystal_fourth_harmonic": harmonic["fourth"],
                "crystal_eighth_harmonic": harmonic["eighth"],
                "interface_sample_count": harmonic["sample_count"],
            })
    return rows


def classify(report: dict[str, object], provenance_complete: bool) -> str:
    cases = report["cases"]
    if not provenance_complete or not report["restart_pass"]:
        return "COMPACT_SUPPORT_OPERATOR_OPERATIONALLY_INCOMPLETE"
    if any(not case["release_targets"]["energy_descent"] for case in cases):
        return "COMPACT_SUPPORT_OPERATOR_ENERGY_FAILURE"
    if not report["tolerance_exact"] or any(
        not case["release_targets"]["kkt"] for case in cases
    ):
        return "COMPACT_SUPPORT_OPERATOR_KKT_FAILURE"
    topology_keys = ("finite", "nonnegative", "phase_sum")
    if any(
        not all(case["release_targets"][key] for key in topology_keys)
        for case in cases
    ):
        return "COMPACT_SUPPORT_OPERATOR_TOPOLOGY_FAILURE"
    scaling_keys = (
        "p95_active", "max_active", "late_early_cost", "support_plateau",
        "late_A2_A0_cost",
    )
    if any(
        not all(case["release_targets"][key] for key in scaling_keys)
        for case in cases
    ):
        return "COMPACT_SUPPORT_OPERATOR_SCALING_FAILURE"
    return "COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qualification", type=Path)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = json.loads(args.qualification.read_text())
    run_record = json.loads((args.run_dir / "status/run-record.json").read_text())
    environment = json.loads((args.run_dir / "status/environment.json").read_text())
    lineage = json.loads((args.run_dir / "status/checkpoint-lineage.json").read_text())
    complete_marker = args.run_dir / "output/qualification/RUN_COMPLETE"
    provenance_complete = (
        complete_marker.exists()
        and run_record["start_checkpoint_sha256"] == lineage["predecessor_checkpoint_sha256"]
        and run_record["scientific_source_commit"] == report["source_commit"]
        and run_record["reconciliation"]["sacct_identity_found"] is False
    )
    classification = classify(report, provenance_complete)
    args.output.mkdir(parents=True, exist_ok=True)
    figures = args.output / "figures"
    checkpoints = args.output / "checkpoints"
    figures.mkdir(exist_ok=True)
    checkpoints.mkdir(exist_ok=True)

    rows = case_rows(report["cases"])
    common = ["case", "dt", "kkt_tolerance", "accepted_step", "physical_time"]
    write_csv(args.output / "energy.csv", common + ["exact_energy", "energy_increment"], rows)
    write_csv(args.output / "support.csv", common + [
        "active_mean", "active_p95", "active_max", "candidate_mean",
        "candidate_p95", "candidate_max", "active_pair_instances",
        "candidate_pair_instances", "entries", "retirements",
    ], rows)
    write_csv(args.output / "kkt.csv", common + [
        "kkt_residual", "minimum_phase_value", "phase_sum_error",
        "line_search_scale",
    ], rows)
    write_csv(args.output / "morphology.csv", common + [
        "grain_count", "boundary_density", "grain_area_mean", "grain_area_std",
        "grain_radius_mean", "grain_radius_std", "crystal_fourth_harmonic",
        "crystal_eighth_harmonic", "interface_sample_count",
    ], rows)
    write_csv(args.output / "timestep.csv", common + ["line_search_scale"], rows)
    write_csv(args.output / "runtime.csv", common + [
        "elapsed_seconds", "runtime_per_accepted_step",
    ], rows)

    decision = {
        "schema": "local-compact-support-gate1-decision-v1",
        "classification": classification,
        "upstream_classification": report["classification"],
        "provenance_complete": provenance_complete,
        "restart_pass": report["restart_pass"],
        "tolerance_exact": report["tolerance_exact"],
        "release_targets_pass": report["release_targets_pass"],
        "short_hpc3_confirmation_required": classification == "COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED",
        "native_qiu_gate_open": False,
    }
    (args.output / "qualification_decision.json").write_text(
        json.dumps(decision, indent=2) + "\n"
    )
    final_manifest = dict(run_record)
    final_manifest.update({
        "state": "COMPLETE", "local_classification": classification,
        "qualification_sha256": sha256(args.qualification),
        "terminal_marker_sha256": sha256(complete_marker),
    })
    (args.output / "run_manifest.json").write_text(
        json.dumps(final_manifest, indent=2) + "\n"
    )
    shutil.copyfile(args.run_dir / "status/environment.json", args.output / "environment.json")
    shutil.copyfile(
        args.run_dir / "status/checkpoint-lineage.json",
        args.output / "checkpoint_lineage.json",
    )
    for item in lineage["checkpoints"]:
        path = checkpoints / f"{item['case']}-step{item['accepted_step']}.json"
        path.write_text(json.dumps(item, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for case in report["cases"]:
        records = case["records"]
        time = [row["physical_time"] for row in records]
        label = f"dt={case['dt']:.3g}, KKT={case['kkt_tolerance']:.0e}"
        axes[0, 0].plot(time, [row["exact_energy"] for row in records], label=label)
        axes[0, 1].plot(time, [row["active_mean"] for row in records])
        axes[1, 0].plot(time, [row["kkt_residual"] for row in records])
        axes[1, 1].plot(time, [row["runtime_per_accepted_step"] for row in records])
    axes[0, 0].set_ylabel("total energy")
    axes[0, 1].set_ylabel("mean exact support")
    axes[1, 0].set_ylabel("KKT residual")
    axes[1, 1].set_ylabel("seconds / accepted step")
    for axis in axes.flat:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=6)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"gate1_histories.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), constrained_layout=True)
    for case in report["cases"]:
        records = case["records"]
        time = [row["physical_time"] for row in records]
        axes[0].plot(time, [row["boundary_density"] for row in records])
        axes[1].plot(time, [row["grain_count"] for row in records])
        axes[2].plot(time, [row["crystal_frame_interface_harmonics"]["fourth"] for row in records])
    axes[0].set_ylabel("boundary density")
    axes[1].set_ylabel("grain count")
    axes[2].set_ylabel("crystal-frame fourth harmonic")
    for axis in axes:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    for suffix in ("png", "pdf"):
        fig.savefig(figures / f"morphology_histories.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = []
    for path in sorted(args.output.rglob("*")):
        if path.is_file() and path.name != "checksum_manifest.json":
            manifest.append({
                "path": str(path.relative_to(args.output)),
                "sha256": sha256(path), "bytes": path.stat().st_size,
            })
    (args.output / "checksum_manifest.json").write_text(
        json.dumps({"schema": "local-gate1-checksums-v1", "files": manifest}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
