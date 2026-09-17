#!/usr/bin/env python3
"""Create tables and publication figures from compact-support qualification."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_figure(fig: plt.Figure, output: Path, stem: str) -> list[Path]:
    files = []
    for suffix in ("png", "pdf"):
        path = output / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        files.append(path)
    plt.close(fig)
    return files


def label(case: dict[str, object]) -> str:
    return f"dt={case['dt']:.3g}, KKT={case['kkt_tolerance']:.0e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qualification", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = json.loads(args.qualification.read_text())
    if report.get("schema") != "compact-support-qualification-v1":
        raise SystemExit("unexpected qualification schema")
    args.output.mkdir(parents=True, exist_ok=True)
    cases = report["cases"]
    artifacts: list[Path] = []

    table_path = args.output / "case_summary.csv"
    columns = (
        "name", "dt", "kkt_tolerance", "steps", "physical_time", "stop_reason",
        "maximum_checkpoint_energy_increment", "maximum_kkt_residual",
        "maximum_exact_active_phases", "maximum_p95_exact_active_phases",
        "late_early_cost_ratio", "late_A2_A0_cost_ratio", "restart_exact",
        "restart_time_exact", "all_release_targets",
    )
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for case in cases:
            row = {key: case.get(key) for key in columns}
            row["all_release_targets"] = all(case["release_targets"].values())
            writer.writerow(row)
    artifacts.append(table_path)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for case in cases:
        records = case["records"]
        time = np.asarray([row["physical_time"] for row in records])
        axes[0, 0].plot(time, [row["active_mean"] for row in records], label=label(case))
        axes[0, 1].plot(time, [row["active_p95"] for row in records])
        axes[1, 0].plot(time, [row["candidate_mean"] for row in records])
        axes[1, 1].plot(time, [row["active_pair_instances"] for row in records])
    axes[0, 0].set_ylabel("mean active phases/cell")
    axes[0, 1].set_ylabel("p95 active phases/cell")
    axes[1, 0].set_ylabel("mean candidate phases/cell")
    axes[1, 1].set_ylabel("active pair instances")
    for axis in axes.flat:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7, ncol=2)
    artifacts.extend(write_figure(fig, args.output, "support_evolution"))

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    for case in cases:
        records = case["records"]
        time = np.asarray([row["physical_time"] for row in records])
        energy = np.asarray([row["exact_energy"] for row in records])
        axes[0].plot(time, energy - energy[0], label=label(case))
        axes[1].plot(time[1:], [row["runtime_per_accepted_step"] for row in records[1:]])
        axes[2].plot(time, [row["boundary_density"] for row in records])
    axes[0].set_ylabel(r"$F(t)-F(0)$")
    axes[1].set_ylabel("seconds/accepted step")
    axes[2].set_ylabel("boundary density")
    for axis in axes:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    artifacts.extend(write_figure(fig, args.output, "energy_cost_morphology"))

    principal = [case for case in cases if case["kkt_tolerance"] == 1e-10]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), constrained_layout=True)
    metrics = (
        ("maximum_p95_exact_active_phases", "max p95 active", 8.0),
        ("late_early_cost_ratio", "late/early cost", 10.0),
        ("late_A2_A0_cost_ratio", "late A2/A0 cost", 20.0),
    )
    for axis, (key, ylabel, target) in zip(axes, metrics, strict=True):
        values = [case[key] for case in principal]
        axis.bar([f"{case['dt']:.3g}" for case in principal], values, color="#35618d")
        axis.axhline(target, color="#b33a3a", linestyle="--", linewidth=1)
        axis.set_xlabel("dt")
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.25)
    artifacts.extend(write_figure(fig, args.output, "release_targets"))

    report_path = args.output / "qualification_report.md"
    lines = [
        "# Compact-support qualification report", "",
        f"Classification: **{report['classification']}**", "",
        f"Slurm job: `{report['job_id']}`  ",
        f"Source: `{report['source_commit']}`  ",
        f"Input evidence: `{sha256(args.qualification)}`", "",
        "| case | stop | max dF | max KKT | max/p95 support | cost late/early | cost A2/A0 | targets |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for case in cases:
        targets = all(case["release_targets"].values())
        lines.append(
            f"| {case['name']} | {case['stop_reason']} | "
            f"{case['maximum_checkpoint_energy_increment']:.3e} | "
            f"{case['maximum_kkt_residual']:.3e} | "
            f"{case['maximum_exact_active_phases']}/{case['maximum_p95_exact_active_phases']:.2f} | "
            f"{case['late_early_cost_ratio']:.2f} | {case['late_A2_A0_cost_ratio']:.2f} | "
            f"{'PASS' if targets else 'FAIL'} |"
        )
    lines.extend([
        "", "## Gate interpretation", "",
        f"Tolerance convergence: **{'PASS' if report['tolerance_exact'] else 'FAIL'}**.  ",
        f"Exact restart: **{'PASS' if report['restart_pass'] else 'FAIL'}**.  ",
        f"All release targets: **{'PASS' if report['release_targets_pass'] else 'FAIL'}**.", "",
        "The Qiu anisotropic production gate opens only for "
        "`COMPACT_SUPPORT_OPERATOR_QUALIFIED`. Scientific failure outputs remain "
        "complete evidence and must not be reclassified as infrastructure failures.",
    ])
    report_path.write_text("\n".join(lines) + "\n")
    artifacts.append(report_path)

    manifest_path = args.output / "artifact_manifest.json"
    manifest = {
        "schema": "compact-support-analysis-manifest-v1",
        "input": {"path": str(args.qualification), "sha256": sha256(args.qualification)},
        "classification": report["classification"],
        "artifacts": [
            {"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in artifacts
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
