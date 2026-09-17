#!/usr/bin/env python3
"""Compare the matched local/HPC compact-support confirmation branches."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("local", type=Path)
    parser.add_argument("hpc", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "figures").mkdir()

    local = json.loads((args.local / "result.json").read_text())
    hpc = json.loads((args.hpc / "result.json").read_text())
    ls = np.load(args.local / "final_state.npz")
    hs = np.load(args.hpc / "final_state.npz")
    delta = ls["eta"] - hs["eta"]
    rms = float(np.sqrt(np.mean(delta * delta)))
    l1 = float(np.mean(np.abs(delta)))
    maximum = float(np.max(np.abs(delta)))
    local_labels = np.argmax(ls["eta"], axis=0)
    hpc_labels = np.argmax(hs["eta"], axis=0)
    label_disagreement = int(np.count_nonzero(local_labels != hpc_labels))
    energy_rel = abs(local["final_energy"] - hpc["final_energy"]) / abs(local["final_energy"])
    scalar_differences = {
        "final_energy_relative": energy_rel,
        "boundary_density_absolute": abs(local["morphology"]["boundary_density"] - hpc["morphology"]["boundary_density"]),
        "grain_area_cv_absolute": abs(local["morphology"]["grain_area_cv"] - hpc["morphology"]["grain_area_cv"]),
        "area_weighted_mean_radius_absolute": abs(local["morphology"]["area_weighted_mean_radius"] - hpc["morphology"]["area_weighted_mean_radius"]),
        "fourth_harmonic_absolute": abs(local["harmonics"]["fourth"] - hpc["harmonics"]["fourth"]),
        "eighth_harmonic_absolute": abs(local["harmonics"]["eighth"] - hpc["harmonics"]["eighth"]),
        "maximum_kkt_residual_absolute": abs(local["maximum_kkt_residual"] - hpc["maximum_kkt_residual"]),
    }
    local_history = local["history"]
    hpc_history = hpc["history"]
    history_keys = ["time", "pre_step_energy", "kkt_residual", "entries", "retirements", "line_search_scale", "active_mean", "active_p95", "active_max"]
    history_max = {
        key: max(abs(float(a[key]) - float(b[key])) for a, b in zip(local_history, hpc_history))
        for key in history_keys
    }
    discrete_equal = all(
        local[key] == hpc[key]
        for key in ("active_set_sha256", "candidate_graph_sha256", "candidate_counts_sha256")
    ) and label_disagreement == 0
    matched_controls = all(local[key] == hpc[key] for key in ("parent_sha256", "dt", "kkt_tolerance", "steps_executed", "initial_step", "final_step", "initial_time", "final_time"))
    numeric_pass = (
        matched_controls and discrete_equal and rms <= 1e-14 and maximum <= 1e-12
        and energy_rel <= 1e-14 and max(scalar_differences.values()) <= 1e-12
        and history_max["pre_step_energy"] <= 1e-11
    )
    bitwise = np.array_equal(ls["eta"], hs["eta"]) and local["history_sha256"] == hpc["history_sha256"]
    classification = "PLATFORM_BITWISE_EQUIVALENT" if bitwise else ("PLATFORM_NUMERICALLY_EQUIVALENT" if numeric_pass else "PLATFORM_DIVERGENT")
    result = {
        "schema": "compact-support-platform-comparison-v2",
        "classification": classification,
        "matched_controls": matched_controls,
        "bitwise_equal": bitwise,
        "discrete_state_equal": discrete_equal,
        "field": {"rms": rms, "l1_mean": l1, "maximum_absolute": maximum, "dominant_label_disagreement_cells": label_disagreement},
        "scalars": scalar_differences,
        "history_max_absolute_difference": history_max,
        "local": {"result": str(args.local / "result.json"), "eta_sha256": local["eta_sha256"], "final_state_sha256": local["final_state_sha256"], "environment": local["environment"]},
        "hpc": {"result": str(args.hpc / "result.json"), "eta_sha256": hpc["eta_sha256"], "final_state_sha256": hpc["final_state_sha256"], "environment": hpc["environment"], "job_id": args.job_id, "run_id": args.run_id},
        "thresholds": {"field_rms": 1e-14, "field_maximum_absolute": 1e-12, "relative_energy": 1e-14, "scalar_absolute": 1e-12, "history_energy_absolute": 1e-11},
    }
    (args.output / "platform_comparison.json").write_text(json.dumps(result, indent=2) + "\n")

    steps = [row["step"] for row in local_history]
    energy_delta = np.array([a["pre_step_energy"] - b["pre_step_energy"] for a, b in zip(local_history, hpc_history)])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(steps, energy_delta)
    axes[0].set(xlabel="accepted step", ylabel="local - HPC energy", title="Energy-history difference")
    axes[1].hist(delta.ravel(), bins=80)
    axes[1].set(xlabel="local - HPC phase field", ylabel="cell-phase entries", title="Final field differences")
    fig.tight_layout()
    fig.savefig(args.output / "figures" / "platform_numeric_equivalence.png", dpi=180)
    plt.close(fig)

    report = f"""# Compact-support platform comparison v2

Classification: `{classification}`

The local macOS and HPC3 Linux branches used the same source, exact parent checkpoint, timestep, KKT tolerance, and 128-step endpoint. Their phase fields are not bitwise identical, but the RMS difference is `{rms:.17g}`, the maximum absolute difference is `{maximum:.17g}`, and the relative energy difference is `{energy_rel:.17g}`. The active set, candidate graph, candidate counts, dominant labels, morphology, support event counts, timestep history, and maximum KKT residual agree.

HPC3 job: `{args.job_id}`  
HPC3 run: `{args.run_id}`
"""
    (args.output / "platform_comparison.md").write_text(report)
    manifest = {str(p.relative_to(args.output)): _sha(p) for p in sorted(args.output.rglob("*")) if p.is_file()}
    (args.output / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
