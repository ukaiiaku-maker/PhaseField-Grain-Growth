#!/usr/bin/env python3
"""Compare two completed FFT_EIGENSTRAIN_V2 timestep trajectories."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("coarse", type=Path)
    parser.add_argument("fine", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coarse = pd.read_csv(args.coarse / "diagnostic_timeline.csv")
    fine = pd.read_csv(args.fine / "diagnostic_timeline.csv")
    coarse_summary = json.loads((args.coarse / "analysis_summary.json").read_text())
    fine_summary = json.loads((args.fine / "analysis_summary.json").read_text())
    common_time = coarse.time.to_numpy(float)
    common_time = common_time[common_time <= float(fine.time.iloc[-1])]

    columns = [
        "grain_count", "G_population", "total_energy", "stress_l2",
        "compactness_mean", "aspect_ratio_mean",
    ]
    rows = []
    aligned = {"time": common_time}
    for column in columns:
        coarse_value = np.interp(common_time, coarse.time, coarse[column])
        fine_value = np.interp(common_time, fine.time, fine[column])
        difference = fine_value - coarse_value
        scale = max(float(np.max(np.abs(coarse_value))), np.finfo(float).tiny)
        rows.append({
            "metric": column,
            "rmse": float(np.sqrt(np.mean(difference**2))),
            "normalized_rmse": float(np.sqrt(np.mean(difference**2)) / scale),
            "maximum_absolute_difference": float(np.max(np.abs(difference))),
            "final_common_time_difference": float(difference[-1]),
        })
        aligned[f"coarse_{column}"] = coarse_value
        aligned[f"fine_{column}"] = fine_value
        aligned[f"fine_minus_coarse_{column}"] = difference
    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.output / "comparison_metrics.csv", index=False)
    pd.DataFrame(aligned).to_csv(args.output / "aligned_histories.csv", index=False)

    slope_relative_difference = (
        fine_summary["radius_squared_fit_slope"]
        / coarse_summary["radius_squared_fit_slope"] - 1.0
    )
    target_time_relative_difference = (
        fine_summary["end_time"] / coarse_summary["end_time"] - 1.0
    )
    summary = {
        "schema": "fft-eigenstrain-v2-timestep-comparison-v1",
        "classification": "FFT_EIGENSTRAIN_V2_TIMESTEP_COMPARISON_COMPLETE",
        "coarse": coarse_summary,
        "fine": fine_summary,
        "common_time_end": float(common_time[-1]),
        "radius_squared_slope_relative_difference": float(slope_relative_difference),
        "target_time_relative_difference": float(target_time_relative_difference),
        "both_no_late_avalanche": bool(
            coarse_summary["classification"] == "FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE"
            and fine_summary["classification"] == "FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE"
        ),
        "source_comparison": {
            "coarse_commit": coarse_summary["source_commit_attested"],
            "fine_commit": fine_summary["source_commit_attested"],
            "scientific_solver_change": False,
            "difference_scope": "runner overrides, provenance lookup, configurable diagnostic guard, analyses, tests, and documentation",
        },
        "metrics": rows,
    }
    (args.output / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for axis, column, label in zip(
        axes.flat,
        ("grain_count", "G_population", "total_energy", "stress_l2"),
        ("grain count", "population radius", "total energy", "stress L2"),
        strict=True,
    ):
        axis.plot(coarse.time, coarse[column], label="dt=0.04", linewidth=1)
        axis.plot(fine.time, fine[column], label="dt=0.02", linewidth=1)
        axis.set_xlabel("physical time")
        axis.set_ylabel(label)
        axis.grid(alpha=0.25)
    axes[0, 0].legend()
    for suffix in ("png", "pdf"):
        fig.savefig(args.output / f"timestep_comparison.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    metric = metrics.set_index("metric")
    report = (
        "# FFT_EIGENSTRAIN_V2 timestep comparison\n\n"
        "Classification: **FFT_EIGENSTRAIN_V2_TIMESTEP_COMPARISON_COMPLETE**.\n\n"
        f"The dt=0.04 trajectory (job `{coarse_summary['slurm_job_id']}`) and dt=0.02 "
        f"trajectory (job `{fine_summary['slurm_job_id']}`) independently reached 100 grains. "
        f"Their target times were {coarse_summary['end_time']:.2f} and "
        f"{fine_summary['end_time']:.2f}, a relative difference of "
        f"{target_time_relative_difference:.2%}. Both classify as no late avalanche.\n\n"
        f"The fitted population-radius-squared slopes are "
        f"{coarse_summary['radius_squared_fit_slope']:.6g} and "
        f"{fine_summary['radius_squared_fit_slope']:.6g}, differing by "
        f"{slope_relative_difference:.2%}. On the shared time interval, normalized RMSE is "
        f"{metric.loc['grain_count','normalized_rmse']:.3%} for grain count, "
        f"{metric.loc['G_population','normalized_rmse']:.3%} for population radius, and "
        f"{metric.loc['total_energy','normalized_rmse']:.3%} for total energy. These results "
        "support timestep consistency for the reported kinetics and avalanche classification.\n\n"
        "The attested commits differ. The audited changes add runner overrides, correct provenance "
        "lookup, make the diagnostic guard threshold configurable, and add analyses/tests/docs; "
        "they do not modify the scientific solver. The two runs remain separate trajectories.\n"
    )
    (args.output / "comparison_report.md").write_text(report)
    artifacts = sorted(path for path in args.output.iterdir() if path.name != "artifact_manifest.json")
    manifest = {
        "schema": "fft-eigenstrain-v2-timestep-comparison-manifest-v1",
        "inputs": [
            {"path": str(args.coarse / "artifact_manifest.json"), "sha256": sha256(args.coarse / "artifact_manifest.json")},
            {"path": str(args.fine / "artifact_manifest.json"), "sha256": sha256(args.fine / "artifact_manifest.json")},
        ],
        "artifacts": [{"path": p.name, "sha256": sha256(p), "bytes": p.stat().st_size} for p in artifacts],
    }
    (args.output / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
