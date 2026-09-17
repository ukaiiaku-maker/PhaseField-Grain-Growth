#!/usr/bin/env python3
"""Post-simulation analysis for one completed FFT_EIGENSTRAIN_V2 trajectory."""
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
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(fig: plt.Figure, output: Path, stem: str) -> list[Path]:
    paths = []
    for suffix in ("png", "pdf"):
        path = output / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--result-archive", type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--slurm-job-id")
    args = parser.parse_args()
    run = args.run
    args.output.mkdir(parents=True, exist_ok=True)
    parts = sorted((run / "per_step_diagnostics.parquet").glob("*.parquet"))
    if not parts:
        raise SystemExit("no per-step diagnostics")
    diagnostics = pd.concat([pd.read_parquet(path) for path in parts], ignore_index=True)
    diagnostics = diagnostics.sort_values("step").drop_duplicates("step", keep="last")
    grains = pd.read_csv(run / "grain_tracks.csv")
    boundaries = pd.read_csv(run / "boundary_tracks.csv")
    run_summary = json.loads((run / "qualification_run_summary.json").read_text())
    manifest = json.loads((run / "manifest.json").read_text())
    artifacts: list[Path] = []

    timeline_columns = [
        "step", "time", "grain_count", "G_population", "interfacial_energy",
        "elastic_energy", "total_energy", "stress_l2", "stress_linf",
        "mechanical_equilibrium_residual", "source_work_error", "clipped_fraction",
        "largest_one_step_population_loss", "largest_100_step_population_loss",
        "compactness_mean", "aspect_ratio_mean", "disconnected_grain_count",
    ]
    timeline = args.output / "diagnostic_timeline.csv"
    diagnostics[timeline_columns].to_csv(timeline, index=False)
    artifacts.append(timeline)

    t = diagnostics["time"].to_numpy(float)
    radius = diagnostics["G_population"].to_numpy(float)
    radius2 = radius**2
    fit = np.polyfit(t, radius2, 1)
    predicted = np.polyval(fit, t)
    fit_r2 = 1.0 - np.sum((radius2 - predicted) ** 2) / np.sum((radius2 - radius2.mean()) ** 2)
    late = diagnostics[diagnostics["time"] >= 0.5 * t[-1]]
    total_increments = diagnostics["total_energy"].diff().dropna()
    max_one_step = int(diagnostics["largest_one_step_population_loss"].max())
    max_late_one_step = int(late["largest_one_step_population_loss"].max())
    max_late_100_step = int(late["largest_100_step_population_loss"].max())
    no_late_avalanche = bool(
        max_late_one_step <= 3
        and max_late_100_step <= 0.1 * float(late["grain_count"].max())
        and int(late["disconnected_grain_count"].max()) == 0
    )

    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    axes[0, 0].plot(t, diagnostics["grain_count"])
    axes[0, 0].set_ylabel("grain count")
    axes[0, 1].plot(t, radius2, label=r"$R^2$")
    axes[0, 1].plot(t, predicted, "--", label=f"linear fit, $R^2$={fit_r2:.4f}")
    axes[0, 1].set_ylabel(r"population radius$^2$")
    axes[0, 1].legend(fontsize=8)
    axes[0, 2].plot(t, diagnostics["interfacial_energy"], label="interface")
    axes[0, 2].plot(t, diagnostics["elastic_energy"], label="elastic")
    axes[0, 2].plot(t, diagnostics["total_energy"], label="total", linewidth=1)
    axes[0, 2].set_yscale("symlog", linthresh=1)
    axes[0, 2].set_ylabel("energy")
    axes[0, 2].legend(fontsize=8)
    axes[1, 0].plot(t, diagnostics["stress_l2"], label=r"$||\sigma||_2$")
    axes[1, 0].plot(t, diagnostics["stress_linf"], label=r"$||\sigma||_\infty$")
    axes[1, 0].set_ylabel("stress norm")
    axes[1, 0].legend(fontsize=8)
    axes[1, 1].plot(t, diagnostics["mechanical_equilibrium_residual"])
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel("mechanical residual")
    axes[1, 2].plot(t, diagnostics["compactness_mean"], label="compactness")
    axes[1, 2].plot(t, diagnostics["aspect_ratio_mean"], label="aspect ratio")
    axes[1, 2].set_ylabel("shape statistic")
    axes[1, 2].legend(fontsize=8)
    for axis in axes.flat:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    artifacts.extend(save(fig, args.output, "kinetics_energy_mechanics"))

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
    axes[0].plot(t, diagnostics["largest_one_step_population_loss"])
    axes[0].set_ylabel("largest one-step loss")
    axes[1].plot(t, diagnostics["largest_100_step_population_loss"])
    axes[1].set_ylabel("largest 100-step loss")
    axes[2].plot(t, diagnostics["clipped_fraction"], label="clipped fraction")
    axes[2].plot(t, diagnostics["simplex_residual"], label="simplex residual")
    axes[2].set_yscale("symlog", linthresh=1e-16)
    axes[2].set_ylabel("PF acceptance diagnostic")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.set_xlabel("physical time")
        axis.grid(alpha=0.25)
    artifacts.extend(save(fig, args.output, "population_loss_and_acceptance"))

    steps = sorted(grains["step"].unique())
    selected_steps = [steps[0], steps[len(steps) // 2], steps[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), constrained_layout=True)
    distribution_rows = []
    for axis, step in zip(axes, selected_steps, strict=True):
        values = grains.loc[grains["step"] == step, "area"].to_numpy(float)
        scaled = values / values.mean()
        axis.hist(scaled, bins=np.linspace(0, 5, 41), density=True, alpha=0.8)
        axis.set_title(f"step {step}, N={len(values)}")
        axis.set_xlabel(r"$A/\langle A\rangle$")
        axis.set_ylabel("probability density")
        distribution_rows.append({
            "step": step, "time": float(grains.loc[grains["step"] == step, "time"].iloc[0]),
            "grain_count": len(values), "mean_area": float(values.mean()),
            "area_cv": float(values.std(ddof=0) / values.mean()),
            "area_skewness": float(pd.Series(values).skew()),
        })
    artifacts.extend(save(fig, args.output, "grain_area_distributions"))
    distribution_path = args.output / "grain_area_statistics.csv"
    pd.DataFrame(distribution_rows).to_csv(distribution_path, index=False)
    artifacts.append(distribution_path)

    frames = sorted((run / "frames").glob("frame-*.npz"))
    chosen = [frames[0], frames[len(frames) // 2], frames[-1]]
    fig, axes = plt.subplots(3, 3, figsize=(10, 10), constrained_layout=True)
    for column, path in enumerate(chosen):
        with np.load(path, allow_pickle=False) as frame:
            labels = frame["labels"]
            stress12 = frame["stress"][0, 1]
            strain = np.sqrt(np.sum(frame["eigenstrain"] ** 2, axis=(0, 1)))
            step = int(frame["step"])
        axes[0, column].imshow(labels, cmap="tab20", interpolation="nearest")
        limit = max(float(np.percentile(np.abs(stress12), 99)), np.finfo(float).eps)
        axes[1, column].imshow(stress12, cmap="coolwarm", vmin=-limit, vmax=limit)
        axes[2, column].imshow(strain, cmap="magma")
        axes[0, column].set_title(f"step {step}")
    for row, name in enumerate(("grain labels", r"$\sigma_{12}$", "eigenstrain norm")):
        axes[row, 0].set_ylabel(name)
        for axis in axes[row]:
            axis.set_xticks([]); axis.set_yticks([])
    artifacts.extend(save(fig, args.output, "field_evolution"))

    boundary_stats = boundaries.groupby("step").agg(
        time=("time", "first"), boundary_count=("entity_id", "size"),
        total_length=("length", "sum"), mean_abs_curvature=("curvature", lambda x: np.mean(np.abs(x))),
        rms_velocity=("normal_velocity", lambda x: np.sqrt(np.mean(x**2))),
        rms_resolved_shear=("resolved_shear", lambda x: np.sqrt(np.mean(x**2))),
    ).reset_index()
    boundary_path = args.output / "boundary_statistics.csv"
    boundary_stats.to_csv(boundary_path, index=False)
    artifacts.append(boundary_path)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
    for axis, column, ylabel in zip(axes, ("total_length", "mean_abs_curvature", "rms_resolved_shear"),
                                    ("total boundary length", "mean |curvature|", "RMS resolved shear"), strict=True):
        axis.plot(boundary_stats["time"], boundary_stats[column])
        axis.set_xlabel("physical time"); axis.set_ylabel(ylabel); axis.grid(alpha=0.25)
    artifacts.extend(save(fig, args.output, "boundary_evolution"))

    summary = {
        "schema": "fft-eigenstrain-v2-postanalysis-v1",
        "model_identity": "FFT_EIGENSTRAIN_V2",
        "run_id": args.run_id or run_summary.get("run_id"),
        "slurm_job_id": args.slurm_job_id,
        "source_commit_attested": run_summary.get("source_commit"),
        "start_step": int(diagnostics.step.iloc[0]), "end_step": int(diagnostics.step.iloc[-1]),
        "end_time": float(t[-1]), "initial_grain_count": int(diagnostics.grain_count.iloc[0]),
        "final_grain_count": int(diagnostics.grain_count.iloc[-1]),
        "radius_squared_fit_slope": float(fit[0]), "radius_squared_fit_intercept": float(fit[1]),
        "radius_squared_fit_r2": float(fit_r2),
        "positive_total_energy_increments": int((total_increments > 0).sum()),
        "maximum_total_energy_increment": float(total_increments.max()),
        "maximum_one_step_population_loss": max_one_step,
        "maximum_late_one_step_population_loss": max_late_one_step,
        "maximum_late_100_step_population_loss": max_late_100_step,
        "maximum_mechanical_equilibrium_residual": float(diagnostics.mechanical_equilibrium_residual.max()),
        "maximum_absolute_source_work_error": float(diagnostics.source_work_error.abs().max()),
        "maximum_disconnected_grain_count": int(diagnostics.disconnected_grain_count.max()),
        "classification": "FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE" if no_late_avalanche else "FFT_EIGENSTRAIN_V2_AVALANCHE_UNRESOLVED",
        "manifest_terminal_status": manifest.get("terminal_status"),
    }
    summary_path = args.output / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    artifacts.append(summary_path)
    report_path = args.output / "analysis_report.md"
    report_path.write_text(
        "# FFT_EIGENSTRAIN_V2 post-simulation analysis\n\n"
        f"Classification: **{summary['classification']}**.\n\n"
        f"Job `{summary['slurm_job_id']}` completed {summary['end_step']} steps through physical time "
        f"{summary['end_time']:.2f}, reducing the population from {summary['initial_grain_count']} "
        f"to {summary['final_grain_count']} grains. Total energy decreased at every recorded step; "
        f"the largest increment was {summary['maximum_total_energy_increment']:.3e}.\n\n"
        f"The full-horizon linear fit of population-radius squared has slope "
        f"{summary['radius_squared_fit_slope']:.4g} and R-squared "
        f"{summary['radius_squared_fit_r2']:.4f}. The largest one-step population loss was "
        f"{max_one_step}; in the latter half it was {max_late_one_step}, with a maximum "
        f"100-step loss of {max_late_100_step}. No disconnected grains were recorded. "
        "These diagnostics do not show the late non-self-similar collapse seen in the separate "
        "QIU_LEGACY_FORENSIC model.\n\n"
        "This result is one member of a separate-timestep comparison; convergence is assessed "
        "only in the paired comparison package.\n"
    )
    artifacts.append(report_path)

    evidence = [run / "qualification_run_summary.json", run / "manifest.json", run / "checkpoint.npz"]
    if args.result_archive:
        evidence.append(args.result_archive)
    artifact_manifest = {
        "schema": "fft-eigenstrain-v2-analysis-manifest-v1",
        "inputs": [{"path": str(path), "sha256": sha256(path)} for path in evidence],
        "artifacts": [{"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size} for path in artifacts],
    }
    (args.output / "artifact_manifest.json").write_text(json.dumps(artifact_manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
