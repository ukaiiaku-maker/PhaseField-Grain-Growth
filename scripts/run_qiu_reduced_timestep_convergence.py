#!/usr/bin/env python3
"""Long reduced-domain FFT-v2 comparison for two external-increment targets."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow.dataset as ds

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.pf.initial_conditions import prepare_initial_condition
from run_qiu_fft_v2_qualification import FFTEigenstrainFrameSimulation


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--max-steps", type=int, default=5000)
    parser.add_argument("--termination-grains", type=int, default=5)
    parser.add_argument("--shape", type=int, default=48)
    parser.add_argument("--grains", type=int, default=18)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source_sha = git_sha(Path(__file__).resolve().parents[1])
    pf = PFConfig(
        shape=(args.shape, args.shape), grid_spacing=1.0, interface_width=4.0,
        time_step=0.04, gb_energy=1.0, intrinsic_mobility=4.0,
        boundary_conditions="periodic", temperature=900.0,
        adaptive_stepping=True, grain_extinction_threshold=0.5,
    )
    initial_parameters = {
        "initial_grains": args.grains, "equilibration_steps": 0,
        "compact_after_equilibration": True,
    }
    initial = args.output / "initial_state.npz"
    prepare_initial_condition(pf, 5101, initial_parameters, initial, source_sha)
    histories: dict[str, object] = {}
    rows: list[dict[str, object]] = []
    for target in (0.02, 0.01):
        label = f"target-{target:g}"
        parameters = {
            **initial_parameters, "initial_state_file": str(initial.resolve()),
            "full_field_source_mapping": "local_sweep",
            "coupled_integration_enabled": True,
            "elastic_shear_modulus": 1.0, "elastic_poisson_ratio": 0.3,
            "elastic_constitutive_state": "plane_strain",
            "external_delta_eta_target": target,
            "coupled_minimum_dt": 1e-8, "coupled_max_rejections": 20,
            "coupled_energy_relative_tolerance": 1e-10,
            "coupled_energy_absolute_tolerance": 1e-10,
            "rigid_shift_search": 2, "max_local_phases": 8,
            "checkpoint_cadence": 250, "energy_diagnostic_cadence": 1,
            "video_frame_cadence": 100,
            "qiu_diagnostics_enabled": True, "qiu_diagnostic_flush_rows": 32,
            "qiu_diagnostic_field_start_step": args.max_steps + 1,
            "qiu_guard_terminate": True, "qiu_guard_clip_warmup_steps": 20,
            "qiu_guard_population_loss_fraction": 1.0,
            "qiu_guard_clip_fraction": 0.02, "qiu_guard_clip_factor": 1.5,
            "qiu_guard_clip_additive": 0.05, "qiu_guard_extinction_count": 10,
        }
        config = ModelConfig(
            regime="FFT_EIGENSTRAIN_V2", seed=5101, pf=pf,
            mechanics_backend="fft_eigenstrain_v2", compatibility_model="off",
            active_modules=("fft_eigenstrain_shear",), output_cadence=100,
            max_steps=args.max_steps, termination_grains=args.termination_grains,
            parameters=parameters,
        )
        run = args.output / label
        started = time.perf_counter()
        simulation = FFTEigenstrainFrameSimulation(config, run, code_sha=source_sha)
        simulation.run()
        elapsed = time.perf_counter() - started
        frame = ds.dataset(
            run / "per_step_diagnostics.parquet", format="parquet"
        ).to_table().to_pandas().sort_values("time")
        histories[label] = frame
        final = frame.iloc[-1]
        rows.append({
            "label": label, "target": target, "steps": int(final["step"]),
            "physical_time": float(final["time"]),
            "grain_count": int(final["grain_count"]),
            "G_population": float(final["G_population"]),
            "interfacial_energy": float(final["interfacial_energy"]),
            "elastic_energy": float(final["elastic_energy"]),
            "total_energy": float(final["total_energy"]),
            "stress_p95": float(final["stress_p95"]),
            "compactness_mean": float(final["compactness_mean"]),
            "max_equilibrium_residual": float(frame["mechanical_equilibrium_residual"].max()),
            "max_abs_work_error": float(frame["source_work_error"].abs().max()),
            "max_rejections": int(frame["coupled_rejection_count"].max()),
            "minimum_used_dt": float(frame["used_dt"].min()),
            "minimum_external_dt_limit": float(frame["external_dt_limit"].min()),
            "external_limit_active_steps": int(np.count_nonzero(
                np.asarray(frame["external_dt_limit"], dtype=float) < 0.04
            )),
            "maximum_100_step_loss": int(frame["largest_100_step_population_loss"].max()),
            "guard_triggered": (run / "diagnostic_capture.json").exists(),
            "wall_seconds": elapsed, "run_path": str(run.resolve()),
        })

    final_time = min(float(row["physical_time"]) for row in rows)
    columns = (
        "G_population", "interfacial_energy", "elastic_energy", "total_energy",
        "stress_p95", "compactness_mean",
    )
    matched: dict[str, dict[str, float]] = {}
    for label, frame in histories.items():
        matched[label] = {
            column: float(np.interp(final_time, frame["time"], frame[column]))
            for column in columns
        }
    coarse = matched["target-0.02"]; fine = matched["target-0.01"]
    relative = {
        column: abs(fine[column] - coarse[column]) /
        max(abs(fine[column]), np.finfo(float).eps) for column in columns
    }
    gates = {
        "target_refinement_exercised": any(
            int(row["external_limit_active_steps"]) > 0 for row in rows
        ),
        "G_within_1pct": relative["G_population"] <= 0.01,
        "interfacial_energy_within_2pct": relative["interfacial_energy"] <= 0.02,
        "elastic_energy_within_2pct": relative["elastic_energy"] <= 0.02,
        "compactness_within_5pct": relative["compactness_mean"] <= 0.05,
        "stress_p95_within_5pct": relative["stress_p95"] <= 0.05,
    }
    plot_dir = args.output / "plots"; plot_dir.mkdir()
    for column in columns:
        figure, axis = plt.subplots(figsize=(6.5, 4.2))
        for label, frame in histories.items():
            axis.plot(frame["time"], frame[column], label=label)
        axis.set_xlabel("physical time"); axis.set_ylabel(column); axis.legend()
        figure.tight_layout(); figure.savefig(plot_dir / f"{column}.png", dpi=180)
        plt.close(figure)
    summary = {
        "schema_version": 1, "source_commit": source_sha,
        "initial_state": str(initial.resolve()), "initial_state_sha256": sha256(initial),
        "shape": [args.shape, args.shape], "initial_grains": args.grains,
        "runs": rows, "matched_time": final_time, "matched_observables": matched,
        "fine_relative_to_coarse": relative, "production_threshold_gates": gates,
        "scope": "long reduced-domain target convergence; not production qualification",
    }
    atomic_write_text(args.output / "timestep_convergence_summary.json", json.dumps(summary, indent=2) + "\n")
    checksums = {
        str(path.relative_to(args.output)): sha256(path)
        for path in sorted(args.output.rglob("*")) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
