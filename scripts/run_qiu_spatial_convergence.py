#!/usr/bin/env python3
"""Run a physical-domain-preserving FFT-v2 grid-refinement comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow.dataset as ds

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.simulation import EventResolvedSimulation


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--physical-size", type=float, default=48.0)
    parser.add_argument("--physical-time", type=float, default=4.0)
    parser.add_argument("--grains", type=int, default=18)
    parser.add_argument("--seed", type=int, default=5101)
    parser.add_argument(
        "--dxs", default="1,0.5,0.25",
        help="Comma-separated grid spacings, ordered coarse to fine.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source_sha = git_sha(Path(__file__).resolve().parents[1])
    histories: dict[str, object] = {}
    rows: list[dict[str, object]] = []

    spacings = tuple(float(value) for value in args.dxs.split(","))
    if len(spacings) < 2 or any(dx <= 0.0 for dx in spacings):
        raise SystemExit("--dxs requires at least two positive grid spacings")
    for dx in spacings:
        label = f"dx-{dx:g}"
        pixels = int(round(args.physical_size / dx))
        dt = 0.04 * dx * dx
        steps = int(round(args.physical_time / dt))
        pf = PFConfig(
            shape=(pixels, pixels), grid_spacing=dx, interface_width=4.0,
            time_step=dt, gb_energy=1.0, intrinsic_mobility=4.0,
            boundary_conditions="periodic", temperature=900.0,
            adaptive_stepping=True, grain_extinction_threshold=0.5,
        )
        # The RNG draws are identical in normalized coordinates at both grids;
        # the compact smoothing width is converted from physical units to cells.
        eta, seeds, orientations = voronoi_polycrystal(
            pf.shape, args.grains, args.seed,
            width=pf.interface_width / (2.0 * pf.grid_spacing), periodic=True,
        )
        initial = args.output / f"initial-{label}.npz"
        np.savez_compressed(
            initial, eta=eta, orientations=orientations, seed_positions=seeds,
            active_original_ids=np.arange(args.grains), equilibration_steps=np.asarray(0),
        )
        parameters = {
            "initial_grains": args.grains, "equilibration_steps": 0,
            "compact_after_equilibration": True,
            "initial_state_file": str(initial.resolve()),
            "full_field_source_mapping": "local_sweep",
            "coupled_integration_enabled": True,
            "elastic_shear_modulus": 1.0, "elastic_poisson_ratio": 0.3,
            "elastic_constitutive_state": "plane_strain",
            "external_delta_eta_target": 0.02,
            "coupled_minimum_dt": 1e-8, "coupled_max_rejections": 20,
            "coupled_energy_relative_tolerance": 1e-10,
            "coupled_energy_absolute_tolerance": 1e-10,
            "rigid_shift_search": 2, "max_local_phases": 8,
            "checkpoint_cadence": steps,
            "energy_diagnostic_cadence": 1,
            "qiu_diagnostics_enabled": True,
            "qiu_diagnostic_flush_rows": 32,
            "qiu_diagnostic_field_start_step": steps + 1,
            "qiu_guard_terminate": False,
            "qiu_guard_clip_warmup_steps": 20,
        }
        config = ModelConfig(
            regime="FFT_EIGENSTRAIN_V2", seed=args.seed, pf=pf,
            mechanics_backend="fft_eigenstrain_v2", compatibility_model="off",
            active_modules=("fft_eigenstrain_shear",), output_cadence=steps,
            max_steps=steps, termination_grains=5, parameters=parameters,
        )
        run = args.output / label
        EventResolvedSimulation(config, run, code_sha=source_sha).run()
        frame = ds.dataset(
            run / "per_step_diagnostics.parquet", format="parquet"
        ).to_table().to_pandas()
        histories[label] = frame
        final = frame.iloc[-1]
        rows.append({
            "label": label, "dx": dx, "pixels": pixels, "dt": dt,
            "steps": steps, "physical_time": float(final["time"]),
            "grain_count": int(final["grain_count"]),
            "G_population": float(final["G_population"]) * dx,
            "interfacial_energy_density": float(final["interfacial_energy"]) / args.physical_size**2,
            "elastic_energy_density": float(final["elastic_energy"]) / args.physical_size**2,
            "total_energy_density": float(final["total_energy"]) / args.physical_size**2,
            "stress_p95": float(final["stress_p95"]),
            "stress_linf": float(final["stress_linf"]),
            "compactness_mean": float(final["compactness_mean"]),
            "compactness_p95": float(final["compactness_p95"]),
            "max_equilibrium_residual": float(frame["mechanical_equilibrium_residual"].max()),
            "max_abs_work_error": float(frame["source_work_error"].abs().max()),
            "max_clipped_fraction": float(frame["clipped_fraction"].max()),
            "initial_state": str(initial.resolve()),
            "initial_state_sha256": sha256(initial), "run_path": str(run.resolve()),
        })

    comparison_keys = (
        "G_population", "interfacial_energy_density", "elastic_energy_density",
        "total_energy_density", "stress_p95", "compactness_mean",
        "compactness_p95",
    )
    adjacent = []
    for coarse, fine in zip(rows[:-1], rows[1:]):
        adjacent.append({
            "coarse": coarse["label"], "fine": fine["label"],
            "relative_differences": {
                key: abs(float(fine[key]) - float(coarse[key])) /
                max(abs(float(fine[key])), np.finfo(float).eps)
                for key in comparison_keys
            },
        })
    relative = adjacent[-1]["relative_differences"]
    gates = {
        "G_within_1pct": relative["G_population"] <= 0.01,
        "interfacial_energy_within_2pct": relative["interfacial_energy_density"] <= 0.02,
        "elastic_energy_within_2pct": relative["elastic_energy_density"] <= 0.02,
        "compactness_mean_within_5pct": relative["compactness_mean"] <= 0.05,
        "stress_p95_within_5pct": relative["stress_p95"] <= 0.05,
    }
    plots = args.output / "plots"; plots.mkdir()
    specs = (
        ("G_population", "population grain size", True),
        ("total_energy", "total energy density", True),
        ("stress_p95", "stress p95", False),
        ("compactness_mean", "mean compactness", False),
    )
    for column, ylabel, scale_length in specs:
        fig, axis = plt.subplots(figsize=(6.5, 4.2))
        for row in rows:
            frame = histories[str(row["label"])]
            values = np.asarray(frame[column], dtype=float)
            if column == "total_energy":
                values /= args.physical_size**2
            elif scale_length:
                values *= float(row["dx"])
            axis.plot(frame["time"], values, label=str(row["label"]))
        axis.set_xlabel("physical time"); axis.set_ylabel(ylabel); axis.legend()
        fig.tight_layout(); fig.savefig(plots / f"{column}.png", dpi=180); plt.close(fig)

    summary = {
        "schema_version": 1, "source_commit": source_sha,
        "study": "fixed physical domain, interface width, grain density, and normalized seeds",
        "physical_size": args.physical_size, "physical_time": args.physical_time,
        "initial_grains": args.grains, "runs": rows,
        "adjacent_grid_comparisons": adjacent,
        "finest_pair_relative_differences": relative,
        "production_threshold_gates": gates,
        "scope": "reduced-domain spatial convergence; not a full production trajectory",
    }
    atomic_write_text(args.output / "spatial_convergence_summary.json", json.dumps(summary, indent=2) + "\n")
    checksums = {
        str(path.relative_to(args.output)): sha256(path)
        for path in sorted(args.output.rglob("*")) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
