#!/usr/bin/env python3
"""Execute the reduced one-change-at-a-time QIU root-cause matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.io.checkpoints import atomic_write_text
from grain_growth_pf.io.provenance import git_sha
from grain_growth_pf.pf.initial_conditions import prepare_initial_condition
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
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--shape", type=int, default=48)
    parser.add_argument("--grains", type=int, default=18)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    sha = git_sha()
    pf = PFConfig(
        shape=(args.shape, args.shape), grid_spacing=1.0,
        interface_width=4.0, time_step=0.04, gb_energy=1.0,
        intrinsic_mobility=4.0, boundary_conditions="periodic",
        temperature=900.0, adaptive_stepping=True,
        grain_extinction_threshold=0.5,
    )
    initial_parameters = {
        "initial_grains": args.grains, "equilibration_steps": 0,
        "compact_after_equilibration": True,
    }
    initial = args.output / "initial_state.npz"
    prepare_initial_condition(pf, 5101, initial_parameters, initial, sha)
    common = {
        **initial_parameters,
        "initial_state_file": str(initial.resolve()),
        "event_domain_length": 12.0,
        "easy_beta": 0.35,
        "checkpoint_cadence": 200,
        "energy_diagnostic_cadence": 1,
        "qiu_diagnostics_enabled": True,
        "qiu_diagnostic_flush_rows": 32,
        "qiu_diagnostic_field_start_step": args.steps + 1,
        "qiu_guard_terminate": False,
        "qiu_guard_clip_fraction": 0.02,
        "qiu_guard_extinction_count": 10,
        "external_delta_eta_target": 0.02,
    }
    variants = {
        "L0": dict(backend="qiu_full_field", source="legacy_midpoint", coupled=False, dt=0.04),
        "L1": dict(backend="qiu_full_field", source="legacy_midpoint", coupled=False, dt=0.02),
        "K": dict(backend="fft_eigenstrain_v2", source="legacy_midpoint", coupled=False, dt=0.04),
        "S": dict(backend="qiu_full_field", source="local_sweep", coupled=False, dt=0.04),
        "T": dict(backend="qiu_full_field", source="legacy_midpoint", coupled=False, dt=0.04, external=True),
        "KS": dict(backend="fft_eigenstrain_v2", source="local_sweep", coupled=False, dt=0.04),
        "KST": dict(backend="fft_eigenstrain_v2", source="local_sweep", coupled=True, dt=0.04),
    }
    matrix_rows = []
    tables = []
    histories = {}
    for name, variant in variants.items():
        parameters = {
            **common,
            "full_field_source_mapping": variant["source"],
            "coupled_integration_enabled": variant["coupled"],
            "full_field_external_timestep_control": variant.get("external", False),
        }
        config = ModelConfig(
            regime=name, seed=5101, pf=replace(pf, time_step=variant["dt"]),
            mechanics_backend=variant["backend"], compatibility_model="off",
            active_modules=(
                ("fft_eigenstrain_shear",)
                if variant["source"] == "local_sweep"
                else ("qiu_reference_shear",)
            ),
            output_cadence=20, max_steps=args.steps, termination_grains=5,
            parameters=parameters,
        )
        run = args.output / name
        started = time.perf_counter()
        simulation = EventResolvedSimulation(config, run, code_sha=sha)
        simulation.run()
        elapsed = time.perf_counter() - started
        frame = ds.dataset(run / "per_step_diagnostics.parquet", format="parquet").to_table().to_pandas()
        frame.insert(0, "variant", name)
        histories[name] = frame
        tables.append(pa.Table.from_pandas(frame, preserve_index=False))
        capture_path = run / "diagnostic_capture.json"
        capture = json.loads(capture_path.read_text()) if capture_path.exists() else None
        final = frame.iloc[-1]
        matrix_rows.append({
            "variant": name, "operator": variant["backend"],
            "source_mapping": variant["source"],
            "integration": (
                "energy_checked_external_limit" if variant["coupled"]
                else "external_limit" if variant.get("external") else "legacy_explicit"
            ),
            "configured_dt": variant["dt"],
            "steps": int(final["step"]), "physical_time": float(final["time"]),
            "final_grains": int(final["grain_count"]),
            "guard_triggered": capture is not None,
            "first_trigger_step": int(capture["step"]) if capture else "",
            "guard_reasons": ";".join(capture["reasons"]) if capture else "",
            "max_100_step_loss": int(frame["largest_100_step_population_loss"].max()),
            "final_interfacial_energy": float(final["interfacial_energy"]),
            "final_elastic_energy": float(final["elastic_energy"]),
            "final_total_energy": float(final["total_energy"]),
            "max_equilibrium_residual": float(frame["mechanical_equilibrium_residual"].max()),
            "max_abs_source_work_error": float(frame["source_work_error"].abs().max()),
            "max_stress_linf": float(frame["stress_linf"].max()),
            "max_eigenstrain_linf": float(frame["eigenstrain_linf"].max()),
            "max_clipped_fraction": float(frame["clipped_fraction"].max()),
            "max_compactness": float(frame["compactness_max"].max()),
            "wall_seconds": elapsed,
            "run_path": str(run.resolve()),
        })

    combined_dir = args.output / "per_step_diagnostics.parquet"
    combined_dir.mkdir()
    for index, table in enumerate(tables):
        pq.write_table(table, combined_dir / f"part-{index:03d}.parquet", compression="zstd")
    matrix_rows.extend([
        {
            "variant": "R", "operator": "qiu_current_geometry_line",
            "source_mapping": "reference_line", "integration": "reference_explicit",
            "guard_triggered": "", "run_path": "",
            "exclusion": "not implemented; distinct archived model, outside selected Path-B production backend",
        },
        {
            **next(row for row in matrix_rows if row["variant"] == "KST"),
            "variant": "F", "run_path": str((args.output / "KST").resolve()),
            "exclusion": "alias of fully corrected FFT_EIGENSTRAIN_V2; no duplicate trajectory launched",
        },
    ])
    columns = sorted({key for row in matrix_rows for key in row})
    with (args.output / "run_matrix.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(matrix_rows)

    plot_specs = (
        ("grain_count", "N", "population.png"),
        ("stress_linf", "max |stress|", "stress.png"),
        ("total_energy", "interfacial + elastic energy", "energy.png"),
        ("compactness_max", "maximum compactness", "compactness.png"),
    )
    plot_dir = args.output / "plots"; plot_dir.mkdir()
    for column, ylabel, filename in plot_specs:
        fig, axis = plt.subplots(figsize=(7.0, 4.5))
        for name, frame in histories.items():
            axis.plot(frame["time"], frame[column], label=name, linewidth=1.0)
        axis.set_xlabel("physical time"); axis.set_ylabel(ylabel)
        axis.legend(ncol=4, fontsize=8); fig.tight_layout()
        fig.savefig(plot_dir / filename, dpi=180); plt.close(fig)

    result = {
        "schema_version": 1, "source_commit": sha,
        "initial_state": str(initial.resolve()), "initial_state_sha256": sha256(initial),
        "shape": list(pf.shape), "initial_grains": args.grains,
        "configured_steps": args.steps, "variants": matrix_rows,
        "classification_scope": "reduced causal control; not production qualification",
    }
    atomic_write_text(args.output / "matrix_summary.json", json.dumps(result, indent=2) + "\n")
    checksums = {
        str(path.relative_to(args.output)): sha256(path)
        for path in sorted(args.output.rglob("*")) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
