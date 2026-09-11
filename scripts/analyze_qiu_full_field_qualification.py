#!/usr/bin/env python3
"""Aggregate completed QIU qualification runs into auditable tables and plots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

from grain_growth_pf.io.checkpoints import atomic_write_text


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_run(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("runs must be LABEL=/absolute/or/relative/path")
    label, raw = value.split("=", 1)
    path = Path(raw).resolve()
    if not label or not path.is_dir():
        raise argparse.ArgumentTypeError(f"invalid run specification {value!r}")
    return label, path


def load(label: str, run: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    diagnostic = run / "per_step_diagnostics.parquet"
    if not diagnostic.is_dir():
        raise FileNotFoundError(diagnostic)
    frame = ds.dataset(diagnostic, format="parquet").to_table().to_pandas()
    frame = frame.sort_values(["step", "time"]).drop_duplicates("step", keep="last")
    frame.insert(0, "qualification_run", label)
    manifest = json.loads((run / "manifest.json").read_text())
    return frame, manifest


def relative_series(frame: pd.DataFrame) -> pd.Series:
    denominator = (
        frame["source_elastic_energy_change"].abs()
        + frame["predicted_source_work"].abs()
        + np.finfo(float).eps
    )
    return frame["source_work_error"].abs() / denominator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    frames: dict[str, pd.DataFrame] = {}
    manifests: dict[str, dict[str, object]] = {}
    matrix: list[dict[str, object]] = []
    for label, path in args.run:
        frame, manifest = load(label, path)
        frames[label], manifests[label] = frame, manifest
        work_relative = relative_series(frame)
        energy = np.asarray(frame["total_energy"], dtype=float)
        scale = np.maximum(np.abs(energy[:-1]), 1.0)
        increases = np.diff(energy) > (1e-10 + 1e-10 * scale)
        capture = path / "diagnostic_capture.json"
        final = frame.iloc[-1]
        matrix.append({
            "label": label, "run_path": str(path),
            "source_commit": manifest.get("git_sha", ""),
            "manifest_status": manifest.get("status", ""),
            "terminal_step": int(final["step"]), "terminal_time": float(final["time"]),
            "terminal_grains": int(final["grain_count"]),
            "terminal_G": float(final["G_population"]),
            "guard_triggered": capture.exists(),
            "guard": capture.read_text().strip() if capture.exists() else "",
            "maximum_equilibrium_residual": float(frame["mechanical_equilibrium_residual"].max()),
            "maximum_relative_work_residual": float(work_relative.max()),
            "p95_relative_work_residual": float(work_relative.quantile(0.95)),
            "complete_energy_increase_count": int(np.count_nonzero(increases)),
            "maximum_stress_linf": float(frame["stress_linf"].max()),
            "maximum_eigenstrain_linf": float(frame["eigenstrain_linf"].max()),
            "maximum_clipped_fraction": float(frame["clipped_fraction"].max()),
            "maximum_one_step_extinctions": int(frame["newly_extinct_phases"].max()),
            "maximum_100_step_loss": int(frame["largest_100_step_population_loss"].max()),
            "maximum_compactness": float(frame["compactness_max"].max()),
            "maximum_aspect_ratio": float(frame["aspect_ratio_max"].max()),
        })

    columns = sorted({key for row in matrix for key in row})
    with (args.output / "run_matrix.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(matrix)
    combined = args.output / "per_step_diagnostics.parquet"; combined.mkdir()
    for index, (label, frame) in enumerate(frames.items()):
        pq.write_table(
            pa.Table.from_pandas(frame, preserve_index=False),
            combined / f"part-{index:03d}-{label}.parquet", compression="zstd",
        )

    plot_dir = args.output / "plots"; plot_dir.mkdir()
    plot_specs = (
        (("grain_count",), "grain count N", "N_vs_time.png", False),
        (("G_population",), "population grain size G", "G_vs_time.png", False),
        (("interfacial_energy", "elastic_energy", "total_energy"), "energy", "energy.png", False),
        (("stress_l2", "stress_linf", "stress_p95"), "stress", "stress_norms.png", True),
        (("eigenstrain_l2", "eigenstrain_linf"), "eigenstrain", "eigenstrain_norms.png", True),
        (("elastic_capillary_ratio_p95", "elastic_capillary_ratio_max"), "|elastic| / |capillary|", "elastic_capillary_ratio.png", True),
        (("clipped_fraction", "newly_extinct_phases"), "clipping / extinctions", "clipping_extinctions.png", False),
        (("compactness_mean", "compactness_p95", "aspect_ratio_mean", "aspect_ratio_p95"), "morphology", "morphology.png", False),
        (("source_work_error",), "source work residual", "work_conjugacy_residual.png", True),
        (("mechanical_equilibrium_residual",), "equilibrium residual", "mechanical_equilibrium_residual.png", True),
        (("used_dt", "external_dt_limit"), "timestep", "timestep.png", True),
    )
    for fields, ylabel, filename, logarithmic in plot_specs:
        figure, axis = plt.subplots(figsize=(7.4, 4.6))
        positive = True
        for label, frame in frames.items():
            for field in fields:
                if field not in frame:
                    continue
                values = np.asarray(frame[field], dtype=float)
                if field == "source_work_error":
                    values = np.abs(values)
                finite = np.isfinite(values)
                positive &= bool(np.all(values[finite] > 0))
                axis.plot(frame.loc[finite, "time"], values[finite], label=f"{label}:{field}", linewidth=1.0)
        if logarithmic and positive:
            axis.set_yscale("log")
        axis.set_xlabel("physical time"); axis.set_ylabel(ylabel)
        axis.legend(fontsize=7, ncol=2); figure.tight_layout()
        figure.savefig(plot_dir / filename, dpi=180); plt.close(figure)

    timestep: dict[str, object] = {"available": False}
    if "corrected" in frames and "refined" in frames:
        base, fine = frames["corrected"], frames["refined"]
        matched_time = min(float(base["time"].max()), float(fine["time"].max()))
        observables = (
            "G_population", "interfacial_energy", "elastic_energy",
            "compactness_mean", "stress_p95",
        )
        values = {
            label: {
                field: float(np.interp(matched_time, frame["time"], frame[field]))
                for field in observables
            }
            for label, frame in (("corrected", base), ("refined", fine))
        }
        relative = {
            field: abs(values["corrected"][field] - values["refined"][field]) /
            max(abs(values["refined"][field]), np.finfo(float).eps)
            for field in observables
        }
        timestep = {
            "available": True, "matched_time": matched_time,
            "observables": values, "relative_differences": relative,
            "gates": {
                "G_within_1pct": relative["G_population"] <= 0.01,
                "interfacial_energy_within_2pct": relative["interfacial_energy"] <= 0.02,
                "elastic_energy_within_2pct": relative["elastic_energy"] <= 0.02,
                "compactness_within_5pct": relative["compactness_mean"] <= 0.05,
                "stress_p95_within_5pct": relative["stress_p95"] <= 0.05,
            },
        }

    summary = {
        "schema_version": 1, "runs": matrix, "timestep_comparison": timestep,
        "plots": {path.name: str(path.resolve()) for path in sorted(plot_dir.glob("*.png"))},
        "scientific_status": (
            "complete_input_set" if timestep["available"] and len(frames) >= 4
            else "partial_qualification_inputs"
        ),
    }
    atomic_write_text(args.output / "qualification_analysis_summary.json", json.dumps(summary, indent=2) + "\n")
    checksums = {
        str(path.relative_to(args.output)): sha256(path)
        for path in sorted(args.output.rglob("*")) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
