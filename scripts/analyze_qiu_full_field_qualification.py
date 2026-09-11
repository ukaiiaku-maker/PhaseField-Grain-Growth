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


def json_safe(value: object) -> object:
    """Replace non-finite numeric values with JSON null recursively."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


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
    attestation = run / "source_attestation.json"
    if attestation.exists():
        verified = json.loads(attestation.read_text())
        manifest["source_attestation"] = verified
        manifest["verified_source_commit"] = verified["verified_source_commit"]
    return frame, manifest


def capture_assessment(run: Path) -> dict[str, object]:
    """Classify a raw guard capture without erasing its provenance."""
    capture = run / "diagnostic_capture.json"
    assessment_path = run / "diagnostic_capture_assessment.json"
    result: dict[str, object] = {
        "raw_capture_exists": capture.exists(),
        "raw_capture_path": str(capture.resolve()) if capture.exists() else "",
        "raw_capture_sha256": sha256(capture) if capture.exists() else "",
        "raw_capture_step": None,
        "accepted_as_transition": capture.exists(),
        "assessment_path": "",
        "assessment_sha256": "",
        "assessment_reason": "",
        "scientific_transition_step": None,
    }
    if capture.exists():
        raw = json.loads(capture.read_text())
        raw_step = raw.get("step")
        if isinstance(raw_step, int) and not isinstance(raw_step, bool):
            result["raw_capture_step"] = raw_step
    if not assessment_path.exists():
        return result
    assessment = json.loads(assessment_path.read_text())
    expected = str(assessment.get("raw_capture_sha256", ""))
    if not capture.exists():
        raise ValueError(f"{assessment_path} assesses a missing raw capture")
    if expected != result["raw_capture_sha256"]:
        raise ValueError(f"{assessment_path} raw capture SHA-256 mismatch")
    if not str(assessment.get("reason", "")).strip():
        raise ValueError(f"{assessment_path} requires a nonempty reason")
    if not isinstance(assessment.get("accepted_as_transition"), bool):
        raise ValueError(f"{assessment_path} requires boolean accepted_as_transition")
    transition_step = assessment.get("scientific_transition_step")
    if transition_step is not None and (
        not isinstance(transition_step, int) or isinstance(transition_step, bool)
    ):
        raise ValueError(f"{assessment_path} scientific_transition_step must be an integer or null")
    result.update({
        "accepted_as_transition": assessment["accepted_as_transition"],
        "assessment_path": str(assessment_path.resolve()),
        "assessment_sha256": sha256(assessment_path),
        "assessment_reason": str(assessment["reason"]),
        "scientific_transition_step": transition_step,
    })
    return result


def _first_step(frame: pd.DataFrame, mask: np.ndarray) -> int | None:
    matches = frame.loc[np.asarray(mask, dtype=bool), "step"]
    return int(matches.iloc[0]) if not matches.empty else None


def transition_timing(frame: pd.DataFrame) -> dict[str, object]:
    """Locate each objective avalanche indicator without conflating them."""
    ordered = frame.sort_values("step").reset_index(drop=True)
    count = np.asarray(ordered["grain_count"], dtype=float)
    loss100 = np.asarray(ordered["largest_100_step_population_loss"], dtype=float)
    preceding_count = count + loss100
    burst_fraction = np.divide(
        loss100,
        preceding_count,
        out=np.zeros_like(loss100),
        where=preceding_count > 0,
    )
    energy = np.asarray(ordered["total_energy"], dtype=float)
    energy_scale = np.maximum(np.abs(energy[:-1]), 1.0)
    increases = np.zeros(len(ordered), dtype=bool)
    increases[1:] = np.diff(energy) > (1e-10 + 1e-10 * energy_scale)
    nonfinite = np.zeros(len(ordered), dtype=bool)
    for field in ("total_energy", "stress_linf", "eigenstrain_linf"):
        nonfinite |= ~np.isfinite(np.asarray(ordered[field], dtype=float))
    population_change = np.zeros(len(ordered), dtype=bool)
    population_change[1:] = np.diff(count) < 0
    timing = {
        "first_population_decrease_step": _first_step(ordered, population_change),
        "first_extinction_step": _first_step(
            ordered, np.asarray(ordered["newly_extinct_phases"], dtype=float) > 0
        ),
        "first_100_step_population_loss_above_10pct_step": _first_step(
            ordered, burst_fraction > 0.10
        ),
        "first_mean_compactness_above_2p5_step": _first_step(
            ordered, np.asarray(ordered["compactness_mean"], dtype=float) > 2.5
        ),
        "first_max_compactness_above_6_step": _first_step(
            ordered, np.asarray(ordered["compactness_max"], dtype=float) > 6.0
        ),
        "first_disconnected_grain_step": _first_step(
            ordered, np.asarray(ordered["disconnected_grain_count"], dtype=float) > 0
        ),
        "first_complete_energy_increase_step": _first_step(ordered, increases),
        "first_nonfinite_critical_field_step": _first_step(ordered, nonfinite),
        "maximum_100_step_population_loss_fraction": float(np.max(burst_fraction)),
    }
    candidates = [
        timing[key] for key in (
            "first_100_step_population_loss_above_10pct_step",
            "first_mean_compactness_above_2p5_step",
            "first_max_compactness_above_6_step",
            "first_nonfinite_critical_field_step",
        ) if timing[key] is not None
    ]
    timing["first_objective_avalanche_indicator_step"] = (
        min(candidates) if candidates else None
    )
    return timing


def pre_extinction_precursor(frame: pd.DataFrame) -> dict[str, object]:
    """Summarize feedback changes that precede the first extinction."""
    ordered = frame.sort_values("step").reset_index(drop=True)
    extinct = np.asarray(ordered["newly_extinct_phases"], dtype=float) > 0
    first_extinction = _first_step(ordered, extinct)
    pre = ordered if first_extinction is None else ordered.loc[ordered["step"] < first_extinction]
    if pre.empty:
        return {"available": False, "first_extinction_step": first_extinction}
    window = min(20, max(1, len(pre) // 2))
    first = pre.iloc[:window]
    last = pre.iloc[-window:]
    result: dict[str, object] = {
        "available": True,
        "first_extinction_step": first_extinction,
        "step_first": int(pre.iloc[0]["step"]),
        "step_last": int(pre.iloc[-1]["step"]),
        "grain_count_change": int(pre.iloc[-1]["grain_count"] - pre.iloc[0]["grain_count"]),
        "maximum_disconnected_grains": int(pre["disconnected_grain_count"].max()),
        "maximum_compactness": float(pre["compactness_max"].max()),
    }
    for field in (
        "source_increment_l2", "stress_linf", "eigenstrain_linf",
        "interfacial_energy", "elastic_energy", "total_energy",
    ):
        initial = float(first[field].median())
        final = float(last[field].median())
        result[f"{field}_initial_20_step_median"] = initial
        result[f"{field}_final_20_step_median"] = final
        result[f"{field}_median_change"] = final - initial
    return result


def relative_series(frame: pd.DataFrame) -> pd.Series:
    denominator = (
        frame["source_elastic_energy_change"].abs()
        + frame["predicted_source_work"].abs()
        + np.finfo(float).eps
    )
    return frame["source_work_error"].abs() / denominator


def interpolate(frame: pd.DataFrame, coordinate: str, fields: tuple[str, ...], value: float) -> dict[str, float]:
    """Interpolate observables after collapsing repeated event coordinates."""
    selected = frame[[coordinate, *fields]].replace([np.inf, -np.inf], np.nan).dropna()
    selected = selected.groupby(coordinate, as_index=False).last().sort_values(coordinate)
    x = np.asarray(selected[coordinate], dtype=float)
    return {
        field: float(np.interp(value, x, np.asarray(selected[field], dtype=float)))
        for field in fields
    }


def relative_values(first: dict[str, float], second: dict[str, float]) -> dict[str, float]:
    return {
        field: abs(first[field] - second[field]) /
        max(abs(second[field]), np.finfo(float).eps)
        for field in first
    }


def avalanche_indicator(frame: pd.DataFrame, capture_exists: bool = False) -> dict[str, object]:
    initial = max(float(frame.iloc[0]["grain_count"]), 1.0)
    burst_fraction = float(frame["largest_100_step_population_loss"].max()) / initial
    nonfinite = any(
        field in frame and not np.all(np.isfinite(np.asarray(frame[field], dtype=float)))
        for field in ("total_energy", "stress_linf", "eigenstrain_linf")
    )
    triggered = bool(
        capture_exists or burst_fraction > 0.10 or
        frame["compactness_mean"].max() > 2.5 or
        frame["compactness_max"].max() > 6.0 or nonfinite
    )
    return {
        "detected": triggered,
        "diagnostic_capture_exists": capture_exists,
        "nonfinite_fields": nonfinite,
        "maximum_100_step_loss_fraction": burst_fraction,
    }


def compare_timestep_runs(
    base: pd.DataFrame, fine: pd.DataFrame, *,
    base_capture: bool = False, fine_capture: bool = False,
) -> dict[str, object]:
    observables = (
        "G_population", "interfacial_energy", "elastic_energy",
        "compactness_mean", "compactness_p95", "stress_p95",
    )
    matched_time = min(float(base["time"].max()), float(fine["time"].max()))
    by_time = {
        label: interpolate(frame, "time", observables, matched_time)
        for label, frame in (("corrected", base), ("refined", fine))
    }
    time_relative = relative_values(by_time["corrected"], by_time["refined"])

    matched_g = min(float(base["G_population"].max()), float(fine["G_population"].max()))
    progress_fields = ("time",) + tuple(field for field in observables if field != "G_population")
    by_progress = {
        label: interpolate(frame, "G_population", progress_fields, matched_g)
        for label, frame in (("corrected", base), ("refined", fine))
    }
    progress_relative = relative_values(by_progress["corrected"], by_progress["refined"])

    overlap_start = max(float(base["time"].min()), float(fine["time"].min()))
    overlap_end = matched_time
    grid = np.unique(np.concatenate((
        np.asarray(base.loc[base["time"].between(overlap_start, overlap_end), "time"], dtype=float),
        np.asarray(fine.loc[fine["time"].between(overlap_start, overlap_end), "time"], dtype=float),
    )))
    history: dict[str, dict[str, float]] = {}
    for field in ("interfacial_energy", "elastic_energy"):
        base_series = base[["time", field]].replace([np.inf, -np.inf], np.nan).dropna()
        base_series = base_series.groupby("time", as_index=False).last().sort_values("time")
        fine_series = fine[["time", field]].replace([np.inf, -np.inf], np.nan).dropna()
        fine_series = fine_series.groupby("time", as_index=False).last().sort_values("time")
        left = np.interp(grid, base_series["time"], base_series[field])
        right = np.interp(grid, fine_series["time"], fine_series[field])
        scale = max(float(np.max(np.abs(right))), np.finfo(float).eps)
        normalized = np.abs(left - right) / scale
        history[field] = {
            "normalized_linf": float(np.max(normalized)),
            "normalized_p95": float(np.quantile(normalized, 0.95)),
        }

    base_avalanche = avalanche_indicator(base, base_capture)
    fine_avalanche = avalanche_indicator(fine, fine_capture)
    gates = {
        "matched_time_G_within_1pct": time_relative["G_population"] <= 0.01,
        "matched_time_interfacial_energy_within_2pct": time_relative["interfacial_energy"] <= 0.02,
        "matched_time_elastic_energy_within_2pct": time_relative["elastic_energy"] <= 0.02,
        "matched_time_compactness_within_5pct": max(
            time_relative["compactness_mean"], time_relative["compactness_p95"]
        ) <= 0.05,
        "matched_time_stress_p95_within_5pct": time_relative["stress_p95"] <= 0.05,
        "matched_progress_interfacial_energy_within_2pct": progress_relative["interfacial_energy"] <= 0.02,
        "matched_progress_elastic_energy_within_2pct": progress_relative["elastic_energy"] <= 0.02,
        "matched_progress_compactness_within_5pct": max(
            progress_relative["compactness_mean"], progress_relative["compactness_p95"]
        ) <= 0.05,
        "matched_progress_stress_p95_within_5pct": progress_relative["stress_p95"] <= 0.05,
        "interfacial_energy_history_within_2pct": history["interfacial_energy"]["normalized_linf"] <= 0.02,
        "elastic_energy_history_within_2pct": history["elastic_energy"]["normalized_linf"] <= 0.02,
        "same_avalanche_classification": base_avalanche["detected"] == fine_avalanche["detected"],
    }
    return {
        "available": True,
        "matched_time": {"time": matched_time, "observables": by_time, "relative_differences": time_relative},
        "matched_grain_size_progress": {
            "G_population": matched_g, "observables": by_progress,
            "relative_differences": progress_relative,
        },
        "interpolated_energy_history": history,
        "avalanche_indicators": {"corrected": base_avalanche, "refined": fine_avalanche},
        "gates": gates,
        "all_gates_pass": all(gates.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    frames: dict[str, pd.DataFrame] = {}
    manifests: dict[str, dict[str, object]] = {}
    captures: dict[str, bool] = {}
    capture_assessments: dict[str, dict[str, object]] = {}
    transition_timings: dict[str, dict[str, object]] = {}
    precursors: dict[str, dict[str, object]] = {}
    matrix: list[dict[str, object]] = []
    for label, path in args.run:
        frame, manifest = load(label, path)
        frames[label], manifests[label] = frame, manifest
        work_relative = relative_series(frame)
        energy = np.asarray(frame["total_energy"], dtype=float)
        scale = np.maximum(np.abs(energy[:-1]), 1.0)
        increases = np.diff(energy) > (1e-10 + 1e-10 * scale)
        capture = path / "diagnostic_capture.json"
        assessment = capture_assessment(path)
        capture_assessments[label] = assessment
        captures[label] = bool(assessment["accepted_as_transition"])
        timing = transition_timing(frame)
        accepted_capture_step = (
            assessment["raw_capture_step"] if captures[label] else None
        )
        explicit_transition_step = assessment["scientific_transition_step"]
        candidates = [
            value for value in (
                timing["first_objective_avalanche_indicator_step"],
                accepted_capture_step,
                explicit_transition_step,
            ) if value is not None
        ]
        timing["selected_transition_step"] = min(candidates) if candidates else None
        transition_timings[label] = timing
        precursors[label] = pre_extinction_precursor(frame)
        final = frame.iloc[-1]
        internal_source = str(manifest.get("git_sha", ""))
        verified_source = str(manifest.get("verified_source_commit", internal_source))
        matrix.append({
            "label": label, "run_path": str(path),
            "source_commit": verified_source,
            "internal_manifest_source_commit": internal_source,
            "source_attestation": str(path / "source_attestation.json")
            if "source_attestation" in manifest else "",
            "manifest_status": manifest.get("status", ""),
            "terminal_step": int(final["step"]), "terminal_time": float(final["time"]),
            "terminal_grains": int(final["grain_count"]),
            "terminal_G": float(final["G_population"]),
            "guard_triggered": capture.exists(),
            "guard_accepted_as_transition": captures[label],
            "guard_assessment_reason": assessment["assessment_reason"],
            "scientific_transition_step": assessment["scientific_transition_step"],
            "selected_transition_step": timing["selected_transition_step"],
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
        timestep = compare_timestep_runs(
            frames["corrected"], frames["refined"],
            base_capture=captures["corrected"], fine_capture=captures["refined"],
        )

        figure, axes = plt.subplots(2, 3, figsize=(12.0, 7.0), sharex=True)
        for axis, field in zip(axes.flat, (
            "grain_count", "G_population", "interfacial_energy",
            "elastic_energy", "compactness_p95", "stress_p95",
        )):
            for label in ("corrected", "refined"):
                axis.plot(frames[label]["time"], frames[label][field], label=label, linewidth=1.0)
            axis.set_title(field); axis.set_xlabel("physical time")
        axes[0, 0].legend(); figure.tight_layout()
        figure.savefig(plot_dir / "timestep_convergence_overlay.png", dpi=180)
        plt.close(figure)

    if "legacy" in frames and "corrected" in frames:
        figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
        for axis, field in zip(axes, ("grain_count", "G_population")):
            for label in ("legacy", "corrected"):
                axis.plot(frames[label]["time"], frames[label][field], label=label, linewidth=1.0)
            axis.set_xlabel("physical time"); axis.set_ylabel(field); axis.legend()
        figure.tight_layout(); figure.savefig(plot_dir / "legacy_vs_corrected.png", dpi=180)
        plt.close(figure)

    seed_labels = [label for label in frames if label == "corrected" or label.startswith("seed")]
    if len(seed_labels) >= 2:
        figure, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)
        for axis, field in zip(axes.flat, (
            "grain_count", "G_population", "compactness_p95", "stress_p95",
        )):
            for label in seed_labels:
                axis.plot(frames[label]["time"], frames[label][field], label=label, linewidth=1.0)
            axis.set_title(field); axis.set_xlabel("physical time")
        axes[0, 0].legend(); figure.tight_layout()
        figure.savefig(plot_dir / "seed_comparison.png", dpi=180); plt.close(figure)

    summary = {
        "schema_version": 1, "runs": matrix, "timestep_comparison": timestep,
        "capture_assessments": capture_assessments,
        "transition_timings": transition_timings,
        "pre_extinction_precursors": precursors,
        "plots": {path.name: str(path.resolve()) for path in sorted(plot_dir.glob("*.png"))},
        "scientific_status": (
            "complete_input_set" if timestep["available"] and len(frames) >= 4
            else "partial_qualification_inputs"
        ),
    }
    atomic_write_text(
        args.output / "qualification_analysis_summary.json",
        json.dumps(json_safe(summary), indent=2, allow_nan=False) + "\n",
    )
    checksums = {
        str(path.relative_to(args.output)): sha256(path)
        for path in sorted(args.output.rglob("*")) if path.is_file()
    }
    atomic_write_text(args.output / "checksums.json", json.dumps(checksums, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
