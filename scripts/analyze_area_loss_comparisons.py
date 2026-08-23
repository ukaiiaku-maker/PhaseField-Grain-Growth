#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


WINDOWS = (
    (190, 175), (175, 160), (190, 160), (160, 140),
    (140, 120), (120, 100), (100, 80), (80, 70),
)


def _runs(root: Path) -> list[Path]:
    for name in ("campaign_manifest.json", "video_manifest.json"):
        path = root / name
        if path.exists():
            manifest = json.loads(path.read_text())
            return [Path(item["path"] if isinstance(item, dict) else item) for item in manifest["runs"]]
    return sorted(path.parent for path in root.glob("*/manifest.json"))


def _series(run: Path) -> pd.DataFrame:
    tracks = pd.read_csv(run / "grain_tracks.csv")
    tracks["radius_area"] = np.sqrt(tracks["area"].to_numpy(float) / np.pi)
    return tracks.groupby(["step", "time"], as_index=False).agg(
        grain_count=("grain_id", "nunique"), mean_radius=("radius_area", "mean")
    ).sort_values("step").drop_duplicates("step", keep="last")


def _window_metrics(series: pd.DataFrame, high: int, low: int) -> dict[str, Any]:
    selected = series[(series.grain_count <= high) & (series.grain_count >= low)]
    if len(selected) < 3:
        return {"samples": len(selected)}
    time = selected.time.to_numpy(float)
    radius = selected.mean_radius.to_numpy(float)
    valid = np.diff(time) > 0.0
    speed = np.diff(radius)[valid] / np.diff(time)[valid]
    absolute = np.abs(speed)
    mean_absolute = float(np.mean(absolute)) if len(absolute) else 0.0
    threshold = 0.05 * mean_absolute
    total_motion = float(np.sum(absolute))
    top = max(1, int(np.ceil(0.05 * len(absolute)))) if len(absolute) else 0
    return {
        "samples": len(selected),
        "time_start": float(time[0]), "time_end": float(time[-1]),
        "duration": float(time[-1] - time[0]),
        "N_start": int(selected.grain_count.iloc[0]),
        "N_end": int(selected.grain_count.iloc[-1]),
        "growth_rate": float(np.polyfit(time, radius, 1)[0]),
        "jerkiness_CV": float(np.std(absolute) / mean_absolute) if mean_absolute else 0.0,
        "stationary_fraction": float(np.mean(absolute <= threshold)) if len(absolute) else np.nan,
        "motion_top_5pct": float(np.sort(absolute)[-top:].sum() / total_motion) if total_motion else 0.0,
    }


def _occupancy(run: Path) -> dict[str, float]:
    path = run / "mechanism_state.csv"
    if not path.exists():
        return {"G_occupancy": 0.0, "T_occupancy": 0.0, "C_occupancy": 0.0}
    state = pd.read_csv(path)
    gb = state[state.entity_type == "GB"]
    tj = state[state.entity_type == "TJ"]
    return {
        "G_occupancy": float(gb.G_pending.mean()) if len(gb) else 0.0,
        "T_occupancy": float(tj.T_pending.mean()) if len(tj) else 0.0,
        "C_occupancy": float(state.C_pending.mean()) if len(state) else 0.0,
    }


def _sink(run: Path) -> dict[str, float]:
    path = run / "defect_inventory.csv"
    if not path.exists():
        return {"f_GB": 0.0, "f_TJ": 0.0, "max_conservation_residual": 0.0}
    final = pd.read_csv(path).iloc[-1]
    return {
        "f_GB": float(final.GB_accommodation_fraction),
        "f_TJ": float(final.TJ_accommodation_fraction),
        "max_conservation_residual": float(final.max_abs_conservation_residual),
    }


def collect(campaign: str, root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary: list[dict[str, Any]] = []
    windows: list[dict[str, Any]] = []
    for run in _runs(root):
        manifest = json.loads((run / "manifest.json").read_text())
        regime = str(manifest["config"]["regime"])
        series = _series(run)
        base = {
            "campaign": campaign, "regime": regime, "path": str(run),
            "ending_step": int(series.step.iloc[-1]),
            "ending_time": float(series.time.iloc[-1]),
            "ending_grains": int(series.grain_count.iloc[-1]),
            **_occupancy(run), **_sink(run),
        }
        first_window: dict[str, Any] = {}
        for high, low in WINDOWS:
            metrics = _window_metrics(series, high, low)
            windows.append({
                "campaign": campaign, "regime": regime,
                "window": f"N{high}_to_{low}", "N_high": high, "N_low": low,
                **metrics,
            })
            if (high, low) == (190, 160):
                first_window = metrics
        summary.append({**base, **{f"N190_160_{key}": value for key, value in first_window.items()}})
    return pd.DataFrame(summary), pd.DataFrame(windows)


def _base_without_c(regime: str) -> str:
    stem = regime.split("_", 1)[0].replace("C", "")
    return stem or "B0"


def _legacy_name(regime: str) -> str:
    stem = regime.split("_", 1)[0]
    return "LEGACY_" + stem


def crossover_table(windows: pd.DataFrame) -> pd.DataFrame:
    lookup = windows.set_index(["campaign", "regime", "window"])
    rows: list[dict[str, Any]] = []
    corrected = windows[windows.campaign == "corrected"].regime.unique()
    for regime in corrected:
        for comparison, campaign, reference in (
            ("corrected_vs_no_C", "non_c", _base_without_c(regime)),
            ("corrected_vs_legacy_C", "legacy", _legacy_name(regime)),
        ):
            previous_sign = 0
            for high, low in WINDOWS:
                window = f"N{high}_to_{low}"
                key = ("corrected", regime, window)
                ref_key = (campaign, reference, window)
                if key not in lookup.index or ref_key not in lookup.index:
                    continue
                rate = float(lookup.loc[key].get("growth_rate", np.nan))
                ref_rate = float(lookup.loc[ref_key].get("growth_rate", np.nan))
                difference = rate - ref_rate
                sign = int(np.sign(difference)) if np.isfinite(difference) else 0
                rows.append({
                    "regime": regime, "comparison": comparison, "reference": reference,
                    "window": window, "corrected_rate": rate, "reference_rate": ref_rate,
                    "rate_ratio": rate / ref_rate if np.isfinite(ref_rate) and ref_rate != 0.0 else np.nan,
                    "rate_difference": difference,
                    "kinetic_crossover": bool(previous_sign and sign and sign != previous_sign),
                })
                if sign:
                    previous_sign = sign
    return pd.DataFrame(rows)


def factorial_effects(summary: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "N190_160_growth_rate", "N190_160_jerkiness_CV", "N190_160_stationary_fraction",
        "N190_160_motion_top_5pct", "G_occupancy", "T_occupancy", "C_occupancy", "f_TJ",
    ]
    rows: list[dict[str, Any]] = []
    corrected = summary[summary.campaign == "corrected"].copy()
    for sink_suffix in ("GB", "GBTJ"):
        design = corrected[corrected.regime.str.endswith("_" + sink_suffix)].copy()
        if len(design) != 8:
            continue
        stem = design.regime.str.split("_").str[0]
        code = pd.DataFrame({letter: stem.str.contains(letter).map({False: -1.0, True: 1.0}) for letter in "GTS"})
        g, t, s = code["G"], code["T"], code["S"]
        matrix = np.column_stack([
            np.ones(len(design)), g, t, s,
            g * t, g * s, t * s, g * t * s,
        ])
        terms = ["intercept", "G", "T", "S", "G:T", "G:S", "T:S", "G:T:S"]
        for metric in metrics:
            values = pd.to_numeric(design.get(metric), errors="coerce").to_numpy(float)
            valid = np.isfinite(values)
            if valid.sum() < 8:
                continue
            coefficients, *_ = np.linalg.lstsq(matrix[valid], values[valid], rcond=None)
            for term, coefficient in zip(terms[1:], coefficients[1:], strict=True):
                rows.append({
                    "sink_model": sink_suffix, "metric": metric, "effect": term,
                    "high_minus_low_effect": 2.0 * float(coefficient),
                })
    return pd.DataFrame(rows)


def minimum_models(summary: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "N190_160_growth_rate", "N190_160_jerkiness_CV",
        "N190_160_stationary_fraction", "N190_160_motion_top_5pct", "C_occupancy", "f_TJ",
    ]
    rows: list[pd.DataFrame] = []
    corrected = summary[summary.campaign == "corrected"].copy()
    for suffix in ("GB", "GBTJ"):
        frame = corrected[corrected.regime.str.endswith("_" + suffix)].copy()
        target_name = "GTSC_" + suffix
        if len(frame) != 8 or target_name not in set(frame.regime):
            continue
        values = frame[metrics].to_numpy(float)
        target = frame.loc[frame.regime == target_name, metrics].iloc[0].to_numpy(float)
        scale = np.nanstd(values, axis=0)
        scale[~np.isfinite(scale) | (scale <= np.finfo(float).tiny)] = 1.0
        frame["target_standardized_RMS"] = np.sqrt(np.nanmean(((values - target) / scale) ** 2, axis=1))
        stems = frame.regime.str.split("_").str[0]
        frame["mechanism_count"] = 1 + stems.map(lambda value: sum(letter in value for letter in "GTS"))
        pareto = []
        for row in frame.itertuples():
            dominated = np.any(
                (frame.mechanism_count <= row.mechanism_count)
                & (frame.target_standardized_RMS <= row.target_standardized_RMS)
                & ((frame.mechanism_count < row.mechanism_count)
                   | (frame.target_standardized_RMS < row.target_standardized_RMS))
            )
            pareto.append(not bool(dominated))
        frame["Pareto_optimal"] = pareto
        frame.insert(0, "sink_model", suffix)
        rows.append(frame[[
            "sink_model", "regime", "mechanism_count", "target_standardized_RMS", "Pareto_optimal", *metrics,
        ]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corrected", required=True, type=Path)
    parser.add_argument("--non-c", required=True, type=Path)
    parser.add_argument("--legacy", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summaries, windows = [], []
    for name, root in (("corrected", args.corrected), ("non_c", args.non_c), ("legacy", args.legacy)):
        summary, topology = collect(name, root)
        summaries.append(summary)
        windows.append(topology)
    summary = pd.concat(summaries, ignore_index=True)
    topology = pd.concat(windows, ignore_index=True)
    summary.to_csv(args.output / "comparative_run_metrics.csv", index=False)
    topology.to_csv(args.output / "common_topology_windows.csv", index=False)
    crossover_table(topology).to_csv(args.output / "kinetic_crossovers.csv", index=False)
    factorial_effects(summary).to_csv(args.output / "corrected_factorial_effects.csv", index=False)
    minimum_models(summary).to_csv(args.output / "revised_minimum_models.csv", index=False)
    (args.output / "comparison_manifest.json").write_text(json.dumps({
        "corrected": str(args.corrected), "non_c": str(args.non_c),
        "legacy": str(args.legacy), "topology_windows": WINDOWS,
    }, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
