#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from grain_growth_pf.analysis.campaign import analyze_run
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks
from grain_growth_pf.analysis.growth_law import (
    fit_growth_law,
    fit_growth_law_fixed_exponent,
)
from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger


def _runs(root: Path) -> list[Path]:
    video = root / "video_manifest.json"
    campaign = root / "campaign_manifest.json"
    if video.exists():
        manifest = json.loads(video.read_text())
        return [
            Path(item["path"]) for item in manifest["runs"]
            if item.get("status") == "completed"
        ]
    if campaign.exists():
        manifest = json.loads(campaign.read_text())
        return [Path(item) for item in manifest["runs"]]
    raise ValueError(f"no video_manifest.json or campaign_manifest.json in {root}")


def _window(radius: pd.DataFrame) -> pd.DataFrame:
    initial = float(radius["grain_count"].iloc[0])
    lower = max(float(radius["grain_count"].iloc[-1]), 0.5 * initial)
    selected = radius[
        (radius["grain_count"] <= 0.95 * initial)
        & (radius["grain_count"] >= lower)
    ]
    if len(selected) < 8:
        selected = radius.iloc[max(0, int(0.1 * len(radius))):]
    return selected


def _series_fit(time: np.ndarray, radius: np.ndarray) -> tuple[float, float, float]:
    dt = time - time[0]
    r0 = radius[0]
    design = np.column_stack((radius - r0, radius**2 - r0**2))
    coefficients, _ = nnls(design, dt)
    prediction = design @ coefficients
    sst = float(np.sum((dt - dt.mean()) ** 2))
    r2 = 1.0 - float(np.sum((dt - prediction) ** 2)) / sst if sst else 1.0
    kt = 1.0 / coefficients[0] if coefficients[0] > 0 else np.inf
    kg = 1.0 / coefficients[1] if coefficients[1] > 0 else np.inf
    return kt, kg, r2


def _frame_shear(run: Path) -> dict[str, float]:
    max_state = max_stress = max_qiu = 0.0
    max_free_volume = 0.0
    rms_state = []
    boundary_rms = []
    stored_energy = []
    nonzero_fraction = []
    free_volume_fraction = []
    frame_count = 0
    for frame in sorted((run / "frames").glob("frame-*.npz")):
        with np.load(frame) as data:
            frame_count += 1
            if "shear_state_max_abs" in data:
                max_state = max(max_state, float(data["shear_state_max_abs"]))
                rms_state.append(float(data["shear_state_rms"]))
            elif "shear" in data:
                shear = np.asarray(data["shear"], dtype=float)
                max_state = max(max_state, float(np.max(np.abs(shear))))
                rms_state.append(float(np.sqrt(np.mean(shear**2))))
            if "boundary_shear_state_rms" in data:
                boundary_rms.append(float(data["boundary_shear_state_rms"]))
            elif "shear" in data:
                shear = np.asarray(data["shear"], dtype=float)
                if "boundary_mask" in data:
                    values = shear[np.asarray(data["boundary_mask"], dtype=bool)]
                else:
                    labels = np.asarray(data["labels"])
                    mask = (
                        (labels != np.roll(labels, 1, axis=0))
                        | (labels != np.roll(labels, 1, axis=1))
                    )
                    values = shear[mask]
                boundary_rms.append(
                    float(np.sqrt(np.mean(values**2))) if values.size else 0.0
                )
                nonzero_fraction.append(
                    float(np.mean(np.abs(values) > 1e-12)) if values.size else 0.0
                )
            if "nonzero_gb_length_fraction" in data:
                if len(nonzero_fraction) < frame_count:
                    nonzero_fraction.append(float(data["nonzero_gb_length_fraction"]))
                else:
                    nonzero_fraction[-1] = float(data["nonzero_gb_length_fraction"])
            if "stored_shear_energy" in data:
                stored_energy.append(float(data["stored_shear_energy"]))
            if "free_volume" in data:
                free_volume = np.asarray(data["free_volume"], dtype=float)
                if "boundary_mask" in data:
                    free_volume = free_volume[
                        np.asarray(data["boundary_mask"], dtype=bool)
                    ]
                if free_volume.size:
                    max_free_volume = max(
                        max_free_volume, float(np.max(np.abs(free_volume)))
                    )
                    free_volume_fraction.append(
                        float(np.mean(np.abs(free_volume) > 1e-12))
                    )
            if "shear_stress_max_abs" in data:
                max_stress = max(max_stress, float(data["shear_stress_max_abs"]))
            if "qiu_shear_stress_max_abs" in data:
                max_qiu = max(max_qiu, float(data["qiu_shear_stress_max_abs"]))

    boundary_path = run / "boundary_tracks.csv"
    boundary_max = np.nan
    if boundary_path.exists():
        boundary = pd.read_csv(boundary_path, usecols=["resolved_shear"])
        if len(boundary):
            boundary_max = float(np.nanmax(np.abs(boundary["resolved_shear"].to_numpy(float))))
    return {
        "frame_count": frame_count,
        "max_abs_shear_state": max_state if frame_count else np.nan,
        "mean_frame_shear_state_rms": float(np.mean(rms_state)) if rms_state else np.nan,
        "max_boundary_shear_state_rms": (
            float(np.max(boundary_rms)) if boundary_rms else np.nan
        ),
        "max_stored_shear_energy": (
            float(np.max(stored_energy)) if stored_energy else np.nan
        ),
        "max_nonzero_gb_length_fraction": (
            float(np.max(nonzero_fraction)) if nonzero_fraction else np.nan
        ),
        "final_nonzero_gb_length_fraction": (
            float(nonzero_fraction[-1]) if nonzero_fraction else np.nan
        ),
        "max_abs_free_volume_deficit": max_free_volume if frame_count else np.nan,
        "max_nonzero_free_volume_gb_fraction": (
            float(np.max(free_volume_fraction)) if free_volume_fraction else np.nan
        ),
        "final_nonzero_free_volume_gb_fraction": (
            float(free_volume_fraction[-1]) if free_volume_fraction else np.nan
        ),
        "max_abs_local_shear_stress": (
            max(max_stress, boundary_max) if np.isfinite(boundary_max) else max_stress
        ),
        "max_abs_qiu_shear_stress": max_qiu if frame_count else np.nan,
        "max_abs_boundary_resolved_shear": boundary_max,
    }


def _activation_work(run: Path) -> dict[str, float]:
    path = run / "activation_work.csv"
    empty = {
        "activation_rows": 0,
        "mean_abs_work_capillary": np.nan,
        "mean_abs_work_shear": np.nan,
        "mean_abs_work_free_volume": np.nan,
        "p90_abs_work_shear": np.nan,
        "p90_abs_work_free_volume": np.nan,
        "tj_activation_rows": 0,
        "mean_abs_tj_work_shear": np.nan,
        "mean_abs_tj_work_free_volume": np.nan,
    }
    if not path.exists() or path.stat().st_size == 0:
        return empty
    frame = pd.read_csv(path)
    if frame.empty:
        return empty
    cap = np.abs(frame["work_capillary"].to_numpy(float))
    shear = np.abs(frame["work_shear"].to_numpy(float))
    free = np.abs(frame["work_free_volume"].to_numpy(float))
    tj = frame["event_type"].eq("tj_compatibility_release").to_numpy(bool)
    return {
        "activation_rows": len(frame),
        "mean_abs_work_capillary": float(np.mean(cap)),
        "mean_abs_work_shear": float(np.mean(shear)),
        "mean_abs_work_free_volume": float(np.mean(free)),
        "p90_abs_work_shear": float(np.quantile(shear, 0.90)),
        "p90_abs_work_free_volume": float(np.quantile(free, 0.90)),
        "tj_activation_rows": int(np.count_nonzero(tj)),
        "mean_abs_tj_work_shear": (
            float(np.mean(shear[tj])) if np.any(tj) else np.nan
        ),
        "mean_abs_tj_work_free_volume": (
            float(np.mean(free[tj])) if np.any(tj) else np.nan
        ),
    }


def _pin_durations(run: Path) -> dict[str, float]:
    empty = {
        "pin_episode_count": 0,
        "median_pin_duration": np.nan,
        "p90_pin_duration": np.nan,
        "max_pin_duration": np.nan,
    }
    path = run / "boundary_tracks.csv"
    if not path.exists() or path.stat().st_size == 0:
        return empty
    frame = pd.read_csv(path, usecols=["time", "entity_id", "blocked"])
    if frame.empty:
        return empty
    durations: list[float] = []
    for _, group in frame.groupby("entity_id", sort=False):
        group = group.sort_values("time")
        time = group["time"].to_numpy(float)
        blocked = group["blocked"].to_numpy(bool)
        if not np.any(blocked):
            continue
        gaps = np.diff(time)
        cadence = float(np.median(gaps[gaps > 0])) if np.any(gaps > 0) else 0.0
        start = None
        for index, value in enumerate(blocked):
            if value and start is None:
                start = index
            if start is not None and (not value or index == len(blocked) - 1):
                end = index - 1 if not value else index
                durations.append(max(0.0, time[end] - time[start] + cadence))
                start = None
    if not durations:
        return empty
    values = np.asarray(durations, dtype=float)
    return {
        "pin_episode_count": len(values),
        "median_pin_duration": float(np.median(values)),
        "p90_pin_duration": float(np.quantile(values, 0.90)),
        "max_pin_duration": float(np.max(values)),
    }


def _shear_release_audit(run: Path) -> dict[str, float]:
    empty = {
        "shear_release_rows": 0,
        "fraction_release_ge_prestate": np.nan,
        "fraction_post_release_state_zero": np.nan,
        "median_abs_shear_state_before_release": np.nan,
        "p90_abs_shear_state_before_release": np.nan,
        "median_shear_release": np.nan,
    }
    work_path = run / "activation_work.csv"
    ledger_path = event_ledger_path(run)
    if not work_path.exists() or not ledger_path.exists():
        return empty
    work = pd.read_csv(work_path)
    if work.empty:
        return empty
    events = read_event_ledger(
        ledger_path,
        columns=(
            "time", "step", "event_type", "entity_id", "shear_state_s",
            "release_Delta_s",
        ),
    )
    if events.empty:
        return empty
    # CSV and Parquet readers can differ by a few ulps when parsing the same
    # event time, so use a sub-picosecond matching key rather than raw floats.
    work["_time_key"] = work["time"].round(12)
    events["_time_key"] = events["time"].round(12)
    keys = ["_time_key", "step", "event_type", "entity_id"]
    work["_occurrence"] = work.groupby(keys, dropna=False).cumcount()
    events["_occurrence"] = events.groupby(keys, dropna=False).cumcount()
    joined = work.merge(events, on=keys + ["_occurrence"], how="inner")
    # TJ compatibility events can be driven by adjacent shear without
    # releasing any GB shear memory.  Their ledger shear state is intentionally
    # absent, so exclude them instead of interpreting NaN as a zeroed state.
    joined = joined[joined["shear_state_s"].notna()]
    if joined.empty:
        return empty
    before = np.abs(joined["shear_state_before_release"].to_numpy(float))
    after = np.abs(joined["shear_state_s"].to_numpy(float))
    release = np.abs(joined["release_Delta_s"].fillna(0.0).to_numpy(float))
    active = before > 1e-12
    if not np.any(active):
        return {**empty, "shear_release_rows": len(joined)}
    return {
        "shear_release_rows": len(joined),
        "fraction_release_ge_prestate": float(
            np.mean(release[active] >= before[active] - 1e-12)
        ),
        "fraction_post_release_state_zero": float(
            np.mean(after[active] <= 1e-12)
        ),
        "median_abs_shear_state_before_release": float(np.median(before[active])),
        "p90_abs_shear_state_before_release": float(np.quantile(before[active], 0.90)),
        "median_shear_release": float(np.median(release[active])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    parser.add_argument(
        "--output",
        default="results/production_summaries/jerky_mechanism_integrity.csv",
    )
    args = parser.parse_args()

    rows = []
    for run in _runs(Path(args.root)):
        manifest = json.loads((run / "manifest.json").read_text())
        regime = manifest["config"]["regime"]
        radius = ensemble_radius(load_tracks(run / "grain_tracks.csv"))
        selected = _window(radius)
        time = selected["time"].to_numpy(float)
        r = selected["R_A"].to_numpy(float)
        fit1 = fit_growth_law_fixed_exponent(time, r, 1.0, transient_fraction=0.0)
        fit2 = fit_growth_law_fixed_exponent(time, r, 2.0, transient_fraction=0.0)
        free = fit_growth_law(time, r, transient_fraction=0.0)
        kt_series, kg_series, series_r2 = _series_fit(time, r)
        standard = analyze_run(run)
        rows.append({
            "regime": regime,
            "N_start": int(selected["grain_count"].iloc[0]),
            "N_end": int(selected["grain_count"].iloc[-1]),
            "K_linear_n1": fit1.coefficient,
            "R2_linear_n1": fit1.r_squared,
            "K_parabolic_n2": fit2.coefficient,
            "R2_parabolic_n2": fit2.r_squared,
            "free_n": free.exponent,
            "R2_free_n": free.r_squared,
            "series_KT": kt_series,
            "series_KG": kg_series,
            "R2_GT_series": series_r2,
            "preferred_simple_scaling": (
                "linear" if fit1.r_squared > fit2.r_squared else "parabolic"
            ),
            "jerkiness_CV": standard.get("jerkiness_CV", np.nan),
            "stationary_fraction": standard.get("stationary_fraction", np.nan),
            "motion_top_1pct": standard.get("motion_top_1pct", np.nan),
            "motion_top_5pct": standard.get("motion_top_5pct", np.nan),
            "motion_top_10pct": standard.get("motion_top_10pct", np.nan),
            "Fano": standard.get("Fano", np.nan),
            "burstiness": standard.get("burstiness", np.nan),
            "reverse_motion_fraction": standard.get("reverse_motion_fraction", np.nan),
            "pinned_fraction": standard.get("pinned_fraction", np.nan),
            "number_of_events": standard.get("number_of_events", 0),
            **_frame_shear(run),
            **_activation_work(run),
            **_pin_durations(run),
            **_shear_release_audit(run),
            "run": str(run),
        })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows).sort_values("regime")
    table.to_csv(output, index=False)
    md = output.with_suffix(".md")
    lines = [
        "# Jerky mechanism integrity analysis",
        "",
        "Primary questions: G should retain approximately parabolic mean scaling while becoming intermittent; T should approach a linear mean kinetic limit; GT is also fit to t=(R-R0)/KT+(R^2-R0^2)/KG. Shear and climb work are reported separately from capillary work for event-resolved GB barrier crossings.",
        "",
        "```text",
        table.to_string(index=False),
        "```",
        "",
    ]
    md.write_text("\n".join(lines), encoding="utf-8")
    print(table.to_string(index=False))
    print(md)


if __name__ == "__main__":
    main()
