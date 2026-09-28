#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, wasserstein_distance

from grain_growth_pf.analysis.long_time import (
    differential_effective_exponent,
    grain_size_moments,
    profile_growth_window,
)
from grain_growth_pf.io.event_ledger import event_ledger_path, iter_event_ledger


PROGRESS_WINDOWS = (
    (1.00, 1.20), (1.15, 1.40), (1.30, 1.60), (1.50, 1.85),
    (1.75, 2.15), (2.00, 2.50), (2.35, 2.90), (2.70, 3.20),
)
STATE_PROGRESS = (1.00, 1.20, 1.40, 1.60, 1.85, 2.15, 2.50, 2.90, 3.20)
SIZE_COLUMNS = ("G_mean", "G_population", "G_area_weighted")


def _runs(root: Path) -> list[tuple[str, Path, dict]]:
    campaign = json.loads((root / "campaign_manifest.json").read_text())
    result = []
    for raw in campaign["runs"]:
        path = Path(raw)
        manifest = json.loads((path / "manifest.json").read_text())
        result.append((manifest["config"]["regime"], path, manifest))
    return result


def _trajectory(run: Path, regime: str) -> pd.DataFrame:
    tracks = pd.read_csv(
        run / "grain_tracks.csv",
        usecols=["time", "step", "grain_id", "area", "radius", "neighbors"],
    )
    rows = []
    for (step, time), group in tracks.groupby(["step", "time"], sort=True):
        moments = grain_size_moments(group["area"].to_numpy(float))
        rows.append({
            "regime": regime, "step": int(step), "time": float(time),
            "grain_count": int(group["grain_id"].nunique()),
            "G_mean": moments.number_mean, "G_population": moments.population,
            "G_area_weighted": moments.area_weighted,
            "mean_side_number": float(group["neighbors"].mean()),
            "total_area": float(group["area"].sum()),
        })
    result = pd.DataFrame(rows).sort_values("step").drop_duplicates("step")
    for name in SIZE_COLUMNS:
        result[f"{name}_over_G0"] = result[name] / result[name].iloc[0]
    rate, differential_n = differential_effective_exponent(
        result["time"].to_numpy(), result["G_population"].to_numpy(),
        window_length=min(21, len(result) if len(result) % 2 else len(result) - 1),
    ) if len(result) >= 5 else (np.full(len(result), np.nan), np.full(len(result), np.nan))
    result["Gdot_smoothed"] = rate
    result["n_eff_differential"] = differential_n
    return result


def _window_fits(trajectory: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    fits, profiles = [], []
    progress = trajectory["G_population_over_G0"]
    windows: list[tuple[str, float, float]] = [
        (f"x_{lo:.2f}_{hi:.2f}", lo, hi) for lo, hi in PROGRESS_WINDOWS
    ]
    windows.extend(
        (f"cumulative_to_{hi:.2f}", 1.0, hi) for hi in STATE_PROGRESS[1:]
    )
    for window_name, lo, hi in windows:
        selected = trajectory[(progress >= lo) & (progress <= hi)]
        for measure in SIZE_COLUMNS:
            base = {
                "regime": trajectory["regime"].iloc[0], "window": window_name,
                "window_start_x": lo, "window_end_x": hi, "size_measure": measure,
                "samples": len(selected),
                "start_grains": int(selected["grain_count"].iloc[0]) if len(selected) else np.nan,
                "end_grains": int(selected["grain_count"].iloc[-1]) if len(selected) else np.nan,
                "grain_loss_events": int(selected["grain_count"].iloc[0] - selected["grain_count"].iloc[-1]) if len(selected) else 0,
            }
            sufficient = bool(
                len(selected) >= 5
                and selected[measure].iloc[-1] / selected[measure].iloc[0] >= 1.04
            )
            if not sufficient:
                fits.append({**base, "available": False})
                continue
            fit = profile_growth_window(
                selected["time"].to_numpy(), selected[measure].to_numpy()
            )
            fits.append({
                **base, "available": True, "start_time": selected["time"].iloc[0],
                "end_time": selected["time"].iloc[-1], "start_x": progress.loc[selected.index[0]],
                "end_x": progress.loc[selected.index[-1]], "n_best": fit.n_best,
                "K_best": fit.coefficient, "normalized_rmse": fit.normalized_rmse,
                "n_profile_width_5pct": fit.profile_width,
                "n_eff_at_or_above_search_bound": fit.hit_search_bound,
                "K2": fit.k2, "K3": fit.k3,
                "mean_Gdot": fit.mean_growth_rate,
                "normalized_Gdot": fit.normalized_growth_rate,
            })
            for exponent, residual in zip(fit.exponents, fit.residual_profile):
                profiles.append({
                    "regime": base["regime"], "window": window_name,
                    "size_measure": measure, "n": exponent,
                    "normalized_rmse": residual,
                })
    return fits, profiles


def _nearest_step(trajectory: pd.DataFrame, progress: float) -> int | None:
    if trajectory["G_population_over_G0"].max() + 1e-12 < progress:
        return None
    index = (trajectory["G_population_over_G0"] - progress).abs().idxmin()
    return int(trajectory.loc[index, "step"])


def _distribution_rows(run: Path, trajectory: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    tracks = pd.read_csv(run / "grain_tracks.csv")
    boundaries = pd.read_csv(run / "boundary_tracks.csv")
    distributions, collapse = [], []
    previous: dict[str, np.ndarray] = {}
    regime = trajectory["regime"].iloc[0]
    for target in STATE_PROGRESS:
        step = _nearest_step(trajectory, target)
        if step is None:
            continue
        grains = tracks[tracks["step"] == step]
        edges = boundaries[boundaries["step"] == step]
        diameter = 2.0 * np.sqrt(grains["area"].to_numpy(float) / np.pi)
        values = {
            "equivalent_diameter": diameter / diameter.mean(),
            "grain_area": grains["area"].to_numpy(float) / grains["area"].mean(),
            "side_number": grains["neighbors"].to_numpy(float),
            "area_weighted_diameter": np.repeat(
                diameter / np.average(diameter, weights=grains["area"]),
                np.maximum(1, np.rint(grains["area"] / grains["area"].min()).astype(int)),
            ),
            "boundary_length": (
                edges["length"].to_numpy(float) / edges["length"].mean()
                if len(edges) else np.asarray([], dtype=float)
            ),
        }
        actual = float(trajectory.loc[trajectory["step"] == step, "G_population_over_G0"].iloc[0])
        for metric, array in values.items():
            finite = array[np.isfinite(array)]
            for value in finite:
                distributions.append({
                    "regime": regime, "target_x": target, "actual_x": actual,
                    "step": step, "metric": metric, "normalized_value": value,
                })
            if metric in previous and len(finite) and len(previous[metric]):
                collapse.append({
                    "regime": regime, "target_x": target, "actual_x": actual,
                    "metric": metric,
                    "ks_distance": ks_2samp(previous[metric], finite).statistic,
                    "wasserstein_distance": wasserstein_distance(previous[metric], finite),
                })
            previous[metric] = finite
    return distributions, collapse


def _mechanism_residence(run: Path, trajectory: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    path = run / "mechanism_state.csv"
    if not path.exists() or path.stat().st_size == 0:
        return [], []
    aggregates = []
    try:
        chunks = pd.read_csv(path, chunksize=250_000)
        for chunk in chunks:
            if not len(chunk):
                continue
            chunk["unblocked"] = (
                (chunk["G_pending"] == 0) & (chunk["T_pending"] == 0)
                & (chunk["C_pending"] == 0) & (chunk["blocked"] == 0)
            ).astype(int)
            chunk["joint_code"] = (
                chunk["G_pending"].astype(int)
                + 2 * chunk["T_pending"].astype(int)
                + 4 * chunk["C_pending"].astype(int)
            )
            base_columns = [
                "G_pending", "T_pending", "C_pending", "blocked", "unblocked",
                "GB_sink_activity", "TJ_sink_activity",
            ]
            grouped = chunk.groupby("step")[base_columns].agg(["sum", "count"])
            grouped.columns = [f"{a}_{b}" for a, b in grouped.columns]
            aggregates.append(grouped.reset_index())
            stage = pd.crosstab(chunk["step"], chunk["climb_stage"]).reset_index()
            aggregates.append(stage.add_prefix("stage_").rename(columns={"stage_step": "step"}))
    except pd.errors.EmptyDataError:
        return [], []
    if not aggregates:
        return [], []
    merged = pd.concat(aggregates, ignore_index=True).groupby("step", as_index=False).sum(numeric_only=True)
    progress_by_step = trajectory.set_index("step")["G_population_over_G0"]
    merged["x"] = merged["step"].map(progress_by_step)
    residence, joint = [], []
    regime = trajectory["regime"].iloc[0]
    for lo, hi in PROGRESS_WINDOWS:
        selected = merged[(merged["x"] >= lo) & (merged["x"] <= hi)]
        if not len(selected):
            continue
        total = float(selected.filter(regex="_count$").max(axis=1).sum())
        def fraction(column: str) -> float:
            return float(selected.get(column, pd.Series(dtype=float)).sum() / total) if total else 0.0
        residence.append({
            "regime": regime, "window_start_x": lo, "window_end_x": hi,
            "G_barrier_waiting": fraction("G_pending_sum"),
            "T_compatibility_waiting": fraction("T_pending_sum"),
            "C_pending": fraction("C_pending_sum"),
            "GB_sink_waiting": fraction("GB_sink_activity_sum"),
            "TJ_sink_waiting": fraction("TJ_sink_activity_sum"),
            "unblocked_capillary_migration": fraction("unblocked_sum"),
            "climb_nucleation": fraction("stage_nucleation"),
            "climb_exchange": fraction("stage_exchange"),
            "climb_transport": fraction("stage_transport"),
            "entity_samples": total,
        })
    # Exact joint occupancy is derived directly in a second bounded-memory pass.
    joint_counts: dict[tuple[int, int], int] = {}
    for chunk in pd.read_csv(path, usecols=["step", "G_pending", "T_pending", "C_pending"], chunksize=250_000):
        codes = chunk["G_pending"].astype(int) + 2 * chunk["T_pending"].astype(int) + 4 * chunk["C_pending"].astype(int)
        counts = pd.DataFrame({"step": chunk["step"], "code": codes}).value_counts()
        for (step, code), count in counts.items():
            joint_counts[(int(step), int(code))] = joint_counts.get((int(step), int(code)), 0) + int(count)
    for lo, hi in PROGRESS_WINDOWS:
        steps = set(trajectory.loc[(trajectory["G_population_over_G0"] >= lo) & (trajectory["G_population_over_G0"] <= hi), "step"].astype(int))
        selected = {(step, code): count for (step, code), count in joint_counts.items() if step in steps}
        total = sum(selected.values())
        for code in range(8):
            joint.append({
                "regime": regime, "window_start_x": lo, "window_end_x": hi,
                "joint_code": code, "G_pending": bool(code & 1),
                "T_pending": bool(code & 2), "C_pending": bool(code & 4),
                "fraction": sum(count for (step, item), count in selected.items() if item == code) / total if total else np.nan,
            })
    return residence, joint


def _frame_rows(run: Path, trajectory: pd.DataFrame) -> list[dict]:
    rows = []
    for path in sorted((run / "frames").glob("frame-*.npz")):
        with np.load(path) as data:
            row = {"regime": trajectory["regime"].iloc[0], "frame_path": str(path)}
            for name in (
                "step", "time", "grain_count", "G_mean", "G_population",
                "G_area_weighted", "G_over_G0", "G_occupancy", "T_occupancy",
                "C_occupancy", "GB_sink_fraction", "TJ_sink_fraction",
                "mean_abs_tau_int_active_length", "mean_abs_p_shear_active_length",
                "mean_chi_s_active_length", "chi_s_p50_active_domain",
                "chi_s_p90_active_domain", "chi_s_p95_active_domain",
                "fraction_active_gb_length_chi_s_near_one",
                "fraction_active_gb_length_chi_s_gt_0p8",
                "fraction_active_gb_length_chi_s_gt_one",
                "stored_shear_energy_per_active_gb_length",
                "active_shear_bearing_gb_fraction", "conservation_residual",
            ):
                row[name] = float(data[name]) if name in data else np.nan
            rows.append(row)
    return rows


def _event_window_rows(run: Path, trajectory: pd.DataFrame, boundary: pd.DataFrame) -> list[dict]:
    counts: dict[int, int] = {}
    type_counts: dict[tuple[int, str], int] = {}
    path = event_ledger_path(run)
    if path.exists():
        for chunk in iter_event_ledger(path, columns=["step", "event_type"]):
            for (step, event_type), count in chunk.groupby(["step", "event_type"]).size().items():
                key = int(step)
                counts[key] = counts.get(key, 0) + int(count)
                type_counts[(key, str(event_type))] = type_counts.get((key, str(event_type)), 0) + int(count)
    rows = []
    regime = trajectory["regime"].iloc[0]
    for lo, hi in PROGRESS_WINDOWS:
        selected = trajectory[(trajectory["G_population_over_G0"] >= lo) & (trajectory["G_population_over_G0"] <= hi)]
        if len(selected) < 2:
            continue
        steps = selected["step"].astype(int).tolist()
        event_values = np.asarray([counts.get(step, 0) for step in steps], dtype=float)
        increments = np.maximum(np.diff(selected["G_population"].to_numpy(float)), 0.0)
        ordered = np.sort(increments)[::-1]
        total_motion = ordered.sum()
        gb = boundary[boundary["step"].isin(steps)].groupby("step")["length"].sum().mean()
        area_loss = max(float(selected["total_area"].iloc[0] - selected["total_area"].iloc[-1]), 0.0)
        rate = event_values.sum() / max(float(selected["time"].iloc[-1] - selected["time"].iloc[0]), np.finfo(float).tiny)
        rows.append({
            "regime": regime, "window_start_x": lo, "window_end_x": hi,
            "event_count": int(event_values.sum()), "event_rate": rate,
            "event_rate_per_GB_length": rate / gb if gb and np.isfinite(gb) else np.nan,
            "event_count_per_area_loss": event_values.sum() / area_loss if area_loss else np.nan,
            "event_count_fano": event_values.var() / event_values.mean() if event_values.mean() else 0.0,
            "motion_top_1pct": ordered[:max(1, int(np.ceil(len(ordered) * .01)))].sum() / total_motion if total_motion else 0.0,
            "motion_top_5pct": ordered[:max(1, int(np.ceil(len(ordered) * .05)))].sum() / total_motion if total_motion else 0.0,
            "motion_top_10pct": ordered[:max(1, int(np.ceil(len(ordered) * .10)))].sum() / total_motion if total_motion else 0.0,
            "stationary_fraction": float(np.mean(increments <= 1e-12)) if len(increments) else np.nan,
            "velocity_autocorrelation": float(np.corrcoef(increments[:-1], increments[1:])[0, 1]) if len(increments) > 2 and np.std(increments) else 0.0,
        })
    return rows


def _activation_work_rows(run: Path, trajectory: pd.DataFrame) -> list[dict]:
    path = run / "activation_work.csv"
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        work = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return []
    if not len(work):
        return []
    progress = trajectory.set_index("step")["G_population_over_G0"]
    known_steps = progress.index.to_numpy(int)
    known_x = progress.to_numpy(float)
    event_steps = pd.to_numeric(work["step"], errors="coerce").fillna(0).to_numpy(int)
    work["x"] = np.interp(event_steps, known_steps, known_x)
    rows = []
    columns = (
        "work_capillary", "work_shear", "work_free_volume",
        "work_total_without_tj_residual", "effective_DeltaG",
    )
    for lo, hi in PROGRESS_WINDOWS:
        selected = work[(work["x"] >= lo) & (work["x"] <= hi)]
        if not len(selected):
            continue
        row = {
            "regime": trajectory["regime"].iloc[0],
            "window_start_x": lo, "window_end_x": hi, "events": len(selected),
        }
        for column in columns:
            values = pd.to_numeric(selected.get(column), errors="coerce").dropna()
            row[f"{column}_mean"] = float(values.mean()) if len(values) else np.nan
            row[f"{column}_p95_abs"] = float(np.quantile(np.abs(values), .95)) if len(values) else np.nan
        rows.append(row)
    return rows


def _event_triggered_rows(run: Path, trajectory: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    event_path = run / "event_trace_events.parquet"
    trace_path = run / "event_traces.parquet"
    event_parts = sorted(event_path.glob("part-*.parquet")) if event_path.is_dir() else []
    trace_parts = sorted(trace_path.glob("part-*.parquet")) if trace_path.is_dir() else []
    if not event_parts or not trace_parts:
        return [], []
    events = pd.read_parquet(
        event_path,
        columns=[
            "event_id", "event_step", "event_type", "trace_start_step",
            "trace_end_step",
        ],
    )
    event_records: dict[str, list[tuple[int, int, int, str]]] = defaultdict(list)
    for event in events.itertuples(index=False):
        event_records[str(event.event_id)].append(
            (
                int(event.event_step), int(event.trace_start_step),
                int(event.trace_end_step), str(event.event_type),
            )
        )
    known_steps = trajectory["step"].to_numpy(int)
    known_x = trajectory["G_population_over_G0"].to_numpy(float)
    sums: dict[tuple[float, float, int, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    import pyarrow.parquet as pq
    for part in trace_parts:
        parquet = pq.ParquetFile(part)
        for batch in parquet.iter_batches(
            batch_size=250_000,
            columns=["step", "event_ids", "entity_type", "local_normal_velocity"],
        ):
            frame = batch.to_pandas()
            for row in frame.itertuples(index=False):
                if not isinstance(row.event_ids, str) or not row.event_ids:
                    continue
                velocity = float(row.local_normal_velocity) if pd.notna(row.local_normal_velocity) else np.nan
                if not np.isfinite(velocity):
                    continue
                for event_id in row.event_ids.split(";"):
                    candidates = [
                        item for item in event_records.get(event_id, [])
                        if item[1] <= int(row.step) <= item[2]
                    ]
                    if not candidates:
                        continue
                    # Event IDs are domain-local and can recur after topology
                    # retirement/recreation.  The trace interval disambiguates
                    # them; nearest-event tie breaking is deterministic.
                    event_step, _, _, _ = min(
                        candidates, key=lambda item: (abs(int(row.step) - item[0]), item[0])
                    )
                    event_x = float(np.interp(event_step, known_steps, known_x))
                    window = next(
                        ((lo, hi) for lo, hi in PROGRESS_WINDOWS if lo <= event_x <= hi),
                        None,
                    )
                    if window is None:
                        continue
                    key = (window[0], window[1], int(row.step) - event_step, str(row.entity_type))
                    sums[key][0] += velocity
                    sums[key][1] += abs(velocity)
                    sums[key][2] += 1.0
    regime = trajectory["regime"].iloc[0]
    rows = []
    for (lo, hi, relative_step, entity_type), (signed, magnitude, count) in sorted(sums.items()):
        rows.append({
            "regime": regime, "window_start_x": lo, "window_end_x": hi,
            "relative_step": relative_step, "entity_type": entity_type,
            "samples": int(count), "mean_local_velocity": signed / count,
            "mean_abs_local_velocity": magnitude / count,
        })
    frame = pd.DataFrame(rows)
    summary = []
    if len(frame):
        for (lo, hi, entity_type), group in frame.groupby(
            ["window_start_x", "window_end_x", "entity_type"]
        ):
            pre = group[group["relative_step"].between(-5, -1)]
            post = group[group["relative_step"].between(1, 5)]
            def weighted_mean(selected: pd.DataFrame) -> float:
                return float(np.average(
                    selected["mean_abs_local_velocity"], weights=selected["samples"]
                )) if len(selected) else np.nan
            summary.append({
                "regime": regime, "window_start_x": lo, "window_end_x": hi,
                "entity_type": entity_type, "pre_mean_abs_velocity": weighted_mean(pre),
                "post_mean_abs_velocity": weighted_mean(post),
                "post_minus_pre": weighted_mean(post) - weighted_mean(pre),
                "sampled_events": sum(
                    1 for event in events.itertuples(index=False)
                    for step in (int(event.event_step),)
                    if lo <= float(np.interp(step, known_steps, known_x)) <= hi
                    and pd.notna(event.event_type)
                ),
            })
    return rows, summary


def _velocity_and_wait_rows(
    run: Path, trajectory: pd.DataFrame,
) -> tuple[list[dict], list[dict], list[dict]]:
    tracks = pd.read_csv(run / "grain_tracks.csv", usecols=["grain_id", "step", "time", "radius"])
    tracks = tracks.sort_values(["grain_id", "time"])
    tracks["dt"] = tracks.groupby("grain_id")["time"].diff()
    tracks["radial_velocity"] = tracks.groupby("grain_id")["radius"].diff() / tracks["dt"]
    progress = trajectory.set_index("step")["G_population_over_G0"]
    tracks["x"] = tracks["step"].map(progress)
    summary, histogram, waits = [], [], []
    regime = trajectory["regime"].iloc[0]
    for lo, hi in PROGRESS_WINDOWS:
        velocities = tracks.loc[
            (tracks["x"] >= lo) & (tracks["x"] <= hi), "radial_velocity"
        ].dropna().to_numpy(float)
        if len(velocities):
            absolute = np.abs(velocities)
            summary.append({
                "regime": regime, "window_start_x": lo, "window_end_x": hi,
                "samples": len(velocities), "velocity_mean": float(velocities.mean()),
                "velocity_std": float(velocities.std()),
                "abs_velocity_p50": float(np.quantile(absolute, .50)),
                "abs_velocity_p90": float(np.quantile(absolute, .90)),
                "abs_velocity_p95": float(np.quantile(absolute, .95)),
                "grain_stationary_fraction": float(np.mean(absolute <= 1e-12)),
            })
            counts, edges = np.histogram(velocities, bins=41, density=True)
            for count, left, right in zip(counts, edges[:-1], edges[1:]):
                histogram.append({
                    "regime": regime, "window_start_x": lo, "window_end_x": hi,
                    "bin_left": left, "bin_right": right, "density": count,
                })
        selected = trajectory[
            (trajectory["G_population_over_G0"] >= lo)
            & (trajectory["G_population_over_G0"] <= hi)
        ]
        if len(selected) < 2:
            continue
        dt = np.diff(selected["time"].to_numpy(float))
        velocity = np.diff(selected["G_population"].to_numpy(float)) / dt
        positive = velocity[velocity > 0.0]
        tolerance = max(float(np.median(positive) * .05) if len(positive) else 0.0, 1e-12)
        stationary = velocity <= tolerance
        start = 0
        while start < len(stationary):
            stop = start + 1
            while stop < len(stationary) and stationary[stop] == stationary[start]:
                stop += 1
            waits.append({
                "regime": regime, "window_start_x": lo, "window_end_x": hi,
                "episode_type": "wait" if stationary[start] else "burst",
                "duration": float(dt[start:stop].sum()),
                "size": float(np.maximum(velocity[start:stop], 0.0) @ dt[start:stop]),
                "intervals": stop - start, "stationary_tolerance": tolerance,
            })
            start = stop
    return summary, histogram, waits


def _plot_lines(frame: pd.DataFrame, x: str, ys: Iterable[str], output: Path, ylabel: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for regime, group in frame.groupby("regime"):
        for y in ys:
            if y in group and group[y].notna().any():
                label = regime if len(tuple(ys)) == 1 else f"{regime} {y}"
                ax.plot(group[x], group[y], marker="o", ms=2, label=label)
    ax.set_xlabel(x)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    root = Path(args.campaign)
    output = Path(args.output) if args.output else root / "analysis"
    output.mkdir(parents=True, exist_ok=True)

    trajectories, fits, profiles = [], [], []
    distributions, collapse, residence, joint = [], [], [], []
    frames, jerkiness, activation_work, event_triggered, event_response = [], [], [], [], []
    velocity_summary, velocity_histogram, burst_wait, run_summary = [], [], [], []
    for regime, run, manifest in _runs(root):
        if manifest.get("status") != "completed" and not args.allow_incomplete:
            raise RuntimeError(f"incomplete run: {run}")
        trajectory = _trajectory(run, regime)
        trajectories.append(trajectory)
        local_fits, local_profiles = _window_fits(trajectory)
        fits.extend(local_fits); profiles.extend(local_profiles)
        local_dist, local_collapse = _distribution_rows(run, trajectory)
        distributions.extend(local_dist); collapse.extend(local_collapse)
        local_residence, local_joint = _mechanism_residence(run, trajectory)
        residence.extend(local_residence); joint.extend(local_joint)
        frames.extend(_frame_rows(run, trajectory))
        boundary = pd.read_csv(run / "boundary_tracks.csv", usecols=["step", "length"])
        jerkiness.extend(_event_window_rows(run, trajectory, boundary))
        activation_work.extend(_activation_work_rows(run, trajectory))
        local_triggered, local_response = _event_triggered_rows(run, trajectory)
        event_triggered.extend(local_triggered); event_response.extend(local_response)
        local_velocity, local_histogram, local_wait = _velocity_and_wait_rows(run, trajectory)
        velocity_summary.extend(local_velocity)
        velocity_histogram.extend(local_histogram)
        burst_wait.extend(local_wait)
        conservation = 0.0
        inventory_path = run / "defect_inventory.csv"
        if inventory_path.exists() and inventory_path.stat().st_size:
            inventory = pd.read_csv(
                inventory_path, usecols=["max_abs_conservation_residual"]
            )
            conservation = float(inventory["max_abs_conservation_residual"].max())
        terminal_tracks = pd.read_csv(
            run / "grain_tracks.csv", usecols=["step", "area", "perimeter"]
        )
        terminal_tracks = terminal_tracks[
            terminal_tracks["step"] == terminal_tracks["step"].max()
        ]
        compactness = terminal_tracks["perimeter"] / np.sqrt(
            4.0 * np.pi * terminal_tracks["area"]
        )
        run_summary.append({
            "regime": regime, "path": str(run), "source_sha": manifest.get("git_sha"),
            "status": manifest.get("status"), "start_N": int(trajectory["grain_count"].iloc[0]),
            "end_N": int(trajectory["grain_count"].iloc[-1]),
            "start_G": trajectory["G_population"].iloc[0],
            "end_G": trajectory["G_population"].iloc[-1],
            "G_over_G0": trajectory["G_population_over_G0"].iloc[-1],
            "steps": int(trajectory["step"].iloc[-1]), "time": trajectory["time"].iloc[-1],
            "target_N100_reached": int(trajectory["grain_count"].iloc[-1]) <= 100,
            "max_abs_conservation_residual": conservation,
            "terminal_mean_compactness": float(compactness.mean()),
            "terminal_max_compactness": float(compactness.max()),
            "movie_path": str(run / "microstructure.mp4"),
        })

    trajectory_frame = pd.concat(trajectories, ignore_index=True)
    fit_frame = pd.DataFrame(fits)
    profile_frame = pd.DataFrame(profiles)
    frame_summary = pd.DataFrame(frames)
    residence_frame = pd.DataFrame(residence)
    if len(residence_frame) and len(frame_summary):
        shear_residence = []
        for row in residence_frame.itertuples(index=False):
            selected = frame_summary[
                (frame_summary["regime"] == row.regime)
                & (frame_summary["G_over_G0"] >= row.window_start_x)
                & (frame_summary["G_over_G0"] <= row.window_end_x)
            ]
            shear_residence.append(
                float(selected["fraction_active_gb_length_chi_s_gt_0p8"].mean())
                if len(selected) else np.nan
            )
        residence_frame["shear_backstress_limited_residence"] = shear_residence
    outputs = {
        "grain_size_trajectory.csv": trajectory_frame,
        "growth_law_windows.csv": fit_frame,
        "growth_law_profiles.csv": profile_frame,
        "normalized_distributions.csv": pd.DataFrame(distributions),
        "distribution_collapse.csv": pd.DataFrame(collapse),
        "mechanism_residence_vs_grain_size.csv": residence_frame,
        "mechanism_joint_occupancy_vs_grain_size.csv": pd.DataFrame(joint),
        "frame_mechanisms_vs_grain_size.csv": frame_summary,
        "jerkiness_vs_grain_size.csv": pd.DataFrame(jerkiness),
        "activation_work_vs_grain_size.csv": pd.DataFrame(activation_work),
        "event_triggered_velocity_vs_grain_size.csv": pd.DataFrame(event_triggered),
        "event_triggered_response_summary.csv": pd.DataFrame(event_response),
        "grain_velocity_summary_vs_grain_size.csv": pd.DataFrame(velocity_summary),
        "grain_velocity_distribution_vs_grain_size.csv": pd.DataFrame(velocity_histogram),
        "burst_wait_distribution_vs_grain_size.csv": pd.DataFrame(burst_wait),
        "long_time_run_summary.csv": pd.DataFrame(run_summary),
    }
    for name, table in outputs.items():
        table.to_csv(output / name, index=False)

    _plot_lines(trajectory_frame, "time", ["G_population"], output / "G_vs_time.png", "G")
    plot = trajectory_frame.copy()
    plot["G2_minus_G02"] = plot.groupby("regime")["G_population"].transform(lambda x: x**2 - x.iloc[0]**2)
    plot["G3_minus_G03"] = plot.groupby("regime")["G_population"].transform(lambda x: x**3 - x.iloc[0]**3)
    _plot_lines(plot, "time", ["G2_minus_G02"], output / "G2_minus_G02_vs_time.png", "G²-G₀²")
    _plot_lines(plot, "time", ["G3_minus_G03"], output / "G3_minus_G03_vs_time.png", "G³-G₀³")
    _plot_lines(trajectory_frame, "G_population", ["Gdot_smoothed"], output / "Gdot_vs_G.png", "dG/dt")
    available = fit_frame[(fit_frame.get("available", False) == True) & fit_frame["window"].str.startswith("x_")]
    _plot_lines(available, "end_x", ["n_best"], output / "n_eff_vs_G_over_G0.png", "n_eff")
    _plot_lines(available, "end_x", ["K2", "K3"], output / "K2_K3_vs_G_over_G0.png", "fixed-law coefficient")
    _plot_lines(available, "end_x", ["K_best"], output / "K_best_vs_G_over_G0.png", "K_n")
    if len(profile_frame):
        representative = profile_frame[profile_frame["size_measure"] == "G_population"]
        _plot_lines(representative, "n", ["normalized_rmse"], output / "profile_residual_vs_n.png", "normalized RMSE")
    if len(frame_summary):
        _plot_lines(frame_summary, "G_over_G0", ["G_occupancy", "T_occupancy", "C_occupancy"], output / "GTC_occupancy_vs_G_over_G0.png", "occupancy")
        _plot_lines(frame_summary, "G_over_G0", ["GB_sink_fraction", "TJ_sink_fraction"], output / "sink_partition_vs_G_over_G0.png", "sink fraction")
        _plot_lines(frame_summary, "G_over_G0", ["mean_chi_s_active_length", "mean_abs_tau_int_active_length", "stored_shear_energy_per_active_gb_length"], output / "shear_vs_G_over_G0.png", "shear metric")
    if len(outputs["distribution_collapse.csv"]):
        _plot_lines(outputs["distribution_collapse.csv"], "actual_x", ["ks_distance", "wasserstein_distance"], output / "distribution_collapse_vs_G_over_G0.png", "distance")
    distribution_frame = outputs["normalized_distributions.csv"]
    if len(distribution_frame):
        selected_distribution = distribution_frame[
            distribution_frame["metric"] == "equivalent_diameter"
        ]
        regimes = list(selected_distribution["regime"].unique())
        fig, axes = plt.subplots(
            max(1, int(np.ceil(len(regimes) / 3))), 3,
            figsize=(12, 3.5 * max(1, int(np.ceil(len(regimes) / 3)))),
            squeeze=False, constrained_layout=True,
        )
        for ax, regime in zip(axes.flat, regimes):
            group = selected_distribution[selected_distribution["regime"] == regime]
            for target in sorted(group["target_x"].unique()):
                values = group.loc[group["target_x"] == target, "normalized_value"]
                density, edges = np.histogram(values, bins=np.linspace(0, 3, 41), density=True)
                ax.plot(0.5 * (edges[:-1] + edges[1:]), density, label=f"x={target:g}")
            ax.set_title(regime); ax.set_xlabel("G_i/<G>"); ax.set_ylabel("density")
        for ax in axes.flat[len(regimes):]:
            ax.axis("off")
        if regimes:
            axes.flat[0].legend(fontsize=7)
        fig.savefig(output / "normalized_grain_size_distributions.png", dpi=180)
        plt.close(fig)
    if len(outputs["jerkiness_vs_grain_size.csv"]):
        _plot_lines(outputs["jerkiness_vs_grain_size.csv"], "window_end_x", ["motion_top_1pct", "motion_top_5pct", "motion_top_10pct", "event_count_fano", "stationary_fraction"], output / "jerkiness_vs_G_over_G0.png", "jerkiness metric")
    if len(outputs["mechanism_residence_vs_grain_size.csv"]):
        _plot_lines(outputs["mechanism_residence_vs_grain_size.csv"], "window_end_x", ["G_barrier_waiting", "T_compatibility_waiting", "C_pending", "GB_sink_waiting", "TJ_sink_waiting", "unblocked_capillary_migration"], output / "mechanism_residence_vs_G_over_G0.png", "residence fraction")
    if len(outputs["event_triggered_velocity_vs_grain_size.csv"]):
        _plot_lines(outputs["event_triggered_velocity_vs_grain_size.csv"], "relative_step", ["mean_abs_local_velocity"], output / "event_triggered_velocity_vs_relative_step.png", "event-linked |v|")
    if len(outputs["activation_work_vs_grain_size.csv"]):
        _plot_lines(outputs["activation_work_vs_grain_size.csv"], "window_end_x", ["work_capillary_p95_abs", "work_shear_p95_abs", "work_free_volume_p95_abs"], output / "activation_work_vs_G_over_G0.png", "p95 |activation work|")
    print(output)


if __name__ == "__main__":
    main()
