#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.dataset as pds

from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger

from analyze_shear_screen import (
    RELEASE_TYPES,
    TOPOLOGY_WINDOWS,
    _event_responses,
    _growth,
    _population_at_steps,
    _run_summary,
    _trace_table,
    _runs,
    _topology_label,
    _window_metrics,
)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return float("nan")
    values, weights = values[valid], weights[valid]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    centers = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
    return float(np.interp(float(q), centers, values))


def arrest_episodes(
    steps: np.ndarray, near_arrest: np.ndarray, curvature: np.ndarray, dt: float,
) -> list[dict[str, float]]:
    steps = np.asarray(steps, dtype=int)
    near_arrest = np.asarray(near_arrest, dtype=bool)
    curvature = np.asarray(curvature, dtype=float)
    episodes: list[dict[str, float]] = []
    start: int | None = None
    for index in range(len(steps)):
        continuous = index > 0 and steps[index] == steps[index - 1] + 1
        if near_arrest[index] and (start is None or continuous):
            if start is None:
                start = index
        elif near_arrest[index]:
            if start is not None:
                episodes.append(_episode_row(steps, curvature, start, index - 1, dt))
            start = index
        elif start is not None:
            episodes.append(_episode_row(steps, curvature, start, index - 1, dt))
            start = None
    if start is not None:
        episodes.append(_episode_row(steps, curvature, start, len(steps) - 1, dt))
    return episodes


def _episode_row(
    steps: np.ndarray, curvature: np.ndarray, start: int, end: int, dt: float,
) -> dict[str, float]:
    count = int(steps[end] - steps[start] + 1)
    return {
        "start_step": int(steps[start]), "end_step": int(steps[end]),
        "duration_steps": count, "duration_time": float(count * dt),
        "entry_curvature": float(curvature[start]),
        "exit_curvature": float(curvature[end]),
    }


def segmented_change_point(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    order = np.argsort(x)
    x, y = x[order], y[order]
    if len(x) < 4:
        return {"critical_Ks": float("nan"), "sse": float("nan")}
    best: dict[str, float] | None = None
    for critical in x[1:-1]:
        design = np.column_stack((np.ones(len(x)), x, np.maximum(0.0, x - critical)))
        coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ coefficients
        candidate = {
            "critical_Ks": float(critical), "sse": float(residual @ residual),
            "intercept": float(coefficients[0]),
            "slope_below": float(coefficients[1]),
            "slope_change_above": float(coefficients[2]),
        }
        if best is None or candidate["sse"] < best["sse"]:
            best = candidate
    assert best is not None
    return best


def _family(regime: str) -> str:
    return re.sub(r"_Ks\d+$", "", regime)


def _read_trace(run: Path) -> pd.DataFrame:
    parquet = run / "event_traces.parquet"
    csv_path = run / "event_traces.csv"
    requested = [
        "entity_type", "entity_id", "step", "time", "local_normal_velocity",
        "curvature", "gb_length", "shear_state_s", "tau_int", "p_cap",
        "p_chem", "p_shear", "p_event", "p_net", "p_applied_total", "chi_s",
        "local_shear_energy",
    ]
    if parquet.exists():
        available = set(pds.dataset(parquet, format="parquet").schema.names)
        frame = pd.read_parquet(
            parquet, columns=[column for column in requested if column in available]
        )
    elif csv_path.exists():
        available = set(pd.read_csv(csv_path, nrows=0).columns)
        frame = pd.read_csv(
            csv_path, usecols=[column for column in requested if column in available]
        )
    else:
        return pd.DataFrame()
    return frame[frame.get("entity_type", "GB").astype(str).eq("GB")].copy()


def _add_force_columns(
    trace: pd.DataFrame, *, beta: float, stiffness: float, gb_energy: float,
) -> pd.DataFrame:
    trace = trace.copy()
    numeric = (
        "step", "time", "curvature", "gb_length", "shear_state_s", "tau_int",
        "p_cap", "p_chem", "p_shear", "p_event", "p_net", "chi_s",
        "p_applied_total", "local_shear_energy",
    )
    for column in numeric:
        if column in trace:
            trace[column] = pd.to_numeric(trace[column], errors="coerce")
    if "p_cap" not in trace:
        trace["p_cap"] = gb_energy * trace["curvature"]
    if "p_chem" not in trace:
        trace["p_chem"] = 0.0
    if "p_shear" not in trace:
        trace["p_shear"] = beta * trace["tau_int"]
    if "p_event" not in trace:
        trace["p_event"] = 0.0
    if "p_net" not in trace:
        trace["p_net"] = trace["p_cap"] + trace["p_chem"] + trace["p_shear"]
    if "p_applied_total" not in trace:
        trace["p_applied_total"] = trace["p_net"] + trace["p_event"]
    if "chi_s" not in trace:
        denominator = trace["p_cap"] + trace["p_chem"]
        trace["chi_s"] = np.where(
            denominator.abs() > 1e-14, -trace["p_shear"] / denominator, np.nan
        )
    if "local_shear_energy" not in trace:
        trace["local_shear_energy"] = 0.5 * stiffness * trace["shear_state_s"] ** 2
    return trace


def _force_windows(
    regime: str, seed: int, stiffness: float, growth: pd.DataFrame,
    trace: pd.DataFrame, dt: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    episodes: list[dict[str, Any]] = []
    if trace.empty:
        return rows, episodes
    trace["topology_window"] = [
        _topology_label(value) for value in _population_at_steps(
            growth, trace["step"].to_numpy(int)
        )
    ]
    trace["active_shear"] = trace["shear_state_s"].abs() > 1e-12
    for high, low in TOPOLOGY_WINDOWS[:3]:
        label = f"N{high}_to_{low}"
        frame = trace[(trace["topology_window"] == label) & trace["active_shear"]].copy()
        valid = frame.dropna(subset=["gb_length", "chi_s", "p_net"])
        if frame.empty or valid.empty:
            rows.append({
                "regime": regime, "family": _family(regime), "seed": seed,
                "shear_stiffness": stiffness, "topology_window": label,
                "available": False,
            })
            continue
        weights = valid["gb_length"].to_numpy(float)
        chi = valid["chi_s"].to_numpy(float)
        tau = valid["tau_int"].abs().to_numpy(float)
        pshear = valid["p_shear"].abs().to_numpy(float)
        energy = valid["local_shear_energy"].to_numpy(float)
        total_weight = float(weights.sum())

        sign_changes = 0
        sign_pairs = 0
        for _, entity in valid.sort_values(["entity_id", "step"]).groupby("entity_id"):
            steps = entity["step"].to_numpy(int)
            signs = np.sign(entity["p_net"].to_numpy(float))
            consecutive = np.diff(steps) == 1
            sign_pairs += int(consecutive.sum())
            sign_changes += int(np.sum(consecutive & (signs[1:] * signs[:-1] < 0.0)))
        rows.append({
            "regime": regime, "family": _family(regime), "seed": seed,
            "shear_stiffness": stiffness, "topology_window": label,
            "available": True, "trace_rows": int(len(frame)),
            "mean_abs_tau_domain": float(np.mean(tau)),
            "mean_abs_tau_length": float(np.average(tau, weights=weights)),
            "mean_abs_p_shear_domain": float(np.mean(pshear)),
            "mean_abs_p_shear_length": float(np.average(pshear, weights=weights)),
            "mean_chi_domain": float(np.mean(chi)),
            "mean_chi_length": float(np.average(chi, weights=weights)),
            "chi_p50_length": weighted_quantile(chi, weights, 0.50),
            "chi_p90_length": weighted_quantile(chi, weights, 0.90),
            "chi_p95_length": weighted_quantile(chi, weights, 0.95),
            "chi_p99_length": weighted_quantile(chi, weights, 0.99),
            "chi_p50_domain": float(np.quantile(chi, 0.50)),
            "chi_p90_domain": float(np.quantile(chi, 0.90)),
            "chi_p95_domain": float(np.quantile(chi, 0.95)),
            "chi_p99_domain": float(np.quantile(chi, 0.99)),
            "fraction_length_chi_gt_0p8": float(weights[chi > 0.8].sum() / total_weight),
            "fraction_length_chi_near_one": float(
                weights[(chi > 0.8) & (chi < 1.2)].sum() / total_weight
            ),
            "fraction_length_chi_gt_one": float(weights[chi > 1.0].sum() / total_weight),
            "fraction_domain_chi_gt_0p8": float(np.mean(chi > 0.8)),
            "fraction_domain_chi_near_one": float(np.mean((chi > 0.8) & (chi < 1.2))),
            "fraction_domain_chi_gt_one": float(np.mean(chi > 1.0)),
            "p_net_p05_length": weighted_quantile(valid["p_net"], weights, 0.05),
            "p_net_p50_length": weighted_quantile(valid["p_net"], weights, 0.50),
            "p_net_p95_length": weighted_quantile(valid["p_net"], weights, 0.95),
            "p_net_p05_domain": float(np.quantile(valid["p_net"], 0.05)),
            "p_net_p50_domain": float(np.quantile(valid["p_net"], 0.50)),
            "p_net_p95_domain": float(np.quantile(valid["p_net"], 0.95)),
            "p_net_sign_change_fraction": sign_changes / sign_pairs if sign_pairs else np.nan,
            "energy_per_active_domain": float(np.mean(energy)),
            "energy_per_active_length": float(np.average(energy, weights=weights)),
            "local_energy_p50": float(np.quantile(energy, 0.50)),
            "local_energy_p90": float(np.quantile(energy, 0.90)),
            "local_energy_p95": float(np.quantile(energy, 0.95)),
            "local_energy_p99": float(np.quantile(energy, 0.99)),
        })

    for entity_id, entity in trace.sort_values(["entity_id", "step"]).groupby("entity_id"):
        near = entity["chi_s"].between(0.8, 1.2).to_numpy(bool)
        for episode in arrest_episodes(
            entity["step"].to_numpy(int), near,
            entity["curvature"].to_numpy(float), dt,
        ):
            episodes.append({
                "regime": regime, "family": _family(regime), "seed": seed,
                "shear_stiffness": stiffness, "entity_id": entity_id, **episode,
            })
    return rows, episodes


def _frame_windows(
    regime: str, seed: int, stiffness: float, run: Path, dx: float,
) -> list[dict[str, Any]]:
    """Aggregate unbiased ordinary-frame scalars in matched topology windows.

    Legacy frames predate the normalized scalar fields.  Their length-based
    quantities are reconstructed exactly from the stored boundary/shear maps;
    domain-count quantities remain unavailable rather than being inferred from
    connected pixels.
    """

    detail: list[dict[str, Any]] = []
    for path in sorted((run / "frames").glob("frame-*.npz")):
        with np.load(path) as frame:
            grain_count = int(frame["grain_count"])
            label = _topology_label(float(grain_count))
            boundary = frame["boundary_mask"].astype(bool)
            shear = frame["shear"].astype(float)
            active = boundary & (np.abs(shear) > 1e-12)
            tau = (
                np.abs(frame["shear_stress"].astype(float))
                if "shear_stress" in frame else np.zeros_like(shear)
            )
            p_shear = (
                np.abs(frame["p_shear"].astype(float))
                if "p_shear" in frame else np.zeros_like(shear)
            )
            chi = (
                frame["chi_s"].astype(float)
                if "chi_s" in frame else np.full_like(shear, np.nan)
            )
            p_net = (
                frame["p_net"].astype(float)
                if "p_net" in frame else np.full_like(shear, np.nan)
            )
            valid_chi = active & np.isfinite(chi)
            valid_p_net = active & np.isfinite(p_net)
            active_length = float(frame["active_shear_length"]) if (
                "active_shear_length" in frame
            ) else float(active.sum() * dx)
            total_length = float(frame["total_gb_length"]) if (
                "total_gb_length" in frame
            ) else float(boundary.sum() * dx)
            total_energy = float(frame["stored_shear_energy"])
            energy_per_length = float(frame["stored_shear_energy_per_active_gb_length"]) if (
                "stored_shear_energy_per_active_gb_length" in frame
            ) else (total_energy / active_length if active_length else 0.0)
            local_energy = 0.5 * stiffness * shear[active] ** 2

            def scalar_or(name: str, fallback: float) -> float:
                return float(frame[name]) if name in frame else fallback

            def pixel_quantile(values: np.ndarray, mask: np.ndarray, q: float) -> float:
                selected = values[mask]
                return float(np.quantile(selected, q)) if selected.size else np.nan

            detail.append({
                "topology_window": label,
                "active_shear_bearing_gb_fraction": (
                    float(frame["active_shear_bearing_gb_fraction"])
                    if "active_shear_bearing_gb_fraction" in frame
                    else (active_length / total_length if total_length else 0.0)
                ),
                "stored_shear_energy_per_active_gb_length": energy_per_length,
                "stored_shear_energy_per_active_domain": (
                    float(frame["stored_shear_energy_per_active_domain"])
                    if "stored_shear_energy_per_active_domain" in frame else np.nan
                ),
                "active_shear_domain_count": (
                    float(frame["active_shear_domain_count"])
                    if "active_shear_domain_count" in frame else np.nan
                ),
                "local_shear_energy_p50": (
                    float(frame["local_shear_energy_p50"])
                    if "local_shear_energy_p50" in frame
                    else (float(np.quantile(local_energy, 0.50)) if local_energy.size else 0.0)
                ),
                "local_shear_energy_p90": (
                    float(frame["local_shear_energy_p90"])
                    if "local_shear_energy_p90" in frame
                    else (float(np.quantile(local_energy, 0.90)) if local_energy.size else 0.0)
                ),
                "local_shear_energy_p95": (
                    float(frame["local_shear_energy_p95"])
                    if "local_shear_energy_p95" in frame
                    else (float(np.quantile(local_energy, 0.95)) if local_energy.size else 0.0)
                ),
                "local_shear_energy_p99": (
                    float(frame["local_shear_energy_p99"])
                    if "local_shear_energy_p99" in frame
                    else (float(np.quantile(local_energy, 0.99)) if local_energy.size else 0.0)
                ),
                "mean_abs_tau_int_active_domain": scalar_or(
                    "mean_abs_tau_int_active_domain",
                    float(tau[active].mean()) if np.any(active) else np.nan,
                ),
                "mean_abs_tau_int_active_length": scalar_or(
                    "mean_abs_tau_int_active_length",
                    float(tau[active].mean()) if np.any(active) else np.nan,
                ),
                "mean_abs_p_shear_active_domain": scalar_or(
                    "mean_abs_p_shear_active_domain",
                    float(p_shear[active].mean()) if np.any(active) else np.nan,
                ),
                "mean_abs_p_shear_active_length": scalar_or(
                    "mean_abs_p_shear_active_length",
                    float(p_shear[active].mean()) if np.any(active) else np.nan,
                ),
                "mean_chi_s_active_domain": scalar_or(
                    "mean_chi_s_active_domain",
                    float(chi[valid_chi].mean()) if np.any(valid_chi) else np.nan,
                ),
                "mean_chi_s_active_length": scalar_or(
                    "mean_chi_s_active_length",
                    float(chi[valid_chi].mean()) if np.any(valid_chi) else np.nan,
                ),
                "chi_s_p50_active_domain": scalar_or(
                    "chi_s_p50_active_domain", pixel_quantile(chi, valid_chi, 0.50)
                ),
                "chi_s_p90_active_domain": scalar_or(
                    "chi_s_p90_active_domain", pixel_quantile(chi, valid_chi, 0.90)
                ),
                "chi_s_p95_active_domain": scalar_or(
                    "chi_s_p95_active_domain", pixel_quantile(chi, valid_chi, 0.95)
                ),
                "chi_s_p99_active_domain": scalar_or(
                    "chi_s_p99_active_domain", pixel_quantile(chi, valid_chi, 0.99)
                ),
                "chi_s_p50_active_length": pixel_quantile(chi, valid_chi, 0.50),
                "chi_s_p90_active_length": pixel_quantile(chi, valid_chi, 0.90),
                "chi_s_p95_active_length": pixel_quantile(chi, valid_chi, 0.95),
                "chi_s_p99_active_length": pixel_quantile(chi, valid_chi, 0.99),
                "fraction_active_gb_length_chi_s_gt_0p8": (
                    float(np.mean(chi[valid_chi] > 0.8)) if np.any(valid_chi) else np.nan
                ),
                "fraction_active_gb_length_chi_s_near_one": scalar_or(
                    "fraction_active_gb_length_chi_s_near_one",
                    float(np.mean((chi[valid_chi] > 0.8) & (chi[valid_chi] < 1.2)))
                    if np.any(valid_chi) else np.nan,
                ),
                "fraction_active_gb_length_chi_s_gt_one": scalar_or(
                    "fraction_active_gb_length_chi_s_gt_one",
                    float(np.mean(chi[valid_chi] > 1.0)) if np.any(valid_chi) else np.nan,
                ),
                "p_net_p05_active_length": pixel_quantile(p_net, valid_p_net, 0.05),
                "p_net_p50_active_length": pixel_quantile(p_net, valid_p_net, 0.50),
                "p_net_p95_active_length": pixel_quantile(p_net, valid_p_net, 0.95),
                "fraction_active_gb_length_p_net_negative": (
                    float(np.mean(p_net[valid_p_net] < 0.0))
                    if np.any(valid_p_net) else np.nan
                ),
            })
    frame = pd.DataFrame(detail)
    rows: list[dict[str, Any]] = []
    for high, low in TOPOLOGY_WINDOWS[:3]:
        label = f"N{high}_to_{low}"
        selected = frame[frame["topology_window"].eq(label)] if not frame.empty else frame
        row: dict[str, Any] = {
            "regime": regime, "family": _family(regime), "seed": seed,
            "shear_stiffness": stiffness, "topology_window": label,
            "available": not selected.empty, "frames": int(len(selected)),
        }
        if not selected.empty:
            for column in selected.columns:
                if column != "topology_window":
                    row[column] = float(selected[column].mean())
        rows.append(row)
    return rows


def _event_triggered_trajectory(
    regime: str, seed: int, stiffness: float, trace: pd.DataFrame,
    index: pd.DataFrame,
) -> list[dict[str, Any]]:
    if trace.empty or index.empty:
        return []
    grouped = {
        str(entity): frame.drop_duplicates("step").set_index("step")
        for entity, frame in trace.groupby("entity_id", sort=False)
    }
    samples: dict[int, list[float]] = {offset: [] for offset in range(-25, 51)}
    selected = index[index["event_type"].astype(str).isin(RELEASE_TYPES)]
    for event in selected.itertuples(index=False):
        entity = grouped.get(str(event.entity_id))
        if entity is None:
            continue
        event_step = int(event.event_step)
        for offset, values in samples.items():
            step = event_step + offset
            if step in entity.index:
                values.append(float(entity.loc[step, "local_normal_velocity"]))
    rows = []
    for offset, values in samples.items():
        array = np.asarray(values, dtype=float)
        if not array.size:
            continue
        rows.append({
            "regime": regime, "family": _family(regime), "seed": seed,
            "shear_stiffness": stiffness, "offset_steps": offset,
            "samples": int(array.size), "mean_vn": float(array.mean()),
            "mean_abs_vn": float(np.abs(array).mean()),
            "median_vn": float(np.median(array)),
            "vn_p25": float(np.quantile(array, 0.25)),
            "vn_p75": float(np.quantile(array, 0.75)),
        })
    return rows


def _target_rows(
    regime: str, seed: int, stiffness: float, growth: pd.DataFrame,
) -> list[dict[str, Any]]:
    rows = []
    for target in (160, 140, 120):
        reached = growth[growth["grain_count"] <= target]
        row = {
            "regime": regime, "family": _family(regime), "seed": seed,
            "shear_stiffness": stiffness, "target_grains": target,
            "reached": not reached.empty,
        }
        if not reached.empty:
            row.update({
                "step": int(reached.iloc[0]["step"]),
                "time": float(reached.iloc[0]["time"]), "censored": False,
            })
        else:
            row.update({
                "step": int(growth.iloc[-1]["step"]),
                "time": float(growth.iloc[-1]["time"]), "censored": True,
            })
        rows.append(row)
    return rows


def _plot_lines(
    frame: pd.DataFrame, metric: str, output: Path, ylabel: str,
    *, group_extra: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    group_columns = ["family"] + ([group_extra] if group_extra else [])
    for keys, group in frame.groupby(group_columns, dropna=False):
        label = " / ".join(map(str, keys if isinstance(keys, tuple) else (keys,)))
        group = group.groupby("shear_stiffness", as_index=False).agg(
            mean=(metric, "mean"), std=(metric, "std"), count=(metric, "count")
        ).sort_values("shear_stiffness")
        ax.plot(group["shear_stiffness"], group["mean"], marker="o", label=label)
        spread = group["std"].fillna(0.0).to_numpy(float)
        if np.any(spread > 0.0):
            ax.fill_between(
                group["shear_stiffness"], group["mean"] - spread,
                group["mean"] + spread, alpha=0.14,
            )
    ax.set(xlabel="shear stiffness Ks", ylabel=ylabel)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaigns", nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--fast", action="store_true",
        help="skip event-response/null trajectories for rapid crossover selection",
    )
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    run_rows: list[dict[str, Any]] = []
    kinetic_rows: list[dict[str, Any]] = []
    force_rows: list[dict[str, Any]] = []
    episode_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []
    frame_rows: list[dict[str, Any]] = []
    response_rows: list[dict[str, Any]] = []
    causal_null_rows: list[dict[str, Any]] = []
    event_trajectory_rows: list[dict[str, Any]] = []

    seen: set[Path] = set()
    for campaign in map(Path, args.campaigns):
        for run in _runs(campaign):
            run = run.resolve()
            if run in seen:
                continue
            seen.add(run)
            manifest = json.loads((run / "manifest.json").read_text())
            config = manifest["config"]
            regime = str(config["regime"])
            if "Ks" not in regime:
                continue
            seed = int(config["seed"])
            stiffness = float(config["parameters"].get("shear_stiffness", 0.0))
            beta = float(config["parameters"].get("easy_beta", 0.35))
            growth = _growth(run)
            ledger_path = event_ledger_path(run)
            events = read_event_ledger(ledger_path) if ledger_path.exists() else pd.DataFrame()
            trace = _add_force_columns(
                _read_trace(run), beta=beta, stiffness=stiffness,
                gb_energy=float(config["pf"]["gb_energy"]),
            )
            minimal_columns = [
                "step", "entity_id", "local_normal_velocity", "shear_state_s", "tau_int",
            ]
            trace_minimal = (
                trace[minimal_columns].copy() if not trace.empty
                else pd.DataFrame(columns=minimal_columns)
            )
            run_summary = _run_summary(run, manifest, trace_minimal, events)
            run_summary.update({
                "regime": regime, "family": _family(regime), "seed": seed,
                "shear_stiffness": stiffness, "path": str(run),
                "status": manifest["status"],
                "source_sha": manifest["git_sha"],
                "end_step": int(growth.iloc[-1]["step"]),
                "end_time": float(growth.iloc[-1]["time"]),
                "end_grains": int(growth.iloc[-1]["grain_count"]),
            })
            run_rows.append(run_summary)
            for row in _window_metrics(regime, growth, events):
                row.update({
                    "family": _family(regime), "seed": seed,
                    "shear_stiffness": stiffness,
                })
                kinetic_rows.append(row)
            rows, episodes = _force_windows(
                regime, seed, stiffness, growth, trace,
                dt=float(config["pf"]["time_step"]),
            )
            force_rows.extend(rows)
            episode_rows.extend(episodes)
            target_rows.extend(_target_rows(regime, seed, stiffness, growth))
            frame_rows.extend(_frame_windows(
                regime, seed, stiffness, run,
                dx=float(config["pf"]["grid_spacing"]),
            ))
            event_index = _trace_table(run, "event_trace_events", [
                "event_id", "event_type", "sink_path", "event_step", "entity_id",
            ])
            response, causal = ([], []) if args.fast else _event_responses(
                regime, growth, trace_minimal, event_index
            )
            for row in response:
                row.update({
                    "family": _family(regime), "seed": seed,
                    "shear_stiffness": stiffness,
                })
            for row in causal:
                row.update({
                    "family": _family(regime), "seed": seed,
                    "shear_stiffness": stiffness,
                })
            response_rows.extend(response)
            causal_null_rows.extend(causal)
            if not args.fast:
                event_trajectory_rows.extend(_event_triggered_trajectory(
                    regime, seed, stiffness, trace_minimal, event_index,
                ))

    runs = pd.DataFrame(run_rows)
    kinetics = pd.DataFrame(kinetic_rows)
    forces = pd.DataFrame(force_rows)
    episodes = pd.DataFrame(episode_rows)
    targets = pd.DataFrame(target_rows)
    frame_metrics = pd.DataFrame(frame_rows)
    responses = pd.DataFrame(response_rows)
    causal_nulls = pd.DataFrame(causal_null_rows)
    event_trajectories = pd.DataFrame(event_trajectory_rows)
    runs.to_csv(output / "transition_run_summary.csv", index=False)
    kinetics.to_csv(output / "transition_kinetics.csv", index=False)
    forces.to_csv(output / "force_balance_windows.csv", index=False)
    episodes.to_csv(output / "arrest_episodes.csv", index=False)
    targets.to_csv(output / "target_times.csv", index=False)
    frame_metrics.to_csv(output / "frame_normalized_energy.csv", index=False)
    responses.to_csv(output / "solver_step_event_responses.csv", index=False)
    causal_nulls.to_csv(output / "solver_step_causal_nulls.csv", index=False)
    event_trajectories.to_csv(output / "event_triggered_vn.csv", index=False)
    censoring = targets.groupby(
        ["family", "shear_stiffness", "target_grains"], as_index=False
    ).agg(
        runs=("censored", "size"), censoring_probability=("censored", "mean"),
        median_time_or_censor=("time", "median"),
    )
    censoring.to_csv(output / "censoring_probability.csv", index=False)

    arrest_summary = (
        episodes.groupby(["regime", "family", "seed", "shear_stiffness"], as_index=False)
        .agg(
            arrest_episodes=("duration_steps", "size"),
            median_arrest_steps=("duration_steps", "median"),
            longest_arrest_steps=("duration_steps", "max"),
            median_arrest_time=("duration_time", "median"),
            longest_arrest_time=("duration_time", "max"),
            mean_entry_curvature=("entry_curvature", "mean"),
            mean_exit_curvature=("exit_curvature", "mean"),
        ) if not episodes.empty else pd.DataFrame()
    )
    arrest_summary.to_csv(output / "arrest_summary.csv", index=False)

    crossover: dict[str, Any] = {}
    for family, frame in kinetics[
        kinetics["topology_window"].eq("N160_to_140")
        & kinetics["available"].eq(True)
    ].groupby("family"):
        frame = frame.groupby("shear_stiffness", as_index=False).agg(Rdot=("Rdot", "mean"))
        crossover[str(family)] = {
            "kinetic": segmented_change_point(
                frame["shear_stiffness"].to_numpy(float),
                np.log10(np.maximum(frame["Rdot"].to_numpy(float), 1e-12)),
            )
        }
        force = frame_metrics[
            (frame_metrics["family"] == family)
            & (frame_metrics["topology_window"] == "N160_to_140")
            & frame_metrics["available"].eq(True)
        ]
        if len(force) >= 4:
            force = force.groupby("shear_stiffness", as_index=False).agg(
                fraction_length_chi_near_one=(
                    "fraction_active_gb_length_chi_s_near_one", "mean"
                )
            )
            crossover[str(family)]["force_balance"] = segmented_change_point(
                force["shear_stiffness"].to_numpy(float),
                force["fraction_length_chi_near_one"].to_numpy(float),
            )
    (output / "crossover_fits.json").write_text(json.dumps(crossover, indent=2) + "\n")

    early = kinetics[kinetics["topology_window"].isin(
        ["N190_to_160", "N160_to_140", "N140_to_120"]
    ) & kinetics["available"].eq(True)]
    _plot_lines(early, "Rdot", output / "Rdot_vs_Ks.png", "Rdot", group_extra="topology_window")
    _plot_lines(targets, "time", output / "target_time_vs_Ks.png", "time to target/censor", group_extra="target_grains")
    if not arrest_summary.empty:
        _plot_lines(arrest_summary, "longest_arrest_steps", output / "longest_arrest_vs_Ks.png", "longest arrest (steps)")
        _plot_lines(arrest_summary, "median_arrest_steps", output / "median_arrest_vs_Ks.png", "median arrest (steps)")
    mandatory_force = forces[forces["topology_window"].eq("N190_to_160") & forces["available"].eq(True)]
    for metric, name, label in (
        ("chi_p50_length", "event_trace_chi_p50_vs_Ks.png", "event-conditioned length-weighted p50 chi_s"),
        ("chi_p90_length", "event_trace_chi_p90_vs_Ks.png", "event-conditioned length-weighted p90 chi_s"),
        ("chi_p95_length", "event_trace_chi_p95_vs_Ks.png", "event-conditioned length-weighted p95 chi_s"),
        ("fraction_length_chi_gt_0p8", "event_trace_fraction_chi_gt_0p8_vs_Ks.png", "event-conditioned GB length fraction chi_s > 0.8"),
        ("fraction_length_chi_gt_one", "event_trace_fraction_chi_gt_one_vs_Ks.png", "event-conditioned GB length fraction chi_s > 1"),
        ("mean_abs_tau_length", "event_trace_mean_abs_tau_vs_Ks.png", "event-conditioned length-weighted mean abs(tau_int)"),
        ("mean_abs_p_shear_length", "event_trace_mean_abs_p_shear_vs_Ks.png", "event-conditioned length-weighted mean abs(p_shear)"),
        ("energy_per_active_length", "event_trace_energy_per_active_length_vs_Ks.png", "event-conditioned local energy, length weighted"),
    ):
        _plot_lines(mandatory_force, metric, output / name, label)
    mandatory_frames = frame_metrics[
        frame_metrics["topology_window"].eq("N190_to_160")
        & frame_metrics["available"].eq(True)
    ]
    for metric, name, label in (
        (
            "stored_shear_energy_per_active_gb_length",
            "stored_energy_per_active_gb_length_vs_Ks.png",
            "stored shear energy / active GB length",
        ),
        (
            "stored_shear_energy_per_active_domain",
            "stored_energy_per_active_domain_vs_Ks.png",
            "stored shear energy / active domain",
        ),
        (
            "active_shear_bearing_gb_fraction",
            "active_shear_bearing_fraction_vs_Ks.png",
            "active shear-bearing GB fraction",
        ),
        ("local_shear_energy_p50", "local_energy_p50_vs_Ks.png", "local shear energy p50"),
        ("local_shear_energy_p90", "local_energy_p90_vs_Ks.png", "local shear energy p90"),
        ("local_shear_energy_p95", "local_energy_p95_vs_Ks.png", "local shear energy p95"),
        ("local_shear_energy_p99", "local_energy_p99_vs_Ks.png", "local shear energy p99"),
        ("chi_s_p50_active_length", "chi_p50_vs_Ks.png", "ordinary-frame length-weighted p50 chi_s"),
        ("chi_s_p90_active_length", "chi_p90_vs_Ks.png", "ordinary-frame length-weighted p90 chi_s"),
        ("chi_s_p95_active_length", "chi_p95_vs_Ks.png", "ordinary-frame length-weighted p95 chi_s"),
        ("fraction_active_gb_length_chi_s_gt_0p8", "fraction_chi_gt_0p8_vs_Ks.png", "ordinary-frame GB length fraction chi_s > 0.8"),
        ("fraction_active_gb_length_chi_s_near_one", "fraction_chi_near_one_vs_Ks.png", "ordinary-frame GB length fraction 0.8 < chi_s < 1.2"),
        ("fraction_active_gb_length_chi_s_gt_one", "fraction_chi_gt_one_vs_Ks.png", "ordinary-frame GB length fraction chi_s > 1"),
        ("mean_abs_tau_int_active_length", "mean_abs_tau_vs_Ks.png", "ordinary-frame length-weighted mean abs(tau_int)"),
        ("mean_abs_p_shear_active_length", "mean_abs_p_shear_vs_Ks.png", "ordinary-frame length-weighted mean abs(p_shear)"),
        ("fraction_active_gb_length_p_net_negative", "fraction_p_net_negative_vs_Ks.png", "ordinary-frame GB length fraction p_net < 0"),
    ):
        if metric in mandatory_frames and mandatory_frames[metric].notna().any():
            _plot_lines(mandatory_frames, metric, output / name, label)
    for metric, name, label in (
        ("p95_abs_tauVtau_over_kBT", "tauV_over_kBT_vs_Ks.png", "p95 |tau Vtau| / kBT"),
        ("f_GBsink", "gb_sink_fraction_vs_Ks.png", "GB sink fraction"),
        ("f_TJsink", "tj_sink_fraction_vs_Ks.png", "TJ sink fraction"),
        ("G_occupancy", "G_occupancy_vs_Ks.png", "direct G occupancy"),
        ("T_occupancy", "T_occupancy_vs_Ks.png", "direct T occupancy"),
        ("C_occupancy", "C_occupancy_vs_Ks.png", "direct C occupancy"),
    ):
        if metric in runs and runs[metric].notna().any():
            _plot_lines(runs, metric, output / name, label)
    if not event_trajectories.empty:
        for family, frame in event_trajectories.groupby("family"):
            fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
            for stiffness, selected in frame.groupby("shear_stiffness"):
                curve = selected.groupby("offset_steps", as_index=False).agg(
                    mean_abs_vn=("mean_abs_vn", "mean")
                )
                ax.plot(
                    curve["offset_steps"], curve["mean_abs_vn"],
                    label=f"Ks={stiffness:g}",
                )
            ax.axvline(0, color="black", linewidth=1, alpha=0.5)
            ax.set(xlabel="solver steps from event", ylabel="mean |vn|")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8)
            fig.savefig(output / f"event_triggered_vn_{family}.png", dpi=180)
            plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
