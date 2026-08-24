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

from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger

from analyze_shear_screen import (
    TOPOLOGY_WINDOWS,
    _growth,
    _population_at_steps,
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
    if parquet.exists():
        frame = pd.read_parquet(parquet)
    elif csv_path.exists():
        frame = pd.read_csv(csv_path)
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
        "local_shear_energy",
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
    trace: pd.DataFrame,
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
            "fraction_length_chi_gt_0p8": float(weights[chi > 0.8].sum() / total_weight),
            "fraction_length_chi_near_one": float(
                weights[(chi > 0.8) & (chi < 1.2)].sum() / total_weight
            ),
            "fraction_length_chi_gt_one": float(weights[chi > 1.0].sum() / total_weight),
            "p_net_p05_length": weighted_quantile(valid["p_net"], weights, 0.05),
            "p_net_p50_length": weighted_quantile(valid["p_net"], weights, 0.50),
            "p_net_p95_length": weighted_quantile(valid["p_net"], weights, 0.95),
            "p_net_sign_change_fraction": sign_changes / sign_pairs if sign_pairs else np.nan,
            "energy_per_active_domain": float(np.mean(energy)),
            "energy_per_active_length": float(np.average(energy, weights=weights)),
            "local_energy_p50": float(np.quantile(energy, 0.50)),
            "local_energy_p90": float(np.quantile(energy, 0.90)),
            "local_energy_p95": float(np.quantile(energy, 0.95)),
            "local_energy_p99": float(np.quantile(energy, 0.99)),
        })

    dt = float(np.median(np.diff(np.sort(trace["time"].dropna().unique()))))
    if not np.isfinite(dt) or dt <= 0.0:
        dt = 0.04
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
        group = group.sort_values("shear_stiffness")
        ax.plot(group["shear_stiffness"], group[metric], marker="o", label=label)
    ax.set(xlabel="shear stiffness Ks", ylabel=ylabel)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaigns", nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    run_rows: list[dict[str, Any]] = []
    kinetic_rows: list[dict[str, Any]] = []
    force_rows: list[dict[str, Any]] = []
    episode_rows: list[dict[str, Any]] = []
    target_rows: list[dict[str, Any]] = []

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
            run_rows.append({
                "regime": regime, "family": _family(regime), "seed": seed,
                "shear_stiffness": stiffness, "path": str(run),
                "status": manifest["status"],
                "source_sha": manifest["git_sha"],
                "end_step": int(growth.iloc[-1]["step"]),
                "end_time": float(growth.iloc[-1]["time"]),
                "end_grains": int(growth.iloc[-1]["grain_count"]),
            })
            for row in _window_metrics(regime, growth, events):
                row.update({
                    "family": _family(regime), "seed": seed,
                    "shear_stiffness": stiffness,
                })
                kinetic_rows.append(row)
            rows, episodes = _force_windows(regime, seed, stiffness, growth, trace)
            force_rows.extend(rows)
            episode_rows.extend(episodes)
            target_rows.extend(_target_rows(regime, seed, stiffness, growth))

    runs = pd.DataFrame(run_rows)
    kinetics = pd.DataFrame(kinetic_rows)
    forces = pd.DataFrame(force_rows)
    episodes = pd.DataFrame(episode_rows)
    targets = pd.DataFrame(target_rows)
    runs.to_csv(output / "transition_run_summary.csv", index=False)
    kinetics.to_csv(output / "transition_kinetics.csv", index=False)
    forces.to_csv(output / "force_balance_windows.csv", index=False)
    episodes.to_csv(output / "arrest_episodes.csv", index=False)
    targets.to_csv(output / "target_times.csv", index=False)

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
        crossover[str(family)] = {
            "kinetic": segmented_change_point(
                frame["shear_stiffness"].to_numpy(float),
                np.log10(np.maximum(frame["Rdot"].to_numpy(float), 1e-12)),
            )
        }
        force = forces[(forces["family"] == family) & (forces["topology_window"] == "N160_to_140")]
        if len(force) >= 4:
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
        ("chi_p50_length", "chi_p50_vs_Ks.png", "length-weighted p50 chi_s"),
        ("chi_p90_length", "chi_p90_vs_Ks.png", "length-weighted p90 chi_s"),
        ("chi_p95_length", "chi_p95_vs_Ks.png", "length-weighted p95 chi_s"),
        ("fraction_length_chi_gt_0p8", "fraction_chi_gt_0p8_vs_Ks.png", "GB length fraction chi_s > 0.8"),
        ("fraction_length_chi_gt_one", "fraction_chi_gt_one_vs_Ks.png", "GB length fraction chi_s > 1"),
        ("mean_abs_tau_length", "mean_abs_tau_vs_Ks.png", "length-weighted mean abs(tau_int)"),
        ("mean_abs_p_shear_length", "mean_abs_p_shear_vs_Ks.png", "length-weighted mean abs(p_shear)"),
        ("energy_per_active_length", "energy_per_active_length_vs_Ks.png", "local energy, length weighted"),
    ):
        _plot_lines(mandatory_force, metric, output / name, label)
    print(output)


if __name__ == "__main__":
    main()
