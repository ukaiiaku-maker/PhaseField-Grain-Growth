#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grain_growth_pf.disconnections.mode import K_B_EV
from grain_growth_pf.io.event_ledger import event_ledger_path, read_event_ledger


TOPOLOGY_WINDOWS = ((190, 160), (160, 140), (140, 120), (120, 100))
TRACE_WINDOWS = (1, 2, 4, 8, 16, 32)
RELEASE_TYPES = {
    "compatibility_release", "tj_compatibility_release",
    "gb_sink_completion", "tj_sink_completion", "climb_quota_completion",
}


def _runs(root: Path) -> list[Path]:
    for name in ("video_manifest.json", "campaign_manifest.json"):
        path = root / name
        if path.exists():
            manifest = json.loads(path.read_text())
            return [
                Path(item["path"] if isinstance(item, dict) else item)
                for item in manifest["runs"]
            ]
    raise ValueError(f"no campaign manifest under {root}")


def _trace_table(run: Path, name: str, columns: list[str]) -> pd.DataFrame:
    parquet = run / f"{name}.parquet"
    csv_path = run / f"{name}.csv"
    if parquet.exists():
        return pd.read_parquet(parquet, columns=columns)
    if csv_path.exists():
        available = pd.read_csv(csv_path, nrows=0).columns
        return pd.read_csv(csv_path, usecols=[column for column in columns if column in available])
    return pd.DataFrame(columns=columns)


def _growth(run: Path) -> pd.DataFrame:
    tracks = pd.read_csv(run / "grain_tracks.csv")
    data = tracks.groupby(["step", "time"], as_index=False).agg(
        mean_radius=("radius", "mean"),
        grain_count=("grain_id", "nunique"),
    ).sort_values("step")
    if len(data) > 1:
        data["radius_rate"] = np.gradient(
            data["mean_radius"].to_numpy(float), data["time"].to_numpy(float)
        )
    else:
        data["radius_rate"] = np.nan
    return data


def _topology_label(count: float) -> str:
    if count > 190:
        return "above_190"
    for high, low in TOPOLOGY_WINDOWS:
        if low <= count <= high:
            return f"N{high}_to_{low}"
    return "below_100"


def _population_at_steps(growth: pd.DataFrame, steps: np.ndarray) -> np.ndarray:
    return np.interp(
        steps.astype(float), growth["step"].to_numpy(float),
        growth["grain_count"].to_numpy(float),
    )


def _window_metrics(
    regime: str, growth: pd.DataFrame, events: pd.DataFrame
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for high, low in TOPOLOGY_WINDOWS:
        label = f"N{high}_to_{low}"
        frame = growth[(growth["grain_count"] <= high) & (growth["grain_count"] >= low)]
        available = bool(
            len(frame) >= 3
            and growth["grain_count"].max() >= high
            and growth["grain_count"].min() <= low
        )
        row: dict[str, Any] = {
            "regime": regime, "topology_window": label, "available": available,
            "points": int(len(frame)),
        }
        if len(frame) >= 3:
            time = frame["time"].to_numpy(float)
            radius = frame["mean_radius"].to_numpy(float)
            rate = np.diff(radius) / np.diff(time)
            finite = rate[np.isfinite(rate)]
            mean_rate = float(np.polyfit(time, radius, 1)[0])
            scale = float(np.median(np.abs(finite))) if finite.size else np.nan
            stationary_threshold = max(1e-12, 0.1 * scale) if np.isfinite(scale) else np.nan
            positive_time = time > 0.0
            exponent = (
                float(np.polyfit(np.log(time[positive_time]), np.log(radius[positive_time]), 1)[0])
                if positive_time.sum() >= 3 else np.nan
            )
            row.update({
                "time_start": float(time[0]), "time_end": float(time[-1]),
                "grain_count_start": int(frame["grain_count"].iloc[0]),
                "grain_count_end": int(frame["grain_count"].iloc[-1]),
                "Rdot": mean_rate,
                "CV_jerk": float(np.std(finite, ddof=1) / abs(np.mean(finite)))
                if finite.size > 1 and np.mean(finite) != 0.0 else np.nan,
                "stationary_fraction": float(np.mean(np.abs(finite) <= stationary_threshold))
                if finite.size else np.nan,
                "local_effective_exponent": exponent,
            })
            if not events.empty and "step" in events:
                selected = events[
                    (pd.to_numeric(events["step"], errors="coerce") >= frame["step"].min())
                    & (pd.to_numeric(events["step"], errors="coerce") <= frame["step"].max())
                ]
                types = selected.get("event_type", pd.Series(dtype=str)).astype(str)
                g_count = int((types == "compatibility_release").sum())
                t_count = int((types == "tj_compatibility_release").sum())
                c_count = int(types.isin({"gb_sink_completion", "tj_sink_completion", "climb_quota_completion"}).sum())
                total = g_count + t_count + c_count
                row.update({
                    "G_releases": g_count, "T_releases": t_count, "C_releases": c_count,
                    "f_G": g_count / total if total else np.nan,
                    "f_T": t_count / total if total else np.nan,
                    "f_C": c_count / total if total else np.nan,
                })
        rows.append(row)
    return rows


def _occupancy(run: Path) -> dict[str, float]:
    path = run / "mechanism_state.csv"
    if not path.exists():
        return {}
    state = pd.read_csv(path)
    if state.empty:
        return {}
    step_times = state.groupby("step")["time"].first().sort_index()
    weights = step_times.diff().shift(-1)
    fallback = float(weights.dropna().median()) if weights.notna().any() else 1.0
    mapped = state["step"].map(weights.fillna(fallback)).to_numpy(float)
    denominator = float(mapped.sum())
    return {
        "G_occupancy": float(np.sum(mapped * state["G_pending"]) / denominator),
        "T_occupancy": float(np.sum(mapped * state["T_pending"]) / denominator),
        "C_occupancy": float(np.sum(mapped * state["C_pending"]) / denominator),
        "blocked_occupancy": float(np.sum(mapped * state["blocked"]) / denominator),
    }


def _event_responses(
    regime: str, growth: pd.DataFrame, trace: pd.DataFrame, index: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    detail: list[dict[str, Any]] = []
    null_rows: list[dict[str, Any]] = []
    if trace.empty or index.empty:
        return detail, null_rows
    trace = trace.sort_values(["entity_id", "step"]).copy()
    trace["topology_window"] = [
        _topology_label(value)
        for value in _population_at_steps(growth, trace["step"].to_numpy(int))
    ]
    grouped = {
        str(entity): frame.set_index("step", drop=False)
        for entity, frame in trace.groupby("entity_id", sort=False)
    }
    selected_events = index[index["event_type"].astype(str).isin(RELEASE_TYPES)]
    for event in selected_events.itertuples(index=False):
        entity = str(event.entity_id)
        frame = grouped.get(entity)
        if frame is None:
            continue
        event_step = int(event.event_step)
        for window in TRACE_WINDOWS:
            if event_step - window not in frame.index or event_step + window not in frame.index:
                continue
            before = frame.loc[event_step - window]
            after = frame.loc[event_step + window]
            if isinstance(before, pd.DataFrame) or isinstance(after, pd.DataFrame):
                continue
            detail.append({
                "regime": regime, "event_id": event.event_id,
                "event_type": event.event_type, "sink_path": event.sink_path,
                "entity_id": entity, "event_step": event_step,
                "window_steps": window,
                "topology_window": _topology_label(float(_population_at_steps(
                    growth, np.asarray([event_step])
                )[0])),
                "signed_velocity_response": float(after.local_normal_velocity - before.local_normal_velocity),
                "absolute_velocity_response": float(abs(after.local_normal_velocity) - abs(before.local_normal_velocity)),
                "shear_state_relaxation": float(abs(after.shear_state_s) - abs(before.shear_state_s)),
                "tau_int_relaxation": float(abs(after.tau_int) - abs(before.tau_int)),
            })

    response = pd.DataFrame(detail)
    if response.empty:
        return detail, null_rows
    rng = np.random.default_rng(20260823)
    for window in TRACE_WINDOWS:
        actual = response[response["window_steps"] == window]
        pools: list[pd.DataFrame] = []
        for entity, frame in grouped.items():
            values = frame[["step", "local_normal_velocity", "topology_window"]].copy()
            values["response"] = (
                values["local_normal_velocity"].shift(-window).abs()
                - values["local_normal_velocity"].shift(window).abs()
            )
            values["entity_id"] = entity
            pools.append(values.dropna(subset=["response"]))
        pool = pd.concat(pools, ignore_index=True) if pools else pd.DataFrame()
        if actual.empty or pool.empty:
            continue
        entity_pools = {
            str(entity): frame["response"].to_numpy(float)
            for entity, frame in pool.groupby("entity_id", sort=False)
        }
        topology_pools = {
            str(label): frame["response"].to_numpy(float)
            for label, frame in pool.groupby("topology_window", sort=False)
        }
        same_entity: list[float] = []
        same_topology: list[float] = []
        for event in actual.itertuples(index=False):
            entity_pool = entity_pools.get(str(event.entity_id), np.asarray([]))
            topology_pool = topology_pools.get(str(event.topology_window), np.asarray([]))
            if entity_pool.size:
                same_entity.append(float(entity_pool[int(rng.integers(entity_pool.size))]))
            if topology_pool.size:
                same_topology.append(float(topology_pool[int(rng.integers(topology_pool.size))]))

        circular_scores: list[float] = []
        block_scores: list[float] = []
        event_counts = {
            str(entity): int(count)
            for entity, count in actual["entity_id"].value_counts().items()
        }
        for _ in range(200):
            circular_sample: list[float] = []
            block_sample: list[float] = []
            for entity, count in event_counts.items():
                entity_pool = entity_pools.get(entity, np.asarray([]))
                if not entity_pool.size:
                    continue
                offset = int(rng.integers(1, entity_pool.size + 1))
                event_indices = np.linspace(0, entity_pool.size - 1, count, dtype=int)
                circular_sample.extend(entity_pool[(event_indices + offset) % entity_pool.size])
                block_size = max(8, 2 * window)
                blocks = [
                    entity_pool[start:start + block_size]
                    for start in range(0, entity_pool.size, block_size)
                ]
                permuted = np.concatenate([blocks[i] for i in rng.permutation(len(blocks))])
                block_sample.extend(permuted[:count].tolist())
            if circular_sample:
                circular_scores.append(float(np.mean(circular_sample)))
            if block_sample:
                block_scores.append(float(np.mean(block_sample)))
        actual_mean = float(actual["absolute_velocity_response"].mean())
        circular_mean = float(np.mean(circular_scores)) if circular_scores else np.nan
        block_mean = float(np.mean(block_scores)) if block_scores else np.nan
        entity_mean = float(np.mean(same_entity)) if same_entity else np.nan
        topology_mean = float(np.mean(same_topology)) if same_topology else np.nan
        null_rows.append({
            "regime": regime, "window_steps": window, "events": int(len(actual)),
            "actual_mean_absolute_velocity_response": actual_mean,
            "circular_shift_null_mean": circular_mean,
            "circular_shift_enrichment": actual_mean - circular_mean,
            "block_shuffle_null_mean": block_mean,
            "block_shuffle_enrichment": actual_mean - block_mean,
            "same_entity_null_mean": entity_mean,
            "same_entity_enrichment": actual_mean - entity_mean,
            "same_topology_window_null_mean": topology_mean,
            "same_topology_window_enrichment": actual_mean - topology_mean,
        })
    return detail, null_rows


def _run_summary(
    run: Path, manifest: dict[str, Any], trace: pd.DataFrame, events: pd.DataFrame,
) -> dict[str, Any]:
    config = manifest["config"]
    regime = str(config["regime"])
    temperature = float(config["pf"]["temperature"])
    row: dict[str, Any] = {
        "regime": regime, "path": str(run), "status": manifest["status"],
        "source_sha": manifest["git_sha"], "temperature": temperature,
        "seed": int(config["seed"]),
        "shear_stiffness": float(config["parameters"].get("shear_stiffness", 0.0)),
        "steps": int(manifest.get("steps_completed", 0)),
        "final_grains": int(manifest.get("final_grains", -1)),
    }
    row.update(_occupancy(run))
    frame_rms: list[float] = []
    frame_tau_max: list[float] = []
    for frame_path in sorted((run / "frames").glob("frame-*.npz")):
        with np.load(frame_path) as frame:
            if "boundary_shear_state_rms" in frame:
                frame_rms.append(float(frame["boundary_shear_state_rms"]))
            if "shear_stress_max_abs" in frame:
                frame_tau_max.append(float(frame["shear_stress_max_abs"]))
    if frame_rms:
        row["RMS_shear_state"] = float(np.sqrt(np.mean(np.square(frame_rms))))
        row["max_frame_RMS_shear_state"] = float(max(frame_rms))
    if frame_tau_max:
        row["max_frame_abs_tau_int"] = float(max(frame_tau_max))
    inventory_path = run / "defect_inventory.csv"
    if inventory_path.exists():
        final = pd.read_csv(inventory_path).iloc[-1]
        accommodated = float(final["N_accommodated_GB"] + final["N_accommodated_TJ"])
        row.update({
            "N_required": float(final["N_required"]),
            "N_accommodated_GB": float(final["N_accommodated_GB"]),
            "N_accommodated_TJ": float(final["N_accommodated_TJ"]),
            "f_GBsink": float(final["N_accommodated_GB"] / accommodated) if accommodated else np.nan,
            "f_TJsink": float(final["N_accommodated_TJ"] / accommodated) if accommodated else np.nan,
            "N_stored": float(final["N_active_deficit"] + final["N_retired"]),
            "conservation_residual": float(final["conservation_residual"]),
            "max_abs_conservation_residual": float(final["max_abs_conservation_residual"]),
        })
    if not trace.empty:
        shear = pd.to_numeric(trace["shear_state_s"], errors="coerce").dropna().to_numpy(float)
        tau = pd.to_numeric(trace["tau_int"], errors="coerce").abs().dropna().to_numpy(float)
        row.update({
            "event_trace_RMS_shear_state": float(np.sqrt(np.mean(shear**2))) if shear.size else np.nan,
            "max_abs_tau_int": float(np.max(tau)) if tau.size else np.nan,
            "p95_abs_tau_int": float(np.quantile(tau, 0.95)) if tau.size else np.nan,
            "p99_abs_tau_int": float(np.quantile(tau, 0.99)) if tau.size else np.nan,
        })
    energy_path = run / "energy.json"
    if energy_path.exists():
        energy = json.loads(energy_path.read_text())
        if energy:
            row["stored_shear_energy_final"] = float(energy[-1].get("stored_shear", 0.0))
            row["stored_shear_energy_max"] = float(max(item.get("stored_shear", 0.0) for item in energy))
    work_path = run / "activation_work.csv"
    if work_path.exists():
        work = pd.read_csv(work_path)
        normalized = pd.to_numeric(work.get("work_shear"), errors="coerce").abs() / (K_B_EV * temperature)
        normalized = normalized.replace([np.inf, -np.inf], np.nan).dropna()
        for percentile in (50, 90, 95, 99):
            row[f"p{percentile}_abs_tauVtau_over_kBT"] = (
                float(np.percentile(normalized, percentile)) if len(normalized) else np.nan
            )
        for factor in (2, 5, 10):
            row[f"fraction_events_changed_gt_{factor}x_by_shear"] = (
                float(np.mean(normalized > np.log(factor))) if len(normalized) else np.nan
            )
    row["release_events"] = int(
        events.get("event_type", pd.Series(dtype=str)).astype(str).isin(RELEASE_TYPES).sum()
    )
    return row


def _plot_responses(summary: pd.DataFrame, windows: pd.DataFrame, output: Path) -> None:
    metrics = [
        "RMS_shear_state", "p95_abs_tau_int", "stored_shear_energy_max",
        "p95_abs_tauVtau_over_kBT", "fraction_events_changed_gt_2x_by_shear",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    for ax, metric in zip(axes.flat, metrics, strict=False):
        for family, frame in summary[summary["regime"].str.contains("Ks")].groupby(
            summary["regime"].str.replace(r"_Ks\d+$", "", regex=True)
        ):
            frame = frame.sort_values("shear_stiffness")
            ax.plot(frame["shear_stiffness"], frame[metric], marker="o", label=family)
        ax.set(xlabel="shear stiffness Ks", ylabel=metric)
        ax.grid(alpha=0.25)
    early = windows[(windows["topology_window"] == "N190_to_160") & windows["regime"].str.contains("Ks")]
    ax = axes.flat[-1]
    for family, frame in early.groupby(early["regime"].str.replace(r"_Ks\d+$", "", regex=True)):
        frame = frame.merge(summary[["regime", "shear_stiffness"]], on="regime")
        frame = frame.sort_values("shear_stiffness")
        ax.plot(frame["shear_stiffness"], frame["Rdot"], marker="o", label=family)
    ax.set(xlabel="shear stiffness Ks", ylabel="Rdot, N=190→160")
    ax.grid(alpha=0.25)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3)
    fig.savefig(output / "shear_response_vs_Ks.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = Path(args.campaign)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, Any]] = []
    windows: list[dict[str, Any]] = []
    responses: list[dict[str, Any]] = []
    nulls: list[dict[str, Any]] = []
    for run in _runs(root):
        manifest = json.loads((run / "manifest.json").read_text())
        regime = str(manifest["config"]["regime"])
        growth = _growth(run)
        events_path = event_ledger_path(run)
        events = read_event_ledger(events_path) if events_path.exists() else pd.DataFrame()
        trace = _trace_table(run, "event_traces", [
            "step", "entity_id", "local_normal_velocity", "shear_state_s", "tau_int",
        ])
        event_index = _trace_table(run, "event_trace_events", [
            "event_id", "event_type", "sink_path", "event_step", "entity_id",
        ])
        summaries.append(_run_summary(run, manifest, trace, events))
        summaries[-1].update({
            "event_trace_rows": int(len(trace)),
            "event_trace_events": int(len(event_index)),
            "duplicate_entity_step_trace_rows": int(
                trace.duplicated(["entity_id", "step"]).sum()
            ) if not trace.empty else 0,
        })
        windows.extend(_window_metrics(regime, growth, events))
        detail, causal = _event_responses(regime, growth, trace, event_index)
        responses.extend(detail)
        nulls.extend(causal)

    summary_frame = pd.DataFrame(summaries)
    window_frame = pd.DataFrame(windows)
    summary_frame.to_csv(output / "shear_run_summary.csv", index=False)
    window_frame.to_csv(output / "topology_window_metrics.csv", index=False)
    pd.DataFrame(responses).to_csv(output / "solver_step_event_responses.csv", index=False)
    pd.DataFrame(nulls).to_csv(output / "solver_step_causal_nulls.csv", index=False)
    _plot_responses(summary_frame, window_frame, output)
    (output / "shear_analysis_manifest.json").write_text(json.dumps({
        "campaign": str(root),
        "topology_windows": [list(window) for window in TOPOLOGY_WINDOWS],
        "event_response_windows_steps": list(TRACE_WINDOWS),
        "causal_nulls": [
            "circular_shift", "block_shuffle", "same_entity", "same_topology_window",
        ],
        "outputs": [
            "shear_run_summary.csv", "topology_window_metrics.csv",
            "solver_step_event_responses.csv", "solver_step_causal_nulls.csv",
            "shear_response_vs_Ks.png",
        ],
    }, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
