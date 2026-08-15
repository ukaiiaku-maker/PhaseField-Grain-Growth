from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from grain_growth_pf.analysis.campaign import _fit_window, analyze_campaign
from grain_growth_pf.analysis.grain_tracks import ensemble_radius, load_tracks


def _save(fig: plt.Figure, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(target.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(target.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _aligned_ensemble(paths: list[Path]) -> pd.DataFrame:
    aligned = None
    for index, path in enumerate(paths):
        data = ensemble_radius(load_tracks(path / "grain_tracks.csv"))
        data = data[["time", "R_A", "grain_count"]].rename(
            columns={"R_A": f"R_{index}", "grain_count": f"N_{index}"}
        )
        aligned = data if aligned is None else aligned.merge(data, on="time", how="inner")
    assert aligned is not None
    aligned = aligned.sort_values("time").reset_index(drop=True)
    r_columns = [column for column in aligned if column.startswith("R_")]
    n_columns = [column for column in aligned if column.startswith("N_")]
    aligned["R_mean"] = aligned[r_columns].mean(axis=1)
    aligned["R_std"] = aligned[r_columns].std(axis=1, ddof=1).fillna(0.0)
    aligned["N_mean"] = aligned[n_columns].mean(axis=1)
    return aligned


def _local_exponent(time: np.ndarray, radius: np.ndarray, half_window: int = 5) -> np.ndarray:
    """Return the best local linearizing power on a bounded profile grid."""
    result = np.full(len(time), np.nan)
    exponents = np.linspace(1.0, 6.0, 101)
    for center in range(half_window, len(time) - half_window):
        selection = slice(center - half_window, center + half_window + 1)
        local_time, local_radius = time[selection], radius[selection]
        errors = []
        for exponent in exponents:
            transformed = local_radius**exponent
            fitted = np.maximum(
                np.polyval(np.polyfit(local_time, transformed, 1), local_time), 1e-30
            ) ** (1.0 / exponent)
            errors.append(
                np.sqrt(np.mean((local_radius - fitted) ** 2))
                / max(np.std(local_radius), 1e-15)
            )
        result[center] = exponents[int(np.argmin(errors))]
    return result


def _kinetics_figure(paths: list[Path], row: pd.Series, target: Path) -> None:
    data = _aligned_ensemble(paths)
    time = data["time"].to_numpy(float)
    radius = data["R_mean"].to_numpy(float)
    spread = data["R_std"].to_numpy(float)
    exponent = float(row["n"])
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    for power, axis, label in zip((1, 2, 3), axes.flat[:3], ("R", r"$R^2$", r"$R^3$")):
        axis.plot(time, radius**power, color="C0", lw=1.8)
        if power == 1:
            axis.fill_between(time, radius - spread, radius + spread, color="C0", alpha=0.2)
        axis.set_ylabel(label)
    transformed = radius**exponent
    axes[1, 0].plot(time, transformed, color="C1", label=f"n={exponent:.2f}")
    start, end, _ = _fit_window(data["N_mean"].to_numpy(float))
    fit_time, fit_values = time[start:end], transformed[start:end]
    coefficient, intercept = np.polyfit(fit_time, fit_values, 1)
    axes[1, 0].plot(fit_time, coefficient * fit_time + intercept, "k--", lw=1, label="fit window")
    axes[1, 0].set_ylabel(r"$R^n$")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].plot(time, _local_exponent(time, radius), color="C2")
    axes[1, 1].axhline(2, color="0.4", ls="--", lw=1)
    axes[1, 1].set_ylabel("local effective n")
    axes[1, 2].plot(time, data["N_mean"], color="C3")
    axes[1, 2].set_ylabel("grain count")
    for axis in axes[1]:
        axis.set_xlabel("time")
    fig.suptitle(f"{row['regime']} at {row['temperature']:g} K — {len(paths)} realizations")
    _save(fig, target)


def _representative_figure(path: Path, target: Path) -> None:
    tracks = load_tracks(path / "grain_tracks.csv")
    counts = tracks.groupby("grain_id").size().sort_values(ascending=False)
    selected = counts.head(8).index
    fig, (area_axis, radius_axis, neighbor_axis) = plt.subplots(1, 3, figsize=(13, 4))
    for grain_id in selected:
        grain = tracks[tracks["grain_id"] == grain_id].sort_values("time")
        area_axis.plot(grain["time"], grain["area"], lw=1, label=str(grain_id))
        radius_axis.plot(grain["time"], grain["radius"], lw=1)
        if len(grain) > 1:
            rate = np.diff(grain["area"].to_numpy(float)) / np.diff(grain["time"].to_numpy(float))
            neighbor_axis.scatter(grain["neighbors"].to_numpy(float)[:-1], rate, s=7, alpha=0.35)
    event_path = path / "events.csv"
    if event_path.exists() and event_path.stat().st_size:
        events = pd.read_csv(event_path)
        if not events.empty and "time" in events:
            for event_time in events["time"].to_numpy(float)[:200]:
                area_axis.axvline(event_time, color="k", alpha=0.035, lw=0.6)
                radius_axis.axvline(event_time, color="k", alpha=0.035, lw=0.6)
    area_axis.set_ylabel("grain area")
    radius_axis.set_ylabel("equivalent radius")
    radius_axis.set_xlabel("time")
    neighbor_axis.set(xlabel="neighbor number", ylabel="area growth rate")
    area_axis.legend(title="grain", ncol=4, frameon=False, fontsize=7)
    fig.suptitle(path.name)
    _save(fig, target)


def _boundary_figure(paths: list[Path], target: Path) -> None:
    frames = []
    for path in paths:
        boundary_path = path / "boundary_tracks.csv"
        if boundary_path.exists() and boundary_path.stat().st_size:
            frames.append(pd.read_csv(boundary_path))
    if not frames:
        return
    boundaries = pd.concat(frames, ignore_index=True).dropna(subset=["curvature", "normal_velocity"])
    if boundaries.empty:
        return
    if len(boundaries) > 30000:
        boundaries = boundaries.sample(30000, random_state=1729)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    colors = np.where(boundaries["blocked"].to_numpy(bool), "C3", "C0")
    axes[0].scatter(boundaries["curvature"], boundaries["normal_velocity"], c=colors, s=5, alpha=0.2)
    axes[0].axhline(0, color="0.5", lw=0.7)
    axes[0].axvline(0, color="0.5", lw=0.7)
    axes[0].set(xlabel="signed curvature", ylabel="signed normal velocity")
    axes[1].hist(boundaries["normal_velocity"], bins=60, color="C0", alpha=0.8)
    axes[1].set(xlabel="normal velocity", ylabel="count")
    _save(fig, target)


def _event_figure(paths: list[Path], target: Path) -> None:
    frames = []
    for path in paths:
        event_path = path / "events.csv"
        if event_path.exists() and event_path.stat().st_size:
            frame = pd.read_csv(event_path)
            if not frame.empty:
                frames.append(frame)
    if not frames:
        return
    events = pd.concat(frames, ignore_index=True)
    times = np.sort(events["time"].to_numpy(float))
    waits = np.diff(times)
    sizes = events.get("packet_size", pd.Series(np.ones(len(events)))).to_numpy(float)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    axes[0].hist(waits[waits > 0], bins=35, density=True, color="C0", alpha=0.8)
    axes[0].set(xlabel="waiting time", ylabel="density")
    axes[1].hist(sizes[np.isfinite(sizes)], bins=25, color="C1", alpha=0.8)
    axes[1].set(xlabel="packet size", ylabel="count")
    ordered = np.sort(sizes[np.isfinite(sizes) & (sizes > 0)])
    if len(ordered):
        axes[2].loglog(ordered, 1.0 - np.arange(len(ordered)) / len(ordered), color="C2")
    axes[2].set(xlabel="event packet/burst size", ylabel="CCDF")
    _save(fig, target)


def plot_campaign(campaign_dir: str | Path, output_dir: str | Path | None = None) -> Path:
    campaign_dir = Path(campaign_dir)
    output = Path(output_dir) if output_dir else campaign_dir / "plots"
    summary_path = campaign_dir / "mechanism_summary.csv"
    summary = pd.read_csv(summary_path) if summary_path.exists() else analyze_campaign(campaign_dir)
    campaign = json.loads((campaign_dir / "campaign_manifest.json").read_text())
    grouped: dict[tuple[str, float], list[Path]] = {}
    for raw in campaign["runs"]:
        path = Path(raw)
        manifest = json.loads((path / "manifest.json").read_text())
        config = manifest["config"]
        grouped.setdefault((config["regime"], float(config["pf"]["temperature"])), []).append(path)
    for key, paths in grouped.items():
        row = summary[(summary["regime"] == key[0]) & (summary["temperature"] == key[1])].iloc[0]
        stem = f"{key[0]}-T{key[1]:g}"
        _kinetics_figure(paths, row, output / f"{stem}-ensemble-kinetics")
        _representative_figure(paths[0], output / f"{stem}-representative-grains")
        _boundary_figure(paths, output / f"{stem}-velocity-curvature")
        _event_figure(paths, output / f"{stem}-event-statistics")

    fig, axis = plt.subplots(figsize=(7, 5))
    for (regime, temperature), paths in grouped.items():
        data = _aligned_ensemble(paths)
        axis.plot(data["time"], data["R_mean"], label=f"{regime}, {temperature:g} K")
    axis.set(xlabel="time", ylabel=r"$R_A$", title="Mechanism comparison")
    axis.legend(frameon=False, fontsize=7, ncol=2)
    _save(fig, output / "mechanism-comparison")

    for regime, group in summary.groupby("regime"):
        if len(group) < 4 or np.any(group["K"] <= 0):
            continue
        ordered = group.sort_values("temperature")
        fig, axis = plt.subplots(figsize=(6, 4))
        axis.errorbar(1.0 / ordered["temperature"], np.log(ordered["K"]),
                      yerr=ordered["K_ci"] / ordered["K"], marker="o", capsize=3)
        axis.set(xlabel=r"$1/T$ (K$^{-1}$)", ylabel=r"$\ln K_n$", title=f"{regime} Arrhenius scaling")
        _save(fig, output / f"{regime}-arrhenius")
    return output
