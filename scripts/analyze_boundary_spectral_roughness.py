#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

TARGET_GRAINS = (175, 150, 125, 100)


def boundary_mask(labels: np.ndarray) -> np.ndarray:
    return (
        (labels != np.roll(labels, 1, axis=0))
        | (labels != np.roll(labels, -1, axis=0))
        | (labels != np.roll(labels, 1, axis=1))
        | (labels != np.roll(labels, -1, axis=1))
    )


def spectral_metrics(labels: np.ndarray) -> tuple[float, float]:
    mask = boundary_mask(labels).astype(float)
    field = mask - mask.mean()
    power = np.abs(np.fft.fft2(field)) ** 2
    power[0, 0] = 0.0
    ky = np.fft.fftfreq(labels.shape[0])[:, None]
    kx = np.fft.fftfreq(labels.shape[1])[None, :]
    kr = np.sqrt(kx * kx + ky * ky)
    total = float(power.sum())
    if total <= 0:
        return np.nan, np.nan
    high_fraction = float(power[kr >= 0.25].sum() / total)
    centroid = float((power * kr).sum() / total)
    return high_fraction, centroid


def compactness(labels: np.ndarray) -> float:
    values = []
    for grain in np.unique(labels):
        phase = labels == grain
        area = float(np.count_nonzero(phase))
        if area <= 0:
            continue
        perimeter = float(
            np.count_nonzero(phase & ~np.roll(phase, 1, axis=0))
            + np.count_nonzero(phase & ~np.roll(phase, -1, axis=0))
            + np.count_nonzero(phase & ~np.roll(phase, 1, axis=1))
            + np.count_nonzero(phase & ~np.roll(phase, -1, axis=1))
        )
        if perimeter > 0:
            values.append(perimeter * perimeter / (4.0 * np.pi * area))
    return float(np.mean(values)) if values else np.nan


def analyze_run(run: Path) -> pd.DataFrame:
    rows = []
    for frame in sorted((run / "frames").glob("frame-*.npz")):
        with np.load(frame) as data:
            labels = data["labels"].astype(np.int64)
            step = int(data["step"])
            time = float(data["time"])
        high, centroid = spectral_metrics(labels)
        rows.append({
            "run": run.name,
            "frame": frame.name,
            "step": step,
            "time": time,
            "grain_count": int(np.unique(labels).size),
            "q_pixel": compactness(labels),
            "high_k_boundary_power_fraction": high,
            "boundary_spectral_centroid": centroid,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("video_campaign")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    root = Path(args.video_campaign)
    all_frames = []
    for run in sorted(root.iterdir()):
        if run.is_dir() and (run / "frames").exists():
            frame = analyze_run(run)
            if not frame.empty:
                all_frames.append(frame)
    if not all_frames:
        raise SystemExit(f"no frame directories found under {root}")
    history = pd.concat(all_frames, ignore_index=True)

    selected = []
    for run, frame in history.groupby("run"):
        for target in TARGET_GRAINS:
            index = (frame["grain_count"] - target).abs().idxmin()
            row = history.loc[index].to_dict()
            row["target_grains"] = target
            selected.append(row)
    table = pd.DataFrame(selected)

    b0_names = [name for name in table["run"].unique() if name.startswith("B0")]
    if b0_names:
        b0 = table[table["run"] == b0_names[0]].set_index("target_grains")
        for i, row in table.iterrows():
            target = int(row["target_grains"])
            if target not in b0.index:
                continue
            ref = b0.loc[target]
            table.loc[i, "q_over_B0"] = row["q_pixel"] / ref["q_pixel"]
            table.loc[i, "high_k_over_B0"] = (
                row["high_k_boundary_power_fraction"] / ref["high_k_boundary_power_fraction"]
            )
            table.loc[i, "spectral_centroid_over_B0"] = (
                row["boundary_spectral_centroid"] / ref["boundary_spectral_centroid"]
            )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.sort_values(["target_grains", "run"]).to_csv(output, index=False)
    history.to_csv(output.with_name(output.stem + "_all_frames.csv"), index=False)

    report = output.with_suffix(".md")
    lines = [
        "# Boundary spectral-roughness audit",
        "",
        "The high-k metric is the fraction of boundary-indicator FFT power above 0.25 cycles/pixel. It is a pixel-scale waviness diagnostic complementary to P^2/(4 pi A).",
        "",
    ]
    for target in TARGET_GRAINS:
        lines.append(f"## target N={target}")
        sub = table[table["target_grains"] == target].sort_values("run")
        for _, row in sub.iterrows():
            lines.append(
                f"- {row['run']}: actual N={int(row['grain_count'])}, "
                f"q/B0={row.get('q_over_B0', np.nan):.4f}, "
                f"high-k/B0={row.get('high_k_over_B0', np.nan):.4f}, "
                f"centroid/B0={row.get('spectral_centroid_over_B0', np.nan):.4f}"
            )
        lines.append("")
    report.write_text("\n".join(lines))
    print(table.sort_values(["target_grains", "run"]).to_string(index=False))
    print(report)


if __name__ == "__main__":
    main()
