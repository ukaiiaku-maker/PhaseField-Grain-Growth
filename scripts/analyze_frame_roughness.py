#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def frame_metrics(path: Path, minimum_area: int) -> dict[str, float]:
    with np.load(path) as data:
        labels = data["labels"].astype(np.int64)
        time = float(data["time"])
        step = int(data["step"])

    nlabels = int(labels.max()) + 1
    area = np.bincount(labels.ravel(), minlength=nlabels).astype(float)
    perimeter = np.zeros(nlabels, dtype=float)

    for neighbor in (np.roll(labels, -1, axis=0), np.roll(labels, -1, axis=1)):
        mask = labels != neighbor
        left = labels[mask]
        right = neighbor[mask]
        np.add.at(perimeter, left, 1.0)
        np.add.at(perimeter, right, 1.0)

    valid = area >= minimum_area
    valid &= perimeter > 0
    q = perimeter[valid] ** 2 / (4.0 * np.pi * area[valid])
    # The lattice estimator is biased above one; matched B0 is therefore the
    # reference. Relative changes in q diagnose extra boundary waviness.
    return {
        "time": time,
        "step": step,
        "grains_used": int(np.count_nonzero(valid)),
        "isoperimetric_q_mean": float(np.mean(q)) if len(q) else np.nan,
        "isoperimetric_q_median": float(np.median(q)) if len(q) else np.nan,
        "isoperimetric_q_p90": float(np.quantile(q, 0.9)) if len(q) else np.nan,
        "boundary_edges_total": float(perimeter.sum() / 2.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Quantify excess grain-boundary waviness from saved PF label frames."
    )
    parser.add_argument("video_root")
    parser.add_argument("--output", default=None)
    parser.add_argument("--minimum-area", type=int, default=25)
    args = parser.parse_args()

    root = Path(args.video_root)
    rows = []
    for run_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        frames = sorted((run_dir / "frames").glob("frame-*.npz"))
        for frame in frames:
            row = frame_metrics(frame, args.minimum_area)
            row["regime"] = run_dir.name.split("-T", 1)[0]
            row["run_dir"] = str(run_dir)
            rows.append(row)
    if not rows:
        raise SystemExit(f"no frame-*.npz files found below {root}")

    table = pd.DataFrame(rows).sort_values(["regime", "time"])
    output = Path(args.output) if args.output else root / "boundary_roughness.csv"
    table.to_csv(output, index=False)

    summary = table.groupby("regime", as_index=False).agg(
        q_mean=("isoperimetric_q_mean", "mean"),
        q_median=("isoperimetric_q_median", "median"),
        q_p90_mean=("isoperimetric_q_p90", "mean"),
        frames=("step", "count"),
    )
    reference = summary.loc[summary["regime"] == "B0_CTL", "q_mean"]
    if len(reference):
        q0 = float(reference.iloc[0])
        summary["q_mean_over_B0"] = summary["q_mean"] / q0
    print(summary.to_string(index=False))
    print(output)


if __name__ == "__main__":
    main()
