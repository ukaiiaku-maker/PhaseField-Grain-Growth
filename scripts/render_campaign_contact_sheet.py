#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _color_table(max_label: int) -> np.ndarray:
    rng = np.random.default_rng(12345)
    colors = rng.random((max(max_label + 1, 256), 3))
    colors[0] = (0.1, 0.1, 0.1)
    return colors


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render the terminal saved frame from every campaign run."
    )
    parser.add_argument("campaign_root")
    parser.add_argument("--output", required=True)
    parser.add_argument("--columns", type=int, default=3)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    root = Path(args.campaign_root)
    runs = []
    for run in sorted(path for path in root.iterdir() if path.is_dir()):
        frames = sorted((run / "frames").glob("frame-*.npz"))
        if frames:
            runs.append((run, frames[-1]))
    if not runs:
        raise SystemExit(f"no frame-*.npz files found below {root}")

    max_label = 0
    for _, frame in runs:
        with np.load(frame) as data:
            max_label = max(max_label, int(data["labels"].max()))
    colors = _color_table(max_label)

    columns = max(1, min(args.columns, len(runs)))
    rows = math.ceil(len(runs) / columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(4.4 * columns, 4.8 * rows), squeeze=False,
        constrained_layout=True,
    )
    for axis, (run, frame) in zip(axes.flat, runs, strict=False):
        with np.load(frame) as data:
            labels = data["labels"].astype(np.int64)
            step = int(data["step"])
            grain_count = int(
                data["grain_count"] if "grain_count" in data else np.unique(labels).size
            )
            stiffness = float(
                data["shear_stiffness"] if "shear_stiffness" in data else np.nan
            )
            seed = int(data["seed"] if "seed" in data else -1)
        axis.imshow(colors[labels % len(colors)], interpolation="nearest", origin="lower")
        axis.set_title(
            f"{run.name.split('-T', 1)[0]}\n"
            f"seed={seed}, step={step}, N={grain_count}, $K_s$={stiffness:g}",
            fontsize=10,
        )
        axis.set_xticks([])
        axis.set_yticks([])
    for axis in axes.flat[len(runs):]:
        axis.set_visible(False)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
