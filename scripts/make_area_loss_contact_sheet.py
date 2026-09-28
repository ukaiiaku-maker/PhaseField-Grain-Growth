#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _runs(root: Path) -> list[Path]:
    name = "campaign_manifest.json" if (root / "campaign_manifest.json").exists() else "video_manifest.json"
    manifest = json.loads((root / name).read_text())
    return [Path(item["path"] if isinstance(item, dict) else item) for item in manifest["runs"]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output", required=True)
    parser.add_argument("--columns", type=int, default=4)
    args = parser.parse_args()
    root = Path(args.campaign)
    runs = _runs(root)
    columns = max(1, args.columns)
    rows = math.ceil(len(runs) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(4 * columns, 4 * rows), squeeze=False)
    rng = np.random.default_rng(12345)
    colors = rng.random((65536, 3))
    colors[0] = 0.1
    for ax, run in zip(axes.flat, runs, strict=False):
        frames = sorted((run / "frames").glob("frame-*.npz"))
        if not frames:
            ax.text(0.5, 0.5, "no frames", ha="center", va="center")
            ax.set_title(run.name)
            ax.axis("off")
            continue
        with np.load(frames[-1]) as data:
            labels = data["labels"]
            pending = data["pending_state"] if "pending_state" in data else data["blocked"]
            step = int(data["step"])
            time = float(data["time"])
            grains = int(data["grain_count"]) if "grain_count" in data else len(np.unique(labels))
        ax.imshow(
            colors[labels.astype(np.int64) % len(colors)],
            origin="lower", interpolation="nearest",
        )
        ax.imshow(
            np.ma.masked_where(pending == 0, pending), origin="lower",
            interpolation="nearest", cmap="autumn", alpha=0.65, vmin=0, vmax=7,
        )
        ax.set_title(f"{run.name}\nstep={step}, t={time:.2f}, N={grains}", fontsize=9)
        ax.axis("off")
    for ax in axes.flat[len(runs):]:
        ax.axis("off")
    fig.suptitle(f"Final long-run states: {root.name}")
    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(destination)


if __name__ == "__main__":
    main()
