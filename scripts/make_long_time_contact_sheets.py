#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TARGETS = (1.00, 1.40, 1.85, 2.50, 3.16)


def _colors(max_label: int) -> np.ndarray:
    rng = np.random.default_rng(12345)
    values = rng.random((max(max_label + 1, 256), 3))
    values[0] = 0.1
    return values


def _frame_metrics(path: Path) -> tuple[float, int, float, np.ndarray]:
    with np.load(path) as data:
        return (
            float(data["G_over_G0"]) if "G_over_G0" in data else np.nan,
            int(data["grain_count"]), float(data["time"]), data["labels"].copy(),
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign")
    parser.add_argument("--output")
    args = parser.parse_args()
    root = Path(args.campaign)
    output = Path(args.output) if args.output else root / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    campaign = json.loads((root / "campaign_manifest.json").read_text())
    runs = [Path(item) for item in campaign["runs"]]

    selected: dict[Path, list[tuple[Path, tuple[float, int, float, np.ndarray]]]] = {}
    max_label = 0
    for run in runs:
        frames = [(path, _frame_metrics(path)) for path in sorted((run / "frames").glob("frame-*.npz"))]
        if not frames:
            continue
        choices = []
        for target in TARGETS:
            available = [item for item in frames if np.isfinite(item[1][0])]
            if not available or max(item[1][0] for item in available) + 1e-12 < target:
                continue
            choice = min(available, key=lambda item: abs(item[1][0] - target))
            choices.append(choice)
            max_label = max(max_label, int(choice[1][3].max()))
        selected[run] = choices
    colors = _colors(max_label)

    rows = len(selected)
    fig, axes = plt.subplots(rows, len(TARGETS), figsize=(3 * len(TARGETS), 3 * rows), squeeze=False, constrained_layout=True)
    for row, (run, choices) in enumerate(selected.items()):
        by_target = {min(TARGETS, key=lambda target: abs(item[1][0] - target)): item for item in choices}
        for column, target in enumerate(TARGETS):
            ax = axes[row, column]
            item = by_target.get(target)
            if item is None:
                ax.text(0.5, 0.5, "not reached", ha="center", va="center")
            else:
                _, (actual, grains, time, labels) = item
                ax.imshow(colors[labels % len(colors)], origin="lower", interpolation="nearest")
                ax.set_title(f"x={actual:.2f}, N={grains}, t={time:.1f}", fontsize=8)
            if column == 0:
                ax.set_ylabel(run.name, fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(output / "progress_matched_contact_sheet.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, rows, figsize=(3 * rows, 3), squeeze=False, constrained_layout=True)
    for column, run in enumerate(selected):
        path = sorted((run / "frames").glob("frame-*.npz"))[-1]
        actual, grains, time, labels = _frame_metrics(path)
        axes[0, column].imshow(colors[labels % len(colors)], origin="lower", interpolation="nearest")
        axes[0, column].set_title(f"{run.name}\nx={actual:.2f}, N={grains}", fontsize=7)
        axes[0, column].set_xticks([]); axes[0, column].set_yticks([])
    fig.savefig(output / "terminal_contact_sheet.png", dpi=160)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
