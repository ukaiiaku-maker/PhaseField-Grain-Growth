#!/usr/bin/env python3
"""Render restart-safe QIU/FFT-v2 field archives and write a frame index."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover(run: Path) -> list[Path]:
    frames = sorted((run / "frames").glob("frame-*.npz"))
    if frames:
        return frames
    fields = sorted((run / "diagnostic_fields").glob("step-*.npz"))
    if fields:
        return fields
    raise SystemExit(f"no movie-capable NPZ archives under {run}")


def load(path: Path) -> dict[str, object]:
    with np.load(path) as data:
        labels = np.asarray(data["labels"], dtype=np.int64)
        stress = np.asarray(data["stress"], dtype=float)
        eigenstrain = np.asarray(data["eigenstrain"], dtype=float)
        return {
            "labels": labels,
            "stress_xy": stress[0, 1],
            "stress_norm": np.sqrt(np.sum(stress * stress, axis=(0, 1))),
            "eigenstrain_xy": eigenstrain[0, 1],
            "step": int(data["step"]), "time": float(data["time"]),
            "grain_count": (
                int(data["grain_count"])
                if "grain_count" in data else int(len(np.unique(labels)))
            ),
        }


def color_table(max_label: int) -> np.ndarray:
    colors = np.random.default_rng(12345).random((max(max_label + 1, 256), 3))
    colors[0] = (0.08, 0.08, 0.08)
    return colors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=130)
    args = parser.parse_args()
    paths = discover(args.run)
    records = [load(path) for path in paths]
    destination = args.output or args.run / "qiu_fields.mp4"
    destination.parent.mkdir(parents=True, exist_ok=True)
    index_path = destination.with_suffix(".frames.csv")
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "frame", "source", "source_sha256", "step", "time", "grain_count"
        ))
        writer.writeheader()
        for index, (path, record) in enumerate(zip(paths, records)):
            writer.writerow({
                "frame": index, "source": str(path.resolve()),
                "source_sha256": sha256(path), "step": record["step"],
                "time": record["time"], "grain_count": record["grain_count"],
            })

    maximum_label = max(int(np.max(record["labels"])) for record in records)
    colors = color_table(maximum_label)
    stress_xy_limit = max(
        float(np.quantile(np.abs(record["stress_xy"]), 0.995)) for record in records
    ) or 1.0
    eigenstrain_limit = max(
        float(np.quantile(np.abs(record["eigenstrain_xy"]), 0.995)) for record in records
    ) or 1.0
    stress_norm_limit = max(
        float(np.quantile(record["stress_norm"], 0.995)) for record in records
    ) or 1.0

    first = records[0]
    figure, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    micro = axes[0, 0].imshow(colors[first["labels"] % len(colors)], origin="lower")
    shear = axes[0, 1].imshow(
        first["stress_xy"], origin="lower", cmap="coolwarm",
        vmin=-stress_xy_limit, vmax=stress_xy_limit,
    )
    eigen = axes[1, 0].imshow(
        first["eigenstrain_xy"], origin="lower", cmap="coolwarm",
        vmin=-eigenstrain_limit, vmax=eigenstrain_limit,
    )
    norm = axes[1, 1].imshow(
        first["stress_norm"], origin="lower", cmap="magma",
        vmin=0.0, vmax=stress_norm_limit,
    )
    titles = ("grain labels", "stress $\\sigma_{yx}$", "eigenstrain $\\epsilon^*_{yx}$", "stress Frobenius norm")
    for axis, title in zip(axes.flat, titles):
        axis.set_title(title); axis.set_xticks([]); axis.set_yticks([])
    figure.colorbar(shear, ax=axes[0, 1], fraction=0.046)
    figure.colorbar(eigen, ax=axes[1, 0], fraction=0.046)
    figure.colorbar(norm, ax=axes[1, 1], fraction=0.046)
    heading = figure.suptitle("")

    def update(index: int):
        record = records[index]
        micro.set_data(colors[record["labels"] % len(colors)])
        shear.set_data(record["stress_xy"])
        eigen.set_data(record["eigenstrain_xy"])
        norm.set_data(record["stress_norm"])
        heading.set_text(
            f"{args.run.name} | frame {index + 1}/{len(records)} | "
            f"step {record['step']} | t={record['time']:.6g} | N={record['grain_count']}"
        )
        return micro, shear, eigen, norm, heading

    update(0)
    movie = animation.FuncAnimation(
        figure, update, frames=len(records), interval=1000 / args.fps, blit=False
    )
    if destination.suffix.lower() == ".mp4" and animation.writers.is_available("ffmpeg"):
        movie.save(destination, writer=animation.FFMpegWriter(fps=args.fps, bitrate=5000), dpi=args.dpi)
    else:
        destination = destination.with_suffix(".gif")
        movie.save(destination, writer=animation.PillowWriter(fps=args.fps), dpi=args.dpi)

    selected = sorted({0, len(records) // 2, len(records) - 1})
    for panel, index in zip(axes.flat, (selected + [selected[-1]])[:4]):
        record = records[index]
        panel.clear(); panel.imshow(colors[record["labels"] % len(colors)], origin="lower")
        panel.set_title(f"step {record['step']}, t={record['time']:.4g}, N={record['grain_count']}")
        panel.set_xticks([]); panel.set_yticks([])
    contact = destination.with_suffix(".contact_sheet.png")
    figure.suptitle(f"{args.run.name}: initial / midpoint / terminal morphology")
    figure.savefig(contact, dpi=args.dpi)
    plt.close(figure)
    metadata = {
        "schema_version": 1, "run": str(args.run.resolve()),
        "movie": str(destination.resolve()), "frame_index": str(index_path.resolve()),
        "contact_sheet": str(contact.resolve()), "frame_count": len(records),
        "first_step": records[0]["step"], "last_step": records[-1]["step"],
    }
    destination.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
