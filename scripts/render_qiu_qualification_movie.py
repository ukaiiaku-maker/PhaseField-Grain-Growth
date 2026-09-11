#!/usr/bin/env python3
"""Render restart-safe QIU/FFT-v2 field archives and write a frame index."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib.animation as animation
from matplotlib.colors import SymLogNorm
import matplotlib.pyplot as plt
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover(run: Path) -> list[Path]:
    # Dense diagnostic fields intentionally supersede compact cadence frames at
    # duplicate steps, while both sources are retained at unique steps.
    by_step: dict[int, Path] = {}
    for path in sorted((run / "frames").glob("frame-*.npz")):
        match = re.search(r"frame-(\d+)", path.name)
        if match:
            by_step[int(match.group(1))] = path
    for path in sorted((run / "diagnostic_fields").glob("step-*.npz")):
        match = re.search(r"step-(\d+)", path.name)
        if match:
            by_step[int(match.group(1))] = path
    if not by_step:
        raise SystemExit(f"no movie-capable NPZ archives under {run}")
    return [by_step[step] for step in sorted(by_step)]


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
    destination = args.output or args.run / "qiu_fields.mp4"
    destination.parent.mkdir(parents=True, exist_ok=True)
    index_path = destination.with_suffix(".frames.csv")
    metadata_records: list[dict[str, object]] = []
    maximum_label = 0
    stress_quantiles: list[float] = []
    eigenstrain_quantiles: list[float] = []
    stress_norm_quantiles: list[float] = []
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "frame", "source", "source_kind", "source_sha256", "step", "time", "grain_count"
        ))
        writer.writeheader()
        for index, path in enumerate(paths):
            record = load(path)
            maximum_label = max(maximum_label, int(np.max(record["labels"])))
            stress_quantiles.append(float(np.quantile(np.abs(record["stress_xy"]), 0.995)))
            eigenstrain_quantiles.append(float(np.quantile(np.abs(record["eigenstrain_xy"]), 0.995)))
            stress_norm_quantiles.append(float(np.quantile(record["stress_norm"], 0.995)))
            row = {
                "frame": index, "source": str(path.resolve()),
                "source_kind": path.parent.name,
                "source_sha256": sha256(path), "step": record["step"],
                "time": record["time"], "grain_count": record["grain_count"],
            }
            metadata_records.append(row)
            writer.writerow(row)

    palette = color_table(maximum_label)
    stress_xy_limit = max(stress_quantiles) or 1.0
    eigenstrain_limit = max(eigenstrain_quantiles) or 1.0
    stress_norm_limit = max(stress_norm_quantiles) or 1.0
    nonzero_stress = [value for value in stress_quantiles if value > 0.0]
    nonzero_eigenstrain = [value for value in eigenstrain_quantiles if value > 0.0]
    typical_stress = float(np.median(nonzero_stress)) if nonzero_stress else stress_xy_limit
    typical_eigenstrain = (
        float(np.median(nonzero_eigenstrain)) if nonzero_eigenstrain else eigenstrain_limit
    )
    stress_scale = SymLogNorm(
        linthresh=max(float(typical_stress) * 0.1, stress_xy_limit * 1e-12),
        vmin=-stress_xy_limit, vmax=stress_xy_limit, base=10,
    )
    eigenstrain_scale = SymLogNorm(
        linthresh=max(float(typical_eigenstrain) * 0.1, eigenstrain_limit * 1e-12),
        vmin=-eigenstrain_limit, vmax=eigenstrain_limit, base=10,
    )

    first = load(paths[0])
    figure, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    micro = axes[0, 0].imshow(palette[first["labels"] % len(palette)], origin="lower")
    shear = axes[0, 1].imshow(
        first["stress_xy"], origin="lower", cmap="coolwarm", norm=stress_scale,
    )
    eigen = axes[1, 0].imshow(
        first["eigenstrain_xy"], origin="lower", cmap="coolwarm", norm=eigenstrain_scale,
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
        record = load(paths[index])
        micro.set_data(palette[record["labels"] % len(palette)])
        shear.set_data(record["stress_xy"])
        eigen.set_data(record["eigenstrain_xy"])
        norm.set_data(record["stress_norm"])
        heading.set_text(
            f"{args.run.name} | frame {index + 1}/{len(paths)} | "
            f"step {record['step']} | t={record['time']:.6g} | N={record['grain_count']}"
        )
        return micro, shear, eigen, norm, heading

    update(0)
    movie = animation.FuncAnimation(
        figure, update, frames=len(paths), interval=1000 / args.fps,
        blit=False, cache_frame_data=False,
    )
    if destination.suffix.lower() == ".mp4" and animation.writers.is_available("ffmpeg"):
        movie.save(destination, writer=animation.FFMpegWriter(fps=args.fps, bitrate=5000), dpi=args.dpi)
    else:
        destination = destination.with_suffix(".gif")
        movie.save(destination, writer=animation.PillowWriter(fps=args.fps), dpi=args.dpi)

    losses = np.asarray([
        int(metadata_records[index - 1]["grain_count"]) - int(metadata_records[index]["grain_count"])
        for index in range(1, len(metadata_records))
    ])
    initial_grains = int(metadata_records[0]["grain_count"])
    if len(losses) and int(np.max(losses)) >= max(10, int(np.ceil(0.05 * initial_grains))):
        transition = int(np.argmax(losses)) + 1
        selected = [max(0, transition - 1), transition, min(len(paths) - 1, transition + 1)]
        roles = ["pre-transition", "transition", "post-transition"]
    else:
        selected = [0, len(paths) // 2, len(paths) - 1]
        roles = ["initial", "midpoint", "terminal"]
    contact_figure, contact_axes = plt.subplots(3, 3, figsize=(12, 11), constrained_layout=True)
    contact_stress = None
    contact_eigenstrain = None
    contact_panels = []
    for column, (role, index) in enumerate(zip(roles, selected)):
        record = load(paths[index])
        contact_axes[0, column].imshow(palette[record["labels"] % len(palette)], origin="lower")
        contact_stress = contact_axes[1, column].imshow(
            record["stress_xy"], origin="lower", cmap="coolwarm", norm=stress_scale,
        )
        contact_eigenstrain = contact_axes[2, column].imshow(
            record["eigenstrain_xy"], origin="lower", cmap="coolwarm", norm=eigenstrain_scale,
        )
        contact_axes[0, column].set_title(
            f"{role}\nstep {record['step']}, t={record['time']:.4g}, N={record['grain_count']}"
        )
        contact_panels.append({
            "role": role, "frame": index, "step": record["step"],
            "time": record["time"], "grain_count": record["grain_count"],
            "source": str(paths[index].resolve()),
        })
    for row, label in enumerate(("grain labels", "stress $\\sigma_{yx}$", "eigenstrain $\\epsilon^*_{yx}$")):
        contact_axes[row, 0].set_ylabel(label)
        for axis in contact_axes[row]:
            axis.set_xticks([]); axis.set_yticks([])
    contact_figure.colorbar(contact_stress, ax=contact_axes[1, :], fraction=0.02)
    contact_figure.colorbar(contact_eigenstrain, ax=contact_axes[2, :], fraction=0.02)
    contact = destination.with_suffix(".contact_sheet.png")
    contact_figure.suptitle(f"{args.run.name}: transition-resolved fields")
    contact_figure.savefig(contact, dpi=args.dpi)
    plt.close(contact_figure)
    plt.close(figure)
    metadata = {
        "schema_version": 1, "run": str(args.run.resolve()),
        "movie": str(destination.resolve()), "frame_index": str(index_path.resolve()),
        "contact_sheet": str(contact.resolve()), "contact_sheet_panels": contact_panels,
        "frame_count": len(paths),
        "compact_frame_count": sum(path.parent.name == "frames" for path in paths),
        "dense_field_count": sum(path.parent.name == "diagnostic_fields" for path in paths),
        "first_step": metadata_records[0]["step"],
        "last_step": metadata_records[-1]["step"],
    }
    destination.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
