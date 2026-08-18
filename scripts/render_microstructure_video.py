#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description="Render saved PF label/blocked-state frames to MP4 or GIF.")
    parser.add_argument("run_dir")
    parser.add_argument("--output")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=140)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    frames = sorted((run_dir / "frames").glob("frame-*.npz"))
    if not frames:
        raise SystemExit(f"no saved frames under {run_dir / 'frames'}")

    first = np.load(frames[0])
    labels0 = first["labels"]
    blocked0 = first["blocked"]
    nlabels = int(max(labels0.max(), 1)) + 1
    rng = np.random.default_rng(12345)
    colors = rng.random((max(nlabels, 256), 3))
    colors[0] = (0.1, 0.1, 0.1)

    fig, ax = plt.subplots(figsize=(6, 6))
    image = ax.imshow(colors[labels0 % len(colors)], interpolation="nearest", origin="lower")
    overlay = ax.imshow(
        np.ma.masked_where(blocked0 == 0, blocked0),
        interpolation="nearest", origin="lower", alpha=0.65, cmap="Reds", vmin=0, vmax=1,
    )
    title = ax.set_title("")
    ax.set_xticks([])
    ax.set_yticks([])

    def update(index: int):
        with np.load(frames[index]) as data:
            labels = data["labels"]
            blocked = data["blocked"]
            step = int(data["step"])
            time = float(data["time"])
        image.set_data(colors[labels % len(colors)])
        overlay.set_data(np.ma.masked_where(blocked == 0, blocked))
        title.set_text(f"{run_dir.name}   step={step}   t={time:.3f}")
        return image, overlay, title

    movie = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000 / args.fps, blit=False)
    requested = Path(args.output) if args.output else run_dir / "microstructure.mp4"
    requested.parent.mkdir(parents=True, exist_ok=True)

    if requested.suffix.lower() == ".mp4" and animation.writers.is_available("ffmpeg"):
        writer = animation.FFMpegWriter(fps=args.fps, bitrate=3000)
        movie.save(requested, writer=writer, dpi=args.dpi)
        output = requested
    else:
        output = requested.with_suffix(".gif")
        writer = animation.PillowWriter(fps=args.fps)
        movie.save(output, writer=writer, dpi=args.dpi)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
