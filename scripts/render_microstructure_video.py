#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def _color_table(max_label: int) -> np.ndarray:
    rng = np.random.default_rng(12345)
    colors = rng.random((max(max_label + 1, 256), 3))
    colors[0] = (0.1, 0.1, 0.1)
    return colors


def _global_limits(frames: list[Path], key: str) -> tuple[float, float]:
    lo = np.inf
    hi = -np.inf
    for frame in frames:
        with np.load(frame) as data:
            a = np.asarray(data[key], dtype=float)
            finite = a[np.isfinite(a)]
            if finite.size:
                lo = min(lo, float(finite.min()))
                hi = max(hi, float(finite.max()))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return -1.0, 1.0
    if lo == hi:
        span = max(abs(lo), 1.0)
        return lo - 0.5 * span, hi + 0.5 * span
    return lo, hi


def _load(frame: Path):
    with np.load(frame) as data:
        return (
            data["labels"].copy(),
            data["blocked"].copy(),
            data["shear"].copy(),
            data["free_volume"].copy(),
            int(data["step"]),
            float(data["time"]),
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render saved PF NPZ frames to a PNG sequence plus MP4/GIF."
    )
    parser.add_argument("run_dir")
    parser.add_argument("--output")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument(
        "--png-frames", action=argparse.BooleanOptionalAction, default=True,
        help="Write a rendered PNG image sequence (default: yes).",
    )
    parser.add_argument(
        "--composite", action=argparse.BooleanOptionalAction, default=True,
        help="Render microstructure, shear, and free-volume panels (default: yes).",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    frames = sorted((run_dir / "frames").glob("frame-*.npz"))
    if not frames:
        raise SystemExit(f"no saved frames under {run_dir / 'frames'}")

    max_label = 0
    for frame in frames:
        with np.load(frame) as data:
            max_label = max(max_label, int(data["labels"].max()))
    colors = _color_table(max_label)

    shear_lo, shear_hi = _global_limits(frames, "shear")
    shear_abs = max(abs(shear_lo), abs(shear_hi), np.finfo(float).eps)
    fv_lo, fv_hi = _global_limits(frames, "free_volume")
    if fv_lo == fv_hi:
        fv_hi = fv_lo + 1.0

    labels0, blocked0, shear0, fv0, step0, time0 = _load(frames[0])

    if args.composite:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
        ax_micro, ax_shear, ax_fv = axes
    else:
        fig, ax_micro = plt.subplots(figsize=(6, 6), constrained_layout=True)
        ax_shear = ax_fv = None

    image = ax_micro.imshow(colors[labels0 % len(colors)], interpolation="nearest", origin="lower")
    overlay = ax_micro.imshow(
        np.ma.masked_where(blocked0 == 0, blocked0),
        interpolation="nearest", origin="lower", alpha=0.65, cmap="Reds", vmin=0, vmax=1,
    )
    ax_micro.set_title("grain structure; blocked GB/TJ in red")
    ax_micro.set_xticks([])
    ax_micro.set_yticks([])

    shear_image = fv_image = None
    if args.composite:
        shear_image = ax_shear.imshow(
            shear0, interpolation="nearest", origin="lower", cmap="coolwarm",
            vmin=-shear_abs, vmax=shear_abs,
        )
        ax_shear.set_title("stored shear state")
        ax_shear.set_xticks([])
        ax_shear.set_yticks([])
        fig.colorbar(shear_image, ax=ax_shear, fraction=0.046, pad=0.04)

        fv_image = ax_fv.imshow(
            fv0, interpolation="nearest", origin="lower", cmap="viridis",
            vmin=fv_lo, vmax=fv_hi,
        )
        ax_fv.set_title("free-volume / climb deficit")
        ax_fv.set_xticks([])
        ax_fv.set_yticks([])
        fig.colorbar(fv_image, ax=ax_fv, fraction=0.046, pad=0.04)

    title = fig.suptitle(f"{run_dir.name}   step={step0}   t={time0:.3f}")

    def update(index: int):
        labels, blocked, shear, free_volume, step, time = _load(frames[index])
        image.set_data(colors[labels % len(colors)])
        overlay.set_data(np.ma.masked_where(blocked == 0, blocked))
        artists = [image, overlay, title]
        if args.composite:
            assert shear_image is not None and fv_image is not None
            shear_image.set_data(shear)
            fv_image.set_data(free_volume)
            artists.extend([shear_image, fv_image])
        title.set_text(f"{run_dir.name}   step={step}   t={time:.3f}")
        return artists

    if args.png_frames:
        png_dir = run_dir / "png_frames"
        png_dir.mkdir(parents=True, exist_ok=True)
        for index, source in enumerate(frames):
            update(index)
            destination = png_dir / f"frame-{index:05d}.png"
            if not destination.exists():
                fig.savefig(destination, dpi=args.dpi)
        print(f"PNG sequence: {png_dir} ({len(frames)} frames)")

    movie = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=1000 / args.fps, blit=False
    )
    requested = Path(args.output) if args.output else run_dir / "microstructure.mp4"
    requested.parent.mkdir(parents=True, exist_ok=True)

    if requested.suffix.lower() == ".mp4" and animation.writers.is_available("ffmpeg"):
        writer = animation.FFMpegWriter(fps=args.fps, bitrate=4500)
        movie.save(requested, writer=writer, dpi=args.dpi)
        output = requested
    else:
        output = requested.with_suffix(".gif")
        writer = animation.PillowWriter(fps=args.fps)
        movie.save(output, writer=writer, dpi=args.dpi)
    plt.close(fig)
    print(f"Movie: {output}")


if __name__ == "__main__":
    main()
