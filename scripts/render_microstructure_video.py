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
            if key not in data:
                continue
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
        labels = data["labels"].copy()
        blocked = data["blocked"].copy()
        shear = data["shear"].copy()
        free_volume = data["free_volume"].copy()
        mobility = (
            data["mobility"].copy()
            if "mobility" in data
            else np.ones(labels.shape, dtype=np.float32)
        )
        pending = data["pending_state"].copy() if "pending_state" in data else blocked.copy()
        climb_stage = (
            data["climb_stage"].copy()
            if "climb_stage" in data else np.zeros(labels.shape, dtype=np.uint8)
        )
        sink_activity = (
            data["sink_activity"].copy()
            if "sink_activity" in data else np.zeros(labels.shape, dtype=np.uint8)
        )
        return (
            labels, blocked, shear, free_volume, mobility, pending, climb_stage,
            sink_activity, int(data["step"]), float(data["time"]),
            float(data["temperature"]) if "temperature" in data else np.nan,
            int(data["seed"]) if "seed" in data else -1,
            float(data["shear_stiffness"]) if "shear_stiffness" in data else np.nan,
            int(data["grain_count"]) if "grain_count" in data else int(len(np.unique(labels))),
            float(data["N_required"]) if "N_required" in data else float(np.sum(free_volume)),
            float(data["N_accommodated_GB"]) if "N_accommodated_GB" in data else 0.0,
            float(data["N_accommodated_TJ"]) if "N_accommodated_TJ" in data else 0.0,
            float(data["conservation_residual"]) if "conservation_residual" in data else 0.0,
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
        help=(
            "Render grain structure, mobility/pinning footprint, shear, and "
            "free-volume panels (default: yes)."
        ),
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

    (
        labels0, blocked0, shear0, fv0, mobility0, pending0, stage0,
        activity0, step0, time0, temperature0, seed0, stiffness0, grains0,
        required0, gb0, tj0, residual0,
    ) = _load(frames[0])

    if args.composite:
        fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
        ax_micro, ax_pending, ax_mobility, ax_shear, ax_fv, ax_sink = axes.flat
    else:
        fig, ax_micro = plt.subplots(figsize=(6, 6), constrained_layout=True)
        ax_pending = ax_mobility = ax_shear = ax_fv = ax_sink = None

    image = ax_micro.imshow(
        colors[labels0 % len(colors)], interpolation="nearest", origin="lower"
    )
    overlay = ax_micro.imshow(
        np.ma.masked_where(blocked0 == 0, blocked0),
        interpolation="nearest", origin="lower", alpha=0.65,
        cmap="Reds", vmin=0, vmax=1,
    )
    ax_micro.set_title("grain structure; blocked GB domains in red")
    ax_micro.set_xticks([])
    ax_micro.set_yticks([])

    pending_image = mobility_image = shear_image = fv_image = None
    stage_image = activity_overlay = None
    if args.composite:
        pending_image = ax_pending.imshow(
            pending0, interpolation="nearest", origin="lower", cmap="tab10",
            vmin=0, vmax=7,
        )
        ax_pending.set_title("pending bits: G=1, T=2, C=4")
        ax_pending.set_xticks([])
        ax_pending.set_yticks([])

        mobility_image = ax_mobility.imshow(
            mobility0, interpolation="nearest", origin="lower", cmap="gray",
            vmin=0.0, vmax=1.0,
        )
        ax_mobility.set_title("mobility scale: black = pinned")
        ax_mobility.set_xticks([])
        ax_mobility.set_yticks([])
        fig.colorbar(mobility_image, ax=ax_mobility, fraction=0.046, pad=0.04)

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

        stage_image = ax_sink.imshow(
            stage0, interpolation="nearest", origin="lower", cmap="viridis",
            vmin=0, vmax=4,
        )
        activity_overlay = ax_sink.imshow(
            np.ma.masked_where(activity0 == 0, activity0),
            interpolation="nearest", origin="lower", cmap="autumn",
            vmin=0, vmax=2, alpha=0.9,
        )
        ax_sink.set_title("climb stage 0–4; sink GB=1, TJ=2")
        ax_sink.set_xticks([])
        ax_sink.set_yticks([])

    title = fig.suptitle(
        f"{run_dir.name}  T={temperature0:g} K  seed={seed0}  Ks={stiffness0:g}  "
        f"step={step0}  t={time0:.3f}  N={grains0}  "
        f"defects req={required0:.3g} GB={gb0:.3g} TJ={tj0:.3g} eps={residual0:.1e}"
    )

    def update(index: int):
        (
            labels, blocked, shear, free_volume, mobility, pending, stage,
            activity, step, time, temperature, seed, stiffness, grains, required,
            gb_sink, tj_sink, residual,
        ) = _load(frames[index])
        image.set_data(colors[labels % len(colors)])
        overlay.set_data(np.ma.masked_where(blocked == 0, blocked))
        artists = [image, overlay, title]
        if args.composite:
            assert mobility_image is not None
            assert pending_image is not None and shear_image is not None and fv_image is not None
            assert stage_image is not None and activity_overlay is not None
            pending_image.set_data(pending)
            mobility_image.set_data(mobility)
            shear_image.set_data(shear)
            fv_image.set_data(free_volume)
            stage_image.set_data(stage)
            activity_overlay.set_data(np.ma.masked_where(activity == 0, activity))
            artists.extend([
                pending_image, mobility_image, shear_image, fv_image,
                stage_image, activity_overlay,
            ])
        title.set_text(
            f"{run_dir.name}  T={temperature:g} K  seed={seed}  Ks={stiffness:g}  "
            f"step={step}  t={time:.3f}  N={grains}  "
            f"defects req={required:.3g} GB={gb_sink:.3g} TJ={tj_sink:.3g} eps={residual:.1e}"
        )
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
