#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render every saved PF run in a video campaign to PNG frames and MP4/GIF."
    )
    parser.add_argument("video_root")
    parser.add_argument("--fps", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument("--no-png-frames", action="store_true")
    parser.add_argument("--microstructure-only", action="store_true")
    args = parser.parse_args()

    root = Path(args.video_root)
    if not root.exists():
        raise SystemExit(f"video campaign does not exist: {root}")

    runs = sorted(path for path in root.iterdir() if (path / "frames").is_dir())
    if not runs:
        raise SystemExit(f"no run directories containing frames/ under {root}")

    failures = []
    for run in runs:
        cmd = [
            sys.executable,
            str(Path(__file__).with_name("render_microstructure_video.py")),
            str(run),
            "--fps", str(args.fps),
            "--dpi", str(args.dpi),
        ]
        if args.no_png_frames:
            cmd.append("--no-png-frames")
        if args.microstructure_only:
            cmd.append("--no-composite")
        print("\n===", run.name, "===", flush=True)
        completed = subprocess.run(cmd)
        if completed.returncode != 0:
            failures.append(run.name)

    print(f"\nRendered {len(runs) - len(failures)}/{len(runs)} runs.")
    if failures:
        print("Failed:", ", ".join(failures))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
