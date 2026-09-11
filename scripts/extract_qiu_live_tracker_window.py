#!/usr/bin/env python3
"""Extract a fixed, append-safe step window from live Qiu CSV trackers."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract(source: Path, destination: Path, first_step: int, last_step: int) -> dict[str, object]:
    if last_step < first_step:
        raise ValueError("last-step must be greater than or equal to first-step")
    source = source.resolve()
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    results: dict[str, object] = {}
    for name in ("grain_tracks.csv", "boundary_tracks.csv"):
        input_path = source / name
        output_path = destination / f"{name}.gz"
        row_count = 0
        steps: set[int] = set()
        with input_path.open("r", newline="", encoding="utf-8") as input_stream:
            reader = csv.DictReader(input_stream)
            if not reader.fieldnames or "step" not in reader.fieldnames:
                raise RuntimeError(f"{input_path} has no step column")
            with gzip.open(output_path, "wt", newline="", encoding="utf-8") as output_stream:
                writer = csv.DictWriter(output_stream, fieldnames=reader.fieldnames)
                writer.writeheader()
                for row in reader:
                    try:
                        step = int(row["step"])
                    except (KeyError, TypeError, ValueError):
                        # A concurrently appended final line may not yet be complete.
                        continue
                    if first_step <= step <= last_step:
                        writer.writerow(row)
                        row_count += 1
                        steps.add(step)
        results[name] = {
            "path": output_path.name,
            "sha256": sha256(output_path),
            "rows": row_count,
            "step_first": min(steps) if steps else None,
            "step_last": max(steps) if steps else None,
            "unique_step_count": len(steps),
        }
    manifest = {
        "schema_version": 1,
        "scientific_status": "live_bounded_tracker_window",
        "source": str(source),
        "requested_step_first": first_step,
        "requested_step_last": last_step,
        "files": results,
        "restart_capable": False,
        "interpretation": (
            "Rows are selected by their completed CSV step value from append-only live "
            "trackers. This window supplements, but does not replace, an atomic checkpoint."
        ),
    }
    (destination / "tracker_window_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("source", type=Path)
    result.add_argument("destination", type=Path)
    result.add_argument("--first-step", required=True, type=int)
    result.add_argument("--last-step", required=True, type=int)
    return result


if __name__ == "__main__":
    args = parser().parse_args()
    print(json.dumps(extract(args.source, args.destination, args.first_step, args.last_step), indent=2))
