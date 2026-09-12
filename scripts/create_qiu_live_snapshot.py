#!/usr/bin/env python3
"""Create a checkpoint-consistent safety archive from a live Qiu run.

The simulation keeps writing beyond its most recent atomic checkpoint.  This
utility uses the recorder counters and tracker byte offsets stored in that
checkpoint, rather than directory listings or filename arithmetic, to exclude
post-checkpoint output from the archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import tarfile
import tempfile
from pathlib import Path


PART_RE = re.compile(r"part-(\d+)\.parquet$")
STEP_RE = re.compile(r"(?:frame-|step-)(\d+)")


def _copy_prefix(source: Path, destination: Path, length: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    remaining = length
    with source.open("rb") as src, destination.open("xb") as dst:
        while remaining:
            chunk = src.read(min(8 * 1024 * 1024, remaining))
            if not chunk:
                raise RuntimeError(f"{source} ended before checkpoint offset {length}")
            dst.write(chunk)
            remaining -= len(chunk)


def _copy_parts(source: Path, destination: Path, closed_count: int) -> int:
    destination.mkdir(parents=True, exist_ok=True)
    copied = 0
    for path in sorted(source.glob("part-*.parquet")):
        match = PART_RE.search(path.name)
        if match and int(match.group(1)) < closed_count:
            shutil.copy2(path, destination / path.name)
            copied += 1
    if copied != closed_count:
        raise RuntimeError(
            f"expected {closed_count} closed parts in {source}, copied {copied}"
        )
    return copied


def _copy_step_files(
    source: Path,
    destination: Path,
    checkpoint_step: int,
    *,
    legacy_dense: bool = False,
) -> list[int]:
    destination.mkdir(parents=True, exist_ok=True)
    steps: list[int] = []
    if not source.exists():
        return steps
    for path in sorted(source.iterdir()):
        if not path.is_file():
            continue
        # macOS archive extraction can materialize AppleDouble sidecars such
        # as ``._frame-0001000.npz``. They are metadata, not restart/movie
        # fields, and their embedded filename must not make them look valid.
        if path.name.startswith("._"):
            continue
        match = STEP_RE.search(path.name)
        if not match:
            continue
        step = int(match.group(1))
        if step > checkpoint_step:
            continue
        if legacy_dense and not (step == 9001 or step % 10 == 0):
            continue
        shutil.copy2(path, destination / path.name)
        steps.append(step)
    return steps


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def create_snapshot(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    checkpoint_path = source / "checkpoint.json"
    checkpoint_before = checkpoint_path.read_bytes()
    checkpoint = json.loads(checkpoint_before)
    step = int(checkpoint["step_number"])
    if step != args.expected_step:
        raise RuntimeError(f"checkpoint is step {step}, expected {args.expected_step}")

    extension = checkpoint.get("extension_state", {})
    recorder = extension.get("qiu_forensics", {})
    scalar_parts = int(recorder["scalar_part"])
    boundary_parts = int(recorder.get("boundary_part", 0))
    event_state = extension.get("event_trace", {})
    event_state_parts = int(event_state.get("state_part_count", 0))
    event_parts = int(event_state.get("event_part_count", 0))

    destination_dir = args.destination.resolve()
    destination_dir.mkdir(parents=True, exist_ok=False)
    archive = destination_dir / args.archive_name
    staging = Path(tempfile.mkdtemp(prefix="qiu-live-snapshot-", dir=args.scratch))
    try:
        payload = staging / "qualification"
        target = payload if args.layout == "corrected" else payload / source.name
        target.mkdir(parents=True)

        (target / "checkpoint.json").write_bytes(checkpoint_before)
        shutil.copy2(source / "checkpoint.npz", target / "checkpoint.npz")
        shutil.copy2(source / "manifest.json", target / "manifest.json")
        _copy_prefix(
            source / "grain_tracks.csv",
            target / "grain_tracks.csv",
            int(checkpoint["grain_tracks_offset"]),
        )
        _copy_prefix(
            source / "boundary_tracks.csv",
            target / "boundary_tracks.csv",
            int(checkpoint["boundary_tracks_offset"]),
        )

        copied_scalar_parts = _copy_parts(
            source / "per_step_diagnostics.parquet",
            target / "per_step_diagnostics.parquet",
            scalar_parts,
        )
        copied_boundary_parts = _copy_parts(
            source / "per_boundary_diagnostics.parquet",
            target / "per_boundary_diagnostics.parquet",
            boundary_parts,
        )
        copied_event_state_parts = _copy_parts(
            source / "event_traces.parquet",
            target / "event_traces.parquet",
            event_state_parts,
        )
        copied_event_parts = _copy_parts(
            source / "event_trace_events.parquet",
            target / "event_trace_events.parquet",
            event_parts,
        )
        # This stream is not currently checkpoint-counted. It is empty for the
        # active jobs, so refuse ambiguity rather than copy live records.
        raw_events = source / "events.parquet"
        if raw_events.exists() and any(raw_events.iterdir()):
            raise RuntimeError("events.parquet is nonempty but has no checkpoint counter")
        (target / "events.parquet").mkdir(exist_ok=True)

        frame_steps = _copy_step_files(
            source / "frames", target / "frames", step
        )
        dense_steps = _copy_step_files(
            source / "diagnostic_fields",
            target / "diagnostic_fields",
            step,
            legacy_dense=args.layout == "legacy",
        )
        for name in ("diagnostic_capture.json", "diagnostic_capture.npz"):
            candidate = source / name
            if candidate.exists():
                shutil.copy2(candidate, target / name)

        checkpoint_after = checkpoint_path.read_bytes()
        if checkpoint_after != checkpoint_before:
            raise RuntimeError("checkpoint advanced while snapshot was assembled")

        temporary_archive = archive.with_suffix(archive.suffix + ".partial")
        with tarfile.open(temporary_archive, "w:gz", compresslevel=6) as stream:
            stream.add(payload, arcname="qualification", recursive=True)
        os.replace(temporary_archive, archive)
        result: dict[str, object] = {
            "archive": str(archive),
            "archive_sha256": _sha256(archive),
            "checkpoint_step": step,
            "checkpoint_time": checkpoint["time"],
            "grain_tracks_offset": checkpoint["grain_tracks_offset"],
            "boundary_tracks_offset": checkpoint["boundary_tracks_offset"],
            "closed_scalar_parts": copied_scalar_parts,
            "closed_boundary_parts": copied_boundary_parts,
            "closed_event_state_parts": copied_event_state_parts,
            "closed_event_parts": copied_event_parts,
            "last_grain_count": recorder.get("population_history", [[None, None]])[-1][1],
            "frame_steps": frame_steps,
            "dense_steps": dense_steps,
        }
        (destination_dir / "assembly.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        return result
    finally:
        shutil.rmtree(staging)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source", type=Path, required=True)
    result.add_argument("--destination", type=Path, required=True)
    result.add_argument("--archive-name", required=True)
    result.add_argument("--expected-step", type=int, required=True)
    result.add_argument("--layout", choices=("corrected", "legacy"), required=True)
    result.add_argument("--scratch", type=Path, required=True)
    return result


if __name__ == "__main__":
    print(json.dumps(create_snapshot(parser().parse_args()), indent=2))
