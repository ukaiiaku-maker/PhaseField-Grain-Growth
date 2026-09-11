#!/usr/bin/env python3
"""Preserve immutable closed live diagnostics beyond the latest checkpoint.

This is a forensic fragment, not a restart archive.  It supplements a separately
verified checkpoint snapshot when a transition occurs between checkpoint
cadences.  Only explicitly bounded, already-closed Parquet parts and atomic NPZ
field files are copied.
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _copy_stable(source: Path, destination: Path) -> None:
    before = source.stat()
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f"source changed during copy: {source}")


def _copy_prefix(source: Path, destination: Path, length: int) -> None:
    before = source.stat()
    destination.parent.mkdir(parents=True, exist_ok=True)
    remaining = length
    with source.open("rb") as src, destination.open("xb") as dst:
        while remaining:
            block = src.read(min(8 * 1024 * 1024, remaining))
            if not block:
                raise RuntimeError(f"{source} ended before checkpoint offset {length}")
            dst.write(block)
            remaining -= len(block)
    after = source.stat()
    if after.st_size < length or before.st_size > after.st_size:
        raise RuntimeError(f"tracker became inconsistent during prefix copy: {source}")


def create_fragment(args: argparse.Namespace) -> dict[str, object]:
    source = args.source.resolve()
    checkpoint_bytes = (source / "checkpoint.json").read_bytes()
    checkpoint = json.loads(checkpoint_bytes)
    checkpoint_step = int(checkpoint["step_number"])
    if checkpoint_step != args.checkpoint_step:
        raise RuntimeError(
            f"checkpoint is step {checkpoint_step}, expected {args.checkpoint_step}"
        )
    if args.through_step <= checkpoint_step:
        raise ValueError("through-step must be later than the checkpoint")

    parts = source / "per_step_diagnostics.parquet"
    selected_parts: list[Path] = []
    for path in sorted(parts.glob("part-*.parquet")):
        match = PART_RE.fullmatch(path.name)
        if match and int(match.group(1)) < args.closed_part_count:
            selected_parts.append(path)
    if len(selected_parts) != args.closed_part_count:
        raise RuntimeError(
            f"expected {args.closed_part_count} closed scalar parts, found "
            f"{len(selected_parts)}"
        )

    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    archive = destination / args.archive_name
    staging = Path(tempfile.mkdtemp(prefix="qiu-transition-fragment-", dir=args.scratch))
    try:
        target = staging / "qualification" / source.name
        target.mkdir(parents=True)
        (target / "checkpoint.json").write_bytes(checkpoint_bytes)
        _copy_stable(source / "checkpoint.npz", target / "checkpoint.npz")
        _copy_stable(source / "manifest.json", target / "manifest.json")
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
        for part in selected_parts:
            _copy_stable(part, target / "per_step_diagnostics.parquet" / part.name)

        field_steps: list[int] = []
        field_files: list[str] = []
        for path in sorted((source / "diagnostic_fields").glob("*.npz")):
            match = STEP_RE.search(path.name)
            if not match:
                continue
            step = int(match.group(1))
            retain = (
                step == args.diagnostic_start
                or step % args.coarse_field_stride == 0
                or step >= args.dense_field_start
            )
            if step <= args.through_step and retain:
                _copy_stable(path, target / "diagnostic_fields" / path.name)
                field_steps.append(step)
                field_files.append(path.name)

        frame_steps: list[int] = []
        for path in sorted((source / "frames").glob("*.npz")):
            match = STEP_RE.search(path.name)
            if match and int(match.group(1)) <= args.through_step:
                _copy_stable(path, target / "frames" / path.name)
                frame_steps.append(int(match.group(1)))
        capture = source / "diagnostic_capture.json"
        if capture.exists():
            _copy_stable(capture, target / capture.name)

        if (source / "checkpoint.json").read_bytes() != checkpoint_bytes:
            raise RuntimeError("checkpoint advanced during fragment assembly")

        temporary = archive.with_suffix(archive.suffix + ".partial")
        with tarfile.open(temporary, "w:gz", compresslevel=6) as output:
            output.add(staging / "qualification", arcname="qualification")
        os.replace(temporary, archive)
        result: dict[str, object] = {
            "archive": str(archive),
            "archive_sha256": _sha256(archive),
            "restart_capable_through_step": checkpoint_step,
            "diagnostics_complete_through_step": args.through_step,
            "closed_scalar_parts": len(selected_parts),
            "field_file_count": len(field_files),
            "field_step_first": min(field_steps) if field_steps else None,
            "field_step_last": max(field_steps) if field_steps else None,
            "field_unique_step_count": len(set(field_steps)),
            "frame_steps": frame_steps,
            "selection": {
                "diagnostic_start": args.diagnostic_start,
                "coarse_field_stride": args.coarse_field_stride,
                "dense_field_start": args.dense_field_start,
            },
        }
        (destination / "assembly.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        return result
    finally:
        shutil.rmtree(staging)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source", required=True, type=Path)
    result.add_argument("--destination", required=True, type=Path)
    result.add_argument("--archive-name", required=True)
    result.add_argument("--checkpoint-step", required=True, type=int)
    result.add_argument("--through-step", required=True, type=int)
    result.add_argument("--closed-part-count", required=True, type=int)
    result.add_argument("--diagnostic-start", type=int, default=9001)
    result.add_argument("--coarse-field-stride", type=int, default=10)
    result.add_argument("--dense-field-start", required=True, type=int)
    result.add_argument("--scratch", required=True, type=Path)
    return result


if __name__ == "__main__":
    print(json.dumps(create_fragment(parser().parse_args()), indent=2))
