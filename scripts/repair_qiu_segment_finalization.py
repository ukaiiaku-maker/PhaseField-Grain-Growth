#!/usr/bin/env python3
"""Recover metadata-only finalization after a completed Qiu segment.

The fetched archive remains immutable. This writes a checksummed sidecar only
after independently validating its scientific checkpoint chain.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def member_bytes(stream: tarfile.TarFile, relative: str) -> bytes:
    for name in (f"./output/{relative}", f"output/{relative}"):
        try:
            member = stream.extractfile(name)
        except KeyError:
            member = None
        if member is not None:
            return member.read()
    raise KeyError(relative)


def repair(archive: Path, output: Path, kind: str, start: int, target: int,
           parent_sha256: str, control: str = "", energy: float = 0.0,
           mobility: float = 0.0) -> dict:
    output.mkdir(parents=True, exist_ok=False)
    with tarfile.open(archive, "r:gz") as stream:
        manifest = json.loads(member_bytes(stream, "checkpoints/checkpoint_manifest.json"))
        for item in manifest["checkpoints"]:
            payload = member_bytes(stream, f"checkpoints/{item['path']}")
            if digest(payload) != item["sha256"] or not item["finite"]:
                raise ValueError(f"invalid archived checkpoint: {item['path']}")
        restart_name = f"checkpoints/native-restart-step{target:06d}.npz"
        restart = member_bytes(stream, restart_name)
        with np.load(io.BytesIO(restart), allow_pickle=False) as state:
            if int(state["accepted_step"]) != target:
                raise ValueError("terminal restart checkpoint step mismatch")
        names = sorted(member.name for member in stream.getmembers() if member.isfile())
        compilation = None
        if kind == "control":
            compilation = json.loads(member_bytes(stream, "checkpoints/compilation_manifest.json"))

    if manifest["restart_used"] is not (start != 0):
        raise ValueError("restart lineage mismatch")
    if kind == "native":
        if manifest["baseline_eligible"] is not True or manifest["logical_trajectory"] != "QIU_SI_NATIVE_BASELINE_V2":
            raise ValueError("native manifest identity mismatch")
        classification = "QIU_SI_NATIVE_BASELINE_SEGMENT_PASSED"
        decision = {
            "schema": "qiu-native-production-segment-v1", "classification": classification,
            "logical_trajectory": "QIU_SI_NATIVE_BASELINE_V2", "segment_start": start,
            "segment_target": target, "parent_checkpoint_sha256": parent_sha256,
            "end_checkpoint": Path(restart_name).name, "end_checkpoint_sha256": digest(restart),
            "execution_threads": 1, "finite": True,
        }
    else:
        if manifest["baseline_eligible"] is not False or manifest["logical_trajectory"] != control:
            raise ValueError("control manifest identity mismatch")
        if compilation is None or compilation["object_mode_or_uncompiled"] != []:
            raise ValueError("control compilation is not nopython complete")
        classification = "QIU_SI_CONTROL_WINDOW_PASSED"
        decision = {
            "schema": "qiu-control-segment-v1", "classification": classification,
            "control": control, "segment_start": start, "segment_target": target,
            "parent_checkpoint_sha256": parent_sha256, "end_checkpoint_sha256": digest(restart),
            "energy_normalization": energy, "mobility_normalization": mobility,
            "execution_threads": 1, "finite": True,
        }
    decision_path = output / "segment_decision.json"
    decision_path.write_text(json.dumps(decision, indent=2) + "\n")
    record = {
        "schema": "qiu-segment-finalization-repair-v1", "classification": classification,
        "original_archive": str(archive), "original_archive_sha256": digest(archive.read_bytes()),
        "archived_file_count": len(names), "decision_sha256": digest(decision_path.read_bytes()),
        "repair_scope": "metadata-only; scientific checkpoints are unchanged",
    }
    (output / "repair_record.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--kind", choices=("native", "control"), required=True)
    parser.add_argument("--start", type=int, required=True); parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--parent-sha256", default=""); parser.add_argument("--control", default="")
    parser.add_argument("--energy", type=float, default=0.0); parser.add_argument("--mobility", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(repair(args.archive, args.output, args.kind, args.start, args.target,
                            args.parent_sha256, args.control, args.energy, args.mobility), indent=2))


if __name__ == "__main__":
    main()
