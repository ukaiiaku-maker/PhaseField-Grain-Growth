#!/usr/bin/env python3
"""Compare checksummed native and A0-port Qiu matched-state NPZ records."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

REQUIRED = (
    "accepted_step", "physical_time", "output_scaling", "orientations",
    "beta", "reference_directions", "line_density", "sigma11", "sigma12",
    "sigma22", "elastic_force", "barrier_eij", "delta_eta_pre",
    "delta_eta_accepted", "energy_components", "phi", "active_support",
    "grain_labels",
)
EXACT = {"accepted_step", "active_support", "grain_labels"}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def compare_state(native_path: Path, port_path: Path, rtol: float, atol: float) -> dict[str, object]:
    result: dict[str, object] = {
        "native_path": str(native_path), "port_path": str(port_path),
        "native_sha256": digest(native_path), "port_sha256": digest(port_path),
        "fields": {},
    }
    passed = True
    with np.load(native_path, allow_pickle=False) as native, np.load(port_path, allow_pickle=False) as port:
        missing = {"native": sorted(set(REQUIRED) - set(native.files)), "port": sorted(set(REQUIRED) - set(port.files))}
        if missing["native"] or missing["port"]:
            raise SystemExit(f"missing required matched-state fields: {missing}")
        result["accepted_step"] = int(np.asarray(native["accepted_step"]).reshape(()))
        result["physical_time"] = float(np.asarray(native["physical_time"]).reshape(()))
        for name in REQUIRED:
            left, right = native[name], port[name]
            same_shape = left.shape == right.shape
            finite = bool(np.all(np.isfinite(left)) and np.all(np.isfinite(right)))
            exact = name in EXACT
            difference = (
                np.not_equal(left, right).astype(float) if same_shape and exact
                else np.abs(left - right) if same_shape else np.array([np.inf])
            )
            scale = np.maximum(np.abs(left), np.abs(right)) if same_shape else np.array([1.0])
            threshold = atol + rtol * scale
            field_pass = bool(
                same_shape and finite
                and (np.array_equal(left, right) if exact else np.all(difference <= threshold))
            )
            result["fields"][name] = {
                "shape": list(left.shape), "max_abs": float(np.max(difference)),
                "max_scaled_error": float(np.max(difference / threshold)),
                "comparison": "exact" if exact else "tolerance", "pass": field_pass,
            }
            passed &= field_pass
    result["pass"] = passed
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("native", type=Path)
    parser.add_argument("port", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("--atol", type=float, default=1e-13)
    args = parser.parse_args()

    if args.native.is_dir() != args.port.is_dir():
        raise SystemExit("native and port inputs must both be files or both be directories")
    if args.native.is_dir():
        native_files = {path.name: path for path in args.native.glob("*.npz")}
        port_files = {path.name: path for path in args.port.glob("*.npz")}
        if not native_files or native_files.keys() != port_files.keys():
            raise SystemExit("matched-state directories must contain the same nonempty NPZ filename set")
        checkpoints = [
            compare_state(native_files[name], port_files[name], args.rtol, args.atol)
            for name in sorted(native_files)
        ]
        checkpoints.sort(key=lambda item: item["accepted_step"])
        failed = [item for item in checkpoints if not item["pass"]]
        report: dict[str, object] = {
            "model": "QIU_SI_REFERENCE_ANISO", "rtol": args.rtol, "atol": args.atol,
            "checkpoints": checkpoints, "checkpoint_count": len(checkpoints),
            "earliest_divergence": None if not failed else {
                "accepted_step": failed[0]["accepted_step"],
                "physical_time": failed[0]["physical_time"],
                "failed_fields": [name for name, value in failed[0]["fields"].items() if not value["pass"]],
            },
            "pass": not failed,
        }
    else:
        state = compare_state(args.native, args.port, args.rtol, args.atol)
        report = {"model": "QIU_SI_REFERENCE_ANISO", "rtol": args.rtol, "atol": args.atol, **state}
    passed = bool(report["pass"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not passed:
        raise SystemExit("QIU_SI_4REF_A0_PORT matched-state correspondence failed")


if __name__ == "__main__":
    main()
