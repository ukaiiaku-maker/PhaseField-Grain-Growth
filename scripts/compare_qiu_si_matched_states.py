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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("native", type=Path)
    parser.add_argument("port", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("--atol", type=float, default=1e-13)
    args = parser.parse_args()

    report: dict[str, object] = {
        "model": "QIU_SI_REFERENCE_ANISO",
        "native_sha256": digest(args.native),
        "port_sha256": digest(args.port),
        "rtol": args.rtol,
        "atol": args.atol,
        "fields": {},
    }
    passed = True
    with np.load(args.native, allow_pickle=False) as native, np.load(args.port, allow_pickle=False) as port:
        missing = {"native": sorted(set(REQUIRED) - set(native.files)), "port": sorted(set(REQUIRED) - set(port.files))}
        if missing["native"] or missing["port"]:
            raise SystemExit(f"missing required matched-state fields: {missing}")
        for name in REQUIRED:
            left, right = native[name], port[name]
            same_shape = left.shape == right.shape
            finite = bool(np.all(np.isfinite(left)) and np.all(np.isfinite(right)))
            difference = np.abs(left - right) if same_shape else np.array([np.inf])
            scale = np.maximum(np.abs(left), np.abs(right)) if same_shape else np.array([1.0])
            threshold = args.atol + args.rtol * scale
            exact = name in EXACT
            field_pass = bool(
                same_shape and finite
                and (np.array_equal(left, right) if exact else np.all(difference <= threshold))
            )
            report["fields"][name] = {
                "shape": list(left.shape),
                "max_abs": float(np.max(difference)),
                "max_scaled_error": float(np.max(difference / threshold)),
                "comparison": "exact" if exact else "tolerance",
                "pass": field_pass,
            }
            passed &= field_pass
    report["pass"] = passed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not passed:
        raise SystemExit("QIU_SI_4REF_A0_PORT matched-state correspondence failed")


if __name__ == "__main__":
    main()
