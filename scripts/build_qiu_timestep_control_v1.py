#!/usr/bin/env python3
"""Build the matched-physical-time half-dt combined Qiu control."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

try:
    from scripts.build_qiu_control_v1 import build as build_control
except ModuleNotFoundError:
    from build_qiu_control_v1 import build as build_control


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(source: Path, instrumentation: Path, output: Path) -> dict:
    manifest = build_control(source, instrumentation, output)
    driver = output / "Bicrystal-4reference-el-pf-control-v1.py"
    text = driver.read_text()
    anchor = "dt = 0.1                            # Time Step"
    if text.count(anchor) != 1:
        raise ValueError("timestep instrumentation anchor is not unique")
    text = text.replace(anchor, "dt = float(os.environ.get('QIU_CONTROL_DT', '0.05'))  # Time Step")
    driver.write_text(text)
    manifest.update({
        "schema": "qiu-timestep-control-v1-build", "control": "QIU_SI_4REF_ANISO_EM_INV_DT_HALF",
        "accepted_dt": 0.05, "matched_reference_steps": 5000, "target_steps": 10000,
        "instrumented_driver_sha256": sha256(driver),
    })
    (output / "control_build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("source",type=Path);parser.add_argument("instrumentation",type=Path);parser.add_argument("output",type=Path);args=parser.parse_args()
    print(json.dumps(build(args.source,args.instrumentation,args.output),indent=2))


if __name__ == "__main__":
    main()
