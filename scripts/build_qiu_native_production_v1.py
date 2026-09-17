#!/usr/bin/env python3
"""Build the exact-restart segmented pristine Qiu-SI production driver."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from scripts.build_qiu_native_preflight_v3 import PRISTINE_DRIVER_SHA256, instrument as preflight_instrument


PRISTINE_FUNCTIONS_SHA256 = "3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"production instrumentation anchor count is {source.count(old)}, expected one: {old[:100]!r}")
    return source.replace(old, new)


def instrument(source: str) -> str:
    generated = preflight_instrument(source)
    generated = replace_once(
        generated,
        "nsteps = int(os.environ.get('QIU_PREFLIGHT_TARGET_STEP', '1000'))",
        "nsteps = int(os.environ.get('QIU_SEGMENT_TARGET_STEP', '200000'))",
    )
    generated = replace_once(
        generated,
        'preflight_output = os.environ.get("QIU_PREFLIGHT_OUTPUT", "preflight-output")\nrestart_path = os.environ.get("QIU_PREFLIGHT_RESTART", "")\nstart_step = 1',
        'preflight_output = os.environ.get("QIU_PRODUCTION_OUTPUT", "production-output")\nrestart_path = os.environ.get("QIU_PRODUCTION_RESTART", "")\ncheckpoint_cadence = int(os.environ.get("QIU_CHECKPOINT_CADENCE", "250"))\nif checkpoint_cadence <= 0: raise ValueError("QIU_CHECKPOINT_CADENCE must be positive")\nstart_step = 1',
    )
    generated = replace_once(
        generated,
        "    if accepted_step in (100, 500):\n        trajectory_state = dict(",
        "    if accepted_step == nsteps:\n        trajectory_state = dict(",
    )
    generated = replace_once(
        generated,
        "    capture_checkpoint = nstep in CHECKPOINT_STEPS",
        "    capture_checkpoint = nstep == nsteps or nstep % checkpoint_cadence == 0",
    )
    generated = generated.replace('manifest["baseline_eligible"]=False', 'manifest["baseline_eligible"]=True')
    generated = generated.replace('manifest["instrumentation"]="qiu-native-preflight-v3"', 'manifest["instrumentation"]="qiu-native-production-v1"\n_manifest["logical_trajectory"]="QIU_SI_NATIVE_BASELINE_V2"\n_manifest["checkpoint_cadence"]=checkpoint_cadence')
    generated = generated.replace('print("QIU_NATIVE_PREFLIGHT_V3_SEGMENT_COMPLETE", nsteps)', 'print("QIU_NATIVE_PRODUCTION_SEGMENT_COMPLETE", start_step, nsteps)')
    return generated


def build(source_dir: Path, instrumentation: Path, output: Path) -> dict:
    driver = source_dir / "Bicrystal-4reference-el-pf.py"
    functions = source_dir / "functions_4ref_new.py"
    if sha256(driver) != PRISTINE_DRIVER_SHA256 or sha256(functions) != PRISTINE_FUNCTIONS_SHA256:
        raise ValueError("pristine source identity mismatch")
    output.mkdir(parents=True, exist_ok=False)
    generated = output / "Bicrystal-4reference-el-pf-production-v1.py"
    generated.write_text(instrument(driver.read_text()))
    shutil.copy2(functions, output / functions.name)
    shutil.copy2(instrumentation, output / instrumentation.name)
    manifest = {
        "schema": "qiu-native-production-v1-build",
        "baseline_eligible": True,
        "logical_trajectory": "QIU_SI_NATIVE_BASELINE_V2",
        "source_variant": "pristine_archive_instrumented_v1",
        "runtime_variant": "anaconda_2022_05_historical_object_mode",
        "execution_threads": 1,
        "pristine_driver_sha256": PRISTINE_DRIVER_SHA256,
        "functions_sha256": PRISTINE_FUNCTIONS_SHA256,
        "instrumented_driver_sha256": sha256(generated),
        "instrumentation_sha256": sha256(output / instrumentation.name),
        "checkpoint_cadence": 250,
        "terminal_step": 200000,
    }
    (output / "production_build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("instrumentation", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.source_dir, args.instrumentation, args.output), indent=2))


if __name__ == "__main__":
    main()
