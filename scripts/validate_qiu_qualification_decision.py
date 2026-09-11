#!/usr/bin/env python3
"""Validate the required structure and terminal evidence in a QIU decision."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


REQUIRED_KEYS = {
    "schema_version", "created_at", "source_branch", "source_commit",
    "historical_source_commit", "historical_run_path", "initial_state_hashes",
    "selected_backend", "model_identity", "classification", "root_causes",
    "rejected_hypotheses", "mathematical_gates", "source_mapping_gates",
    "timestep_gates", "production_runs", "seed_results", "excluded_runs",
    "tests", "commit_shas", "remaining_limitations",
}
ALLOWED_CLASSIFICATIONS = {
    "NUMERICAL_OR_COUPLING_DEFECT_FIXED",
    "PHYSICAL_QIU_INSTABILITY_QUALIFIED",
    "QIU_REFERENCE_MODEL_MISIDENTIFIED",
}
REQUIRED_RUN_ROLES = {
    "legacy_seed5101", "corrected_seed5101", "refined_seed5101",
    "corrected_seed5102", "corrected_seed5103",
}


def validate(decision: dict[str, object]) -> list[str]:
    errors: list[str] = []
    missing = sorted(REQUIRED_KEYS - decision.keys())
    if missing:
        errors.append(f"missing required keys: {', '.join(missing)}")
    if decision.get("schema_version") != 1:
        errors.append("schema_version must equal 1")
    if decision.get("classification") not in ALLOWED_CLASSIFICATIONS:
        errors.append("classification is not an allowed mission outcome")
    if decision.get("selected_backend") not in {"FFT_EIGENSTRAIN_V2", "QIU_REFERENCE_V2"}:
        errors.append("selected_backend must be an explicit qualified identity")
    source_commit = str(decision.get("source_commit", ""))
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        errors.append("source_commit must be a full 40-character Git SHA")

    runs = decision.get("production_runs", [])
    if not isinstance(runs, list):
        errors.append("production_runs must be a list")
    else:
        roles = {str(run.get("role")) for run in runs if isinstance(run, dict)}
        absent = sorted(REQUIRED_RUN_ROLES - roles)
        if absent:
            errors.append(f"missing production run roles: {', '.join(absent)}")
        for run in runs:
            if not isinstance(run, dict):
                errors.append("every production run must be an object")
                continue
            if run.get("pooling_allowed") is not True:
                errors.append(f"production run {run.get('role')} is not marked pooling_allowed=true")
            if run.get("terminal_status") not in {"completed", "diagnostic_capture"}:
                errors.append(f"production run {run.get('role')} lacks a terminal status")

    tests = decision.get("tests", {})
    if not isinstance(tests, dict):
        errors.append("tests must be an object")
    elif any(int(tests.get(key, -1)) != 0 for key in ("failures", "errors", "skipped")):
        errors.append("final tests must have zero failures, errors, and skips")
    elif int(tests.get("total", 0)) <= 0:
        errors.append("final test total must be positive")

    for key in ("mathematical_gates", "source_mapping_gates", "timestep_gates"):
        gates = decision.get(key)
        if not isinstance(gates, dict) or not gates:
            errors.append(f"{key} must be a nonempty object")
    for key in ("root_causes", "rejected_hypotheses", "commit_shas", "remaining_limitations"):
        if not isinstance(decision.get(key), list):
            errors.append(f"{key} must be a list")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("decision", type=Path)
    args = parser.parse_args()
    decision = json.loads(args.decision.read_text())
    errors = validate(decision)
    if errors:
        raise SystemExit("invalid qualification decision:\n- " + "\n- ".join(errors))
    print(f"valid qualification decision: {args.decision.resolve()}")


if __name__ == "__main__":
    main()
