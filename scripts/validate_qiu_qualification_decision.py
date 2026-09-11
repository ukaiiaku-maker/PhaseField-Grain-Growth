#!/usr/bin/env python3
"""Validate the required structure and terminal evidence in a QIU decision."""

from __future__ import annotations

import argparse
import hashlib
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
    "artifacts", "canonical_preservation", "pull_request",
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
REQUIRED_SEEDS = {5101, 5102, 5103}
REQUIRED_MOVIE_ROLES = {"legacy", "corrected_seed5101", "corrected_additional_seed"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_sha(value: object, length: int) -> bool:
    return bool(re.fullmatch(rf"[0-9a-f]{{{length}}}", str(value)))


def _validate_artifact(
    record: object, label: str, errors: list[str], *, check_paths: bool,
) -> None:
    if not isinstance(record, dict):
        errors.append(f"artifact {label} must be an object")
        return
    raw_path = str(record.get("path", ""))
    if not raw_path or not Path(raw_path).is_absolute():
        errors.append(f"artifact {label} path must be absolute")
    if not _is_sha(record.get("sha256"), 64):
        errors.append(f"artifact {label} requires a SHA-256")
    if check_paths and raw_path:
        path = Path(raw_path)
        if not path.is_file():
            errors.append(f"artifact {label} does not exist: {path}")
        elif _is_sha(record.get("sha256"), 64) and _sha256(path) != record["sha256"]:
            errors.append(f"artifact {label} SHA-256 mismatch: {path}")


def validate(decision: dict[str, object], *, check_paths: bool = False) -> list[str]:
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
            if not _is_sha(run.get("source_commit"), 40):
                errors.append(f"production run {run.get('role')} requires a full source commit")
            if not _is_sha(run.get("archive_sha256"), 64):
                errors.append(f"production run {run.get('role')} requires an archive SHA-256")
            local_path = str(run.get("local_path", ""))
            if not local_path or not Path(local_path).is_absolute():
                errors.append(f"production run {run.get('role')} requires an absolute local path")
            elif check_paths and not Path(local_path).is_dir():
                errors.append(f"production run {run.get('role')} local path does not exist")
            if int(run.get("terminal_step", 0)) <= 0:
                errors.append(f"production run {run.get('role')} requires a positive terminal step")

    tests = decision.get("tests", {})
    if not isinstance(tests, dict):
        errors.append("tests must be an object")
    elif any(int(tests.get(key, -1)) != 0 for key in ("failures", "errors", "skipped")):
        errors.append("final tests must have zero failures, errors, and skips")
    elif int(tests.get("total", 0)) <= 0:
        errors.append("final test total must be positive")
    else:
        _validate_artifact(tests.get("junit"), "tests.junit", errors, check_paths=check_paths)

    for key in ("mathematical_gates", "source_mapping_gates", "timestep_gates"):
        gates = decision.get(key)
        if not isinstance(gates, dict) or not gates:
            errors.append(f"{key} must be a nonempty object")
            continue
        for name, gate in gates.items():
            passed = gate if isinstance(gate, bool) else (
                gate.get("passed") if isinstance(gate, dict) else None
            )
            if passed is not True:
                errors.append(f"{key}.{name} is not recorded as passed")

    for key in (
        "root_causes", "rejected_hypotheses", "excluded_runs",
        "commit_shas", "remaining_limitations",
    ):
        value = decision.get(key)
        if not isinstance(value, list) or not value:
            errors.append(f"{key} must be a nonempty list")
    commits = decision.get("commit_shas", [])
    if isinstance(commits, list):
        for commit in commits:
            if not _is_sha(commit, 40):
                errors.append("every commit_shas entry must be a full Git SHA")

    seeds = decision.get("seed_results")
    if not isinstance(seeds, list):
        errors.append("seed_results must be a list")
    else:
        covered = {
            int(result.get("seed")) for result in seeds
            if isinstance(result, dict) and str(result.get("seed", "")).isdigit()
        }
        missing_seeds = sorted(REQUIRED_SEEDS - covered)
        if missing_seeds:
            errors.append(f"seed_results missing seeds: {', '.join(map(str, missing_seeds))}")
        for result in seeds:
            if not isinstance(result, dict):
                errors.append("every seed_results entry must be an object")
            elif result.get("terminal") is not True:
                errors.append(f"seed result {result.get('seed')} is not terminal")

    initial_hashes = decision.get("initial_state_hashes")
    if not isinstance(initial_hashes, dict):
        errors.append("initial_state_hashes must be an object")
    else:
        for seed in REQUIRED_SEEDS:
            if not _is_sha(initial_hashes.get(str(seed)), 64):
                errors.append(f"initial_state_hashes missing valid SHA-256 for seed {seed}")

    artifacts = decision.get("artifacts")
    if not isinstance(artifacts, dict):
        errors.append("artifacts must be an object")
    else:
        for key in ("report", "run_matrix", "analysis_summary"):
            _validate_artifact(artifacts.get(key), key, errors, check_paths=check_paths)
        for key, minimum in (("movies", 3), ("plots", 1)):
            records = artifacts.get(key)
            if not isinstance(records, list) or len(records) < minimum:
                errors.append(f"artifacts.{key} must contain at least {minimum} records")
                continue
            for index, record in enumerate(records):
                _validate_artifact(record, f"{key}[{index}]", errors, check_paths=check_paths)
        movies = artifacts.get("movies", [])
        if isinstance(movies, list):
            roles = {str(movie.get("role")) for movie in movies if isinstance(movie, dict)}
            missing_roles = sorted(REQUIRED_MOVIE_ROLES - roles)
            if missing_roles:
                errors.append(f"movie artifacts missing roles: {', '.join(missing_roles)}")

    canonical = decision.get("canonical_preservation")
    if not isinstance(canonical, dict) or canonical.get("verified_unchanged") is not True:
        errors.append("canonical_preservation must record verified_unchanged=true")
    elif isinstance(canonical, dict):
        _validate_artifact(
            canonical.get("verification_manifest"),
            "canonical_preservation.verification_manifest",
            errors,
            check_paths=check_paths,
        )

    pull_request = decision.get("pull_request")
    if not isinstance(pull_request, dict):
        errors.append("pull_request must be an object")
    else:
        if pull_request.get("state") not in {"OPEN", "DRAFT"}:
            errors.append("pull_request must be OPEN or DRAFT")
        if pull_request.get("merged") is not False:
            errors.append("pull_request must record merged=false")
        if not re.fullmatch(r"https://github\.com/[^/]+/[^/]+/pull/[0-9]+", str(pull_request.get("url", ""))):
            errors.append("pull_request requires a GitHub pull URL")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("decision", type=Path)
    args = parser.parse_args()
    decision = json.loads(args.decision.read_text())
    errors = validate(decision, check_paths=True)
    if errors:
        raise SystemExit("invalid qualification decision:\n- " + "\n- ".join(errors))
    print(f"valid qualification decision: {args.decision.resolve()}")


if __name__ == "__main__":
    main()
