#!/usr/bin/env python3
"""Create a checksummed source attestation from one immutable HPC3 ledger."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from grain_growth_pf.io.checkpoints import atomic_write_text


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_source_commit(entrypoint: list[object]) -> str:
    candidates = [
        str(value) for value in entrypoint
        if re.fullmatch(r"[0-9a-f]{40}", str(value))
    ]
    if len(candidates) != 1:
        raise ValueError(f"expected exactly one full source commit in entrypoint, found {candidates}")
    return candidates[0]


def create_attestation(ledger_path: Path) -> dict[str, object]:
    ledger = json.loads(ledger_path.read_text())
    archive = Path(ledger["input_bundle"])
    archive_hash = sha256(archive)
    if archive_hash != ledger["input_bundle_checksum"]:
        raise ValueError("HPC3 input archive no longer matches its ledger checksum")

    source_bundle = Path(ledger["local_source_path"]) / "source.bundle"
    source_hash = sha256(source_bundle)
    expected_source_hash = ledger["input_checksums"]["source.bundle"]
    if source_hash != expected_source_hash:
        raise ValueError("source.bundle no longer matches its immutable input checksum")
    entrypoint = list(ledger["manifest"]["entrypoint"])
    source_commit = extract_source_commit(entrypoint)
    heads = subprocess.run(
        ["git", "bundle", "list-heads", str(source_bundle)], text=True,
        capture_output=True, check=True,
    ).stdout.splitlines()
    if not any(line.split()[0] == source_commit for line in heads if line.split()):
        raise ValueError("entrypoint source commit is absent from source.bundle")
    return {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "source_identity_attested_independently_of_application_cwd",
        "hpc3_run_id": ledger["run_id"],
        "slurm_job_id": ledger["job_id"],
        "verified_source_commit": source_commit,
        "source_bundle_sha256": source_hash,
        "hpc3_input_archive_sha256": archive_hash,
        "entrypoint": entrypoint,
        "bundle_heads": heads,
        "verification": [
            "local HPC3 input archive equals the durable ledger checksum",
            "source.bundle equals the archived per-input checksum",
            "exactly one 40-character commit is asserted by the immutable entrypoint",
            "the asserted commit is present in the verified Git bundle",
            "the HPC wrapper checks out that commit detached and rejects dirty state before execution",
        ],
        "pooling_allowed": False,
        "reason": (
            "This attests source identity only. It does not mark a running, failed, "
            "or partial trajectory scientifically poolable."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ledger", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    attestation = create_attestation(args.ledger)
    atomic_write_text(args.output, json.dumps(attestation, indent=2) + "\n")
    print(args.output.resolve())


if __name__ == "__main__":
    main()
