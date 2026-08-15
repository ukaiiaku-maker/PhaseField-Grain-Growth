from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import scipy


def git_sha(root: str | Path = ".") -> str:
    result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, text=True,
                            capture_output=True, check=False)
    if result.returncode != 0:
        return "UNCOMMITTED"
    sha = result.stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain"], cwd=root, text=True,
                           capture_output=True, check=False)
    return sha + ("-dirty" if dirty.stdout.strip() else "")


def software_versions() -> dict[str, str]:
    return {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__}


def canonical_hash(data: dict[str, Any]) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def file_sha256(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(path: str | Path, config: dict[str, Any], status: str = "started",
                   extra: dict[str, Any] | None = None,
                   code_sha: str | None = None) -> dict[str, Any]:
    # A long-running job must retain its launch revision even if the checkout
    # advances before finalization. Callers may therefore pin the captured SHA.
    code_sha = code_sha or git_sha()
    target = Path(path)
    previous: dict[str, Any] = {}
    if target.exists():
        try:
            previous = json.loads(target.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            previous = {}
    manifest = {**previous,
        "config": config,
        "config_sha256": canonical_hash({"config": config, "git_sha": code_sha}),
        "git_sha": code_sha,
        "software": software_versions(),
        "status": status,
    }
    manifest.update(extra or {})
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
