from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from typing import Any

import numpy as np


def atomic_write_text(path: str | Path, text: str) -> None:
    """Replace a text file only after its complete contents reach disk."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=destination.parent,
            prefix=f".{destination.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def atomic_savez_compressed(path: str | Path, **arrays: np.ndarray) -> None:
    """Write a compressed NumPy archive without exposing a partial archive."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b", dir=destination.parent,
            prefix=f".{destination.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def save_checkpoint(path: str | Path, eta: np.ndarray, metadata: dict[str, Any],
                    rng: np.random.Generator) -> None:
    state = json.dumps(rng.bit_generator.state)
    atomic_savez_compressed(
        path, eta=eta, metadata=np.asarray(json.dumps(metadata)), rng_state=np.asarray(state)
    )


def load_checkpoint(path: str | Path, rng: np.random.Generator) -> tuple[np.ndarray, dict[str, Any]]:
    with np.load(path, allow_pickle=False) as data:
        eta = data["eta"].copy()
        metadata = json.loads(str(data["metadata"]))
        rng.bit_generator.state = json.loads(str(data["rng_state"]))
    return eta, metadata
