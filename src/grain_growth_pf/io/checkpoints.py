from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def save_checkpoint(path: str | Path, eta: np.ndarray, metadata: dict[str, Any],
                    rng: np.random.Generator) -> None:
    state = json.dumps(rng.bit_generator.state)
    np.savez_compressed(path, eta=eta, metadata=json.dumps(metadata), rng_state=state)


def load_checkpoint(path: str | Path, rng: np.random.Generator) -> tuple[np.ndarray, dict[str, Any]]:
    with np.load(path, allow_pickle=False) as data:
        eta = data["eta"].copy()
        metadata = json.loads(str(data["metadata"]))
        rng.bit_generator.state = json.loads(str(data["rng_state"]))
    return eta, metadata

