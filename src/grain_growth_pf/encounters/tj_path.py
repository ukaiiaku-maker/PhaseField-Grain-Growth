from __future__ import annotations

import numpy as np


def periodic_path_increment(previous: tuple[float, float], current: tuple[float, float],
                            box: tuple[float, float]) -> float:
    delta = np.asarray(current, float) - np.asarray(previous, float)
    size = np.asarray(box, float)
    delta -= np.round(delta / size) * size
    return float(np.linalg.norm(delta))

