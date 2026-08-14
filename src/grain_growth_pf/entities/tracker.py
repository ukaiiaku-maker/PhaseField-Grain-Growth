from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import NDArray

from .gb_segment import GBSegment
from .grain import Grain
from .triple_junction import TripleJunction


@dataclass
class GeometrySnapshot:
    grains: dict[int, Grain]
    boundaries: dict[str, GBSegment]
    triple_junctions: dict[str, TripleJunction]


def _periodic_centroid(mask: NDArray[np.bool_]) -> tuple[float, float]:
    coords = np.argwhere(mask)
    if not len(coords):
        return (float("nan"), float("nan"))
    result = []
    for axis, size in enumerate(mask.shape):
        theta = 2 * np.pi * coords[:, axis] / size
        mean = np.angle(np.mean(np.exp(1j * theta))) % (2 * np.pi)
        result.append(float(mean * size / (2 * np.pi)))
    return result[0], result[1]


class EntityTracker:
    """Persistent label-keyed grain, GB-domain, and TJ tracker.

    Pixel masks determine geometry but never own kinetic state. Pair/triplet keys
    preserve state across motion; a key that disappears is retired rather than
    transferred. Boundary domains are deterministic chunks of physical length.
    """

    def __init__(self, orientations: NDArray[np.float64], dx: float = 1.0,
                 domain_length: float = 8.0, periodic: bool = True):
        self.orientations = np.asarray(orientations, dtype=float)
        self.dx = dx
        self.domain_length = domain_length
        self.periodic = periodic
        self.grains: dict[int, Grain] = {}
        self.boundaries: dict[str, GBSegment] = {}
        self.triple_junctions: dict[str, TripleJunction] = {}

    def update(self, labels: NDArray[np.integer]) -> GeometrySnapshot:
        labels = np.asarray(labels)
        seen_grains: dict[int, Grain] = {}
        pair_edges: dict[tuple[int, int], list[tuple[float, float]]] = {}
        perimeter: dict[int, float] = {}

        for axis in (0, 1):
            shifted = np.roll(labels, -1, axis=axis)
            diff = labels != shifted
            if not self.periodic:
                slicer = [slice(None), slice(None)]
                slicer[axis] = -1
                diff[tuple(slicer)] = False
            for y, x in np.argwhere(diff):
                a, b = int(labels[y, x]), int(shifted[y, x])
                pair = tuple(sorted((a, b)))
                pair_edges.setdefault(pair, []).append((float(y), float(x)))
                perimeter[a] = perimeter.get(a, 0.0) + self.dx
                perimeter[b] = perimeter.get(b, 0.0) + self.dx

        for gid in np.unique(labels):
            gid = int(gid)
            mask = labels == gid
            grain = self.grains.get(gid, Grain(gid, float(self.orientations[gid])))
            grain.area = float(mask.sum() * self.dx**2)
            grain.equivalent_radius = float(np.sqrt(grain.area / np.pi))
            grain.centroid = _periodic_centroid(mask) if self.periodic else tuple(np.argwhere(mask).mean(axis=0))
            grain.neighbors = set()
            grain.perimeter = perimeter.get(gid, 0.0)
            seen_grains[gid] = grain

        seen_boundaries: dict[str, GBSegment] = {}
        for pair, raw_points in pair_edges.items():
            for gid in pair:
                if gid in seen_grains:
                    seen_grains[gid].neighbors.add(pair[1] if gid == pair[0] else pair[0])
            points = np.asarray(raw_points)
            n_domains = max(1, int(np.ceil(len(points) * self.dx / self.domain_length)))
            order = np.lexsort((points[:, 1], points[:, 0]))
            for sid, indices in enumerate(np.array_split(order, n_domains)):
                key = f"gb:{pair[0]}-{pair[1]}:{sid}"
                seg = self.boundaries.get(key, GBSegment(pair[0], pair[1], sid))
                seg.points = points[indices]
                seg.length = float(len(indices) * self.dx)
                seg.misorientation = float(abs(self.orientations[pair[0]] - self.orientations[pair[1]]))
                seg.age += 1
                seen_boundaries[key] = seg

        # Robust topology from every 2x2 cell; the triplet key suppresses
        # threshold flicker and a circular mean localizes periodic junctions.
        triplet_points: dict[tuple[int, int, int], list[tuple[float, float]]] = {}
        ny, nx = labels.shape
        y_range = range(ny if self.periodic else ny - 1)
        x_range = range(nx if self.periodic else nx - 1)
        for y in y_range:
            for x in x_range:
                local = {int(labels[y, x]), int(labels[(y + 1) % ny, x]),
                         int(labels[y, (x + 1) % nx]), int(labels[(y + 1) % ny, (x + 1) % nx])}
                for tri in combinations(sorted(local), 3):
                    triplet_points.setdefault(tri, []).append((y + 0.5, x + 0.5))
        seen_tj: dict[str, TripleJunction] = {}
        for tri, points in triplet_points.items():
            position = tuple(np.asarray(points).mean(axis=0))
            key = "tj:" + "-".join(map(str, tri))
            tj = self.triple_junctions.get(key, TripleJunction(tri, position))
            if tj.age:
                delta = np.asarray(position) - np.asarray(tj.position)
                if self.periodic:
                    box = np.asarray(labels.shape, dtype=float)
                    delta -= np.round(delta / box) * box
                tj.travel_distance += float(np.linalg.norm(delta) * self.dx)
            tj.position = position
            tj.adjoining_boundaries = {
                k for k, b in seen_boundaries.items()
                if set((b.grain_i, b.grain_j)).issubset(tri)
            }
            tj.age += 1
            seen_tj[key] = tj

        self.grains, self.boundaries, self.triple_junctions = seen_grains, seen_boundaries, seen_tj
        return GeometrySnapshot(self.grains, self.boundaries, self.triple_junctions)

