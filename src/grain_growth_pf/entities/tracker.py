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
        edge_points, edge_a, edge_b = [], [], []
        for axis in (0, 1):
            shifted = np.roll(labels, -1, axis=axis)
            diff = labels != shifted
            if not self.periodic:
                slicer = [slice(None), slice(None)]
                slicer[axis] = -1
                diff[tuple(slicer)] = False
            points = np.argwhere(diff)
            edge_points.append(points.astype(float))
            edge_a.append(labels[diff].astype(int))
            edge_b.append(shifted[diff].astype(int))
        points = np.concatenate(edge_points) if edge_points else np.empty((0, 2), float)
        side_a = np.concatenate(edge_a) if edge_a else np.empty(0, int)
        side_b = np.concatenate(edge_b) if edge_b else np.empty(0, int)
        low, high = np.minimum(side_a, side_b), np.maximum(side_a, side_b)
        max_label = int(labels.max())
        label_base = max_label + 1
        perimeter_counts = (
            np.bincount(side_a, minlength=max_label + 1)
            + np.bincount(side_b, minlength=max_label + 1)
        )

        grain_ids, grain_counts = np.unique(labels, return_counts=True)
        flat_labels = labels.ravel()
        yy, xx = np.indices(labels.shape)
        if self.periodic:
            centroids = []
            for coordinate, size in ((yy, labels.shape[0]), (xx, labels.shape[1])):
                theta = 2.0 * np.pi * coordinate.ravel() / size
                real = np.bincount(flat_labels, weights=np.cos(theta), minlength=max_label + 1)
                imag = np.bincount(flat_labels, weights=np.sin(theta), minlength=max_label + 1)
                angle = np.angle(real + 1j * imag) % (2.0 * np.pi)
                centroids.append(angle * size / (2.0 * np.pi))
        else:
            centroids = [
                np.bincount(flat_labels, weights=coordinate.ravel(), minlength=max_label + 1)
                / np.maximum(np.bincount(flat_labels, minlength=max_label + 1), 1)
                for coordinate in (yy, xx)
            ]
        for gid, count in zip(grain_ids, grain_counts):
            gid = int(gid)
            grain = self.grains.get(gid, Grain(gid, float(self.orientations[gid])))
            grain.area = float(count * self.dx**2)
            grain.equivalent_radius = float(np.sqrt(grain.area / np.pi))
            grain.centroid = (float(centroids[0][gid]), float(centroids[1][gid]))
            grain.neighbors = set()
            grain.perimeter = float(perimeter_counts[gid] * self.dx)
            seen_grains[gid] = grain

        seen_boundaries: dict[str, GBSegment] = {}
        keys_by_pair: dict[tuple[int, int], set[str]] = {}
        if len(low):
            pair_codes = low * label_base + high
            edge_order = np.argsort(pair_codes, kind="stable")
            sorted_pair_codes = pair_codes[edge_order]
            pair_starts = np.r_[0, np.flatnonzero(np.diff(sorted_pair_codes)) + 1]
            pair_ends = np.r_[pair_starts[1:], len(edge_order)]
            group_order = np.argsort(edge_order[pair_starts])
        else:
            pair_starts = pair_ends = group_order = np.empty(0, int)
        for group in group_order:
            code = int(sorted_pair_codes[pair_starts[group]])
            pair = (code // label_base, code % label_base)
            raw_points = points[edge_order[pair_starts[group]:pair_ends[group]]]
            for gid in pair:
                if gid in seen_grains:
                    seen_grains[gid].neighbors.add(pair[1] if gid == pair[0] else pair[0])
            domain_points = np.asarray(raw_points)
            n_domains = max(1, int(np.ceil(len(domain_points) * self.dx / self.domain_length)))
            order = np.lexsort((domain_points[:, 1], domain_points[:, 0]))
            for sid, indices in enumerate(np.array_split(order, n_domains)):
                key = f"gb:{pair[0]}-{pair[1]}:{sid}"
                seg = self.boundaries.get(key, GBSegment(pair[0], pair[1], sid))
                seg.points = domain_points[indices]
                seg.length = float(len(indices) * self.dx)
                seg.misorientation = float(abs(self.orientations[pair[0]] - self.orientations[pair[1]]))
                seg.age += 1
                seen_boundaries[key] = seg
                keys_by_pair.setdefault(pair, set()).add(key)

        # Robust topology from every 2x2 cell; the triplet key suppresses
        # threshold flicker and a circular mean localizes periodic junctions.
        ny, nx = labels.shape
        y_stop, x_stop = (ny, nx) if self.periodic else (ny - 1, nx - 1)
        base = labels[:y_stop, :x_stop]
        down = np.roll(labels, -1, axis=0)[:y_stop, :x_stop]
        right = np.roll(labels, -1, axis=1)[:y_stop, :x_stop]
        diagonal = np.roll(np.roll(labels, -1, axis=0), -1, axis=1)[:y_stop, :x_stop]
        local = np.sort(np.stack((base, down, right, diagonal), axis=-1), axis=-1).reshape(-1, 4)
        combination_indices = np.asarray(tuple(combinations(range(4), 3)), dtype=int)
        candidates = local[:, combination_indices].reshape(-1, 3)
        cell_indices = np.repeat(np.arange(len(local)), len(combination_indices))
        valid = (candidates[:, 0] != candidates[:, 1]) & (candidates[:, 1] != candidates[:, 2])
        candidates, cell_indices = candidates[valid], cell_indices[valid]
        if len(candidates):
            triplet_codes = (
                (candidates[:, 0] * label_base + candidates[:, 1]) * label_base
                + candidates[:, 2]
            )
            cell_triplet_codes = cell_indices * label_base**3 + triplet_codes
            _, unique_cell_indices = np.unique(cell_triplet_codes, return_index=True)
            unique_cell_indices.sort()
            candidates = candidates[unique_cell_indices]
            cell_indices = cell_indices[unique_cell_indices]
            triplet_codes = triplet_codes[unique_cell_indices]
        cell_y, cell_x = np.divmod(cell_indices, x_stop)
        candidate_positions = np.column_stack((cell_y + 0.5, cell_x + 0.5))
        if len(candidates):
            triplet_order = np.argsort(triplet_codes, kind="stable")
            sorted_triplet_codes = triplet_codes[triplet_order]
            triplet_starts = np.r_[0, np.flatnonzero(np.diff(sorted_triplet_codes)) + 1]
            triplet_ends = np.r_[triplet_starts[1:], len(triplet_order)]
            triplet_group_order = np.argsort(triplet_order[triplet_starts])
        else:
            triplet_starts = triplet_ends = triplet_group_order = np.empty(0, int)
        seen_tj: dict[str, TripleJunction] = {}
        for group in triplet_group_order:
            indices = triplet_order[triplet_starts[group]:triplet_ends[group]]
            raw_tri = candidates[indices[0]]
            tri = tuple(map(int, raw_tri))
            tri_points = candidate_positions[indices]
            position = tuple(tri_points.mean(axis=0))
            key = "tj:" + "-".join(map(str, tri))
            tj = self.triple_junctions.get(key, TripleJunction(tri, position))
            if tj.age:
                delta = np.asarray(position) - np.asarray(tj.position)
                if self.periodic:
                    box = np.asarray(labels.shape, dtype=float)
                    delta -= np.round(delta / box) * box
                tj.travel_distance += float(np.linalg.norm(delta) * self.dx)
            tj.position = position
            tj.adjoining_boundaries = set().union(*(
                keys_by_pair.get(tuple(sorted(pair)), set()) for pair in combinations(tri, 2)
            ))
            tj.age += 1
            seen_tj[key] = tj

        self.grains, self.boundaries, self.triple_junctions = seen_grains, seen_boundaries, seen_tj
        return GeometrySnapshot(self.grains, self.boundaries, self.triple_junctions)
