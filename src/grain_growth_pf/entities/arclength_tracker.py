from __future__ import annotations

from itertools import combinations

import numpy as np
from numpy.typing import NDArray

from .gb_segment import GBSegment
from .tracker import EntityTracker, GeometrySnapshot


_NEIGHBORS = tuple(
    (dy, dx)
    for dy in (-1, 0, 1)
    for dx in (-1, 0, 1)
    if dy or dx
)


def _wrapped(point: tuple[int, int], shape: tuple[int, int], periodic: bool) -> tuple[int, int] | None:
    y, x = point
    if periodic:
        return y % shape[0], x % shape[1]
    if 0 <= y < shape[0] and 0 <= x < shape[1]:
        return y, x
    return None


def _periodic_delta(a: np.ndarray, b: np.ndarray, shape: tuple[int, int], periodic: bool) -> np.ndarray:
    delta = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    if periodic:
        box = np.asarray(shape, dtype=float)
        delta -= np.round(delta / box) * box
    return delta


def _component_order(points: np.ndarray, shape: tuple[int, int], periodic: bool) -> list[np.ndarray]:
    """Return 8-connected GB pixel chains ordered approximately by arclength.

    The old tracker lexicographically sorted all pixels belonging to one grain
    pair. That can place remote pieces of the same boundary in one kinetic
    domain. Here each connected component is treated independently and walked
    along nearest connected pixels before it is split into finite physical
    arclength domains.
    """
    if not len(points):
        return []
    unique = np.unique(np.rint(points).astype(int), axis=0)
    lookup = {tuple(point): index for index, point in enumerate(unique)}
    adjacency: list[list[int]] = [[] for _ in range(len(unique))]
    for index, point in enumerate(unique):
        y, x = map(int, point)
        for dy, dx in _NEIGHBORS:
            neighbor = _wrapped((y + dy, x + dx), shape, periodic)
            if neighbor is None:
                continue
            other = lookup.get(neighbor)
            if other is not None and other != index:
                adjacency[index].append(other)

    components: list[list[int]] = []
    unseen = set(range(len(unique)))
    while unseen:
        start = min(unseen)
        stack = [start]
        unseen.remove(start)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(component)

    ordered_components: list[np.ndarray] = []
    for component in components:
        component_set = set(component)
        endpoints = [
            index for index in component
            if sum(neighbor in component_set for neighbor in adjacency[index]) <= 1
        ]
        start = min(endpoints or component, key=lambda index: tuple(unique[index]))
        remaining = set(component)
        remaining.remove(start)
        order = [start]
        previous_direction: np.ndarray | None = None
        current = start
        while remaining:
            local = [neighbor for neighbor in adjacency[current] if neighbor in remaining]
            if local:
                if previous_direction is None:
                    chosen = min(local, key=lambda index: tuple(unique[index]))
                else:
                    def continuation(index: int) -> tuple[float, tuple[int, int]]:
                        direction = -_periodic_delta(
                            unique[current], unique[index], shape, periodic
                        )
                        norm = float(np.linalg.norm(direction))
                        score = (
                            float(previous_direction @ direction / norm)
                            if norm > 0 else -np.inf
                        )
                        return (-score, tuple(unique[index]))
                    chosen = min(local, key=continuation)
            else:
                chosen = min(
                    remaining,
                    key=lambda index: (
                        float(np.linalg.norm(_periodic_delta(
                            unique[index], unique[current], shape, periodic
                        ))),
                        tuple(unique[index]),
                    ),
                )
            direction = -_periodic_delta(unique[current], unique[chosen], shape, periodic)
            norm = float(np.linalg.norm(direction))
            if norm > 0:
                previous_direction = direction / norm
            order.append(chosen)
            remaining.remove(chosen)
            current = chosen
        ordered_components.append(unique[np.asarray(order, dtype=int)].astype(float))

    ordered_components.sort(key=lambda component: tuple(component[0]))
    return ordered_components


class ArclengthEntityTracker(EntityTracker):
    """Entity tracker whose GB kinetic domains are connected arclength intervals."""

    def update(self, labels: NDArray[np.integer]) -> GeometrySnapshot:
        previous = dict(self.boundaries)
        base = super().update(labels)
        shape = tuple(map(int, labels.shape))

        pair_points: dict[tuple[int, int], list[np.ndarray]] = {}
        for segment in base.boundaries.values():
            pair = tuple(sorted((segment.grain_i, segment.grain_j)))
            pair_points.setdefault(pair, []).append(np.asarray(segment.points, dtype=float))

        boundaries: dict[str, GBSegment] = {}
        keys_by_pair: dict[tuple[int, int], set[str]] = {}
        for pair in sorted(pair_points):
            raw = np.concatenate(pair_points[pair]) if pair_points[pair] else np.empty((0, 2))
            components = _component_order(raw, shape, self.periodic)
            chunks: list[np.ndarray] = []
            for component in components:
                physical_length = len(component) * self.dx
                n_domains = max(1, int(np.ceil(physical_length / self.domain_length)))
                chunks.extend(
                    chunk for chunk in np.array_split(component, n_domains) if len(chunk)
                )

            for sid, chunk in enumerate(chunks):
                key = f"gb:{pair[0]}-{pair[1]}:{sid}"
                if key in previous:
                    segment = previous[key]
                    # EntityTracker.update already advanced age on reused keys.
                else:
                    segment = GBSegment(pair[0], pair[1], sid)
                    segment.age = 1
                segment.grain_i, segment.grain_j, segment.segment_id = pair[0], pair[1], sid
                segment.points = np.asarray(chunk, dtype=float)
                segment.length = float(len(chunk) * self.dx)
                segment.misorientation = float(
                    abs(self.orientations[pair[0]] - self.orientations[pair[1]])
                )
                boundaries[key] = segment
                keys_by_pair.setdefault(pair, set()).add(key)

        # Reconnect each TJ to the nearest connected arclength domain for each
        # of its three grain-pair arms.
        for tj in base.triple_junctions.values():
            adjoining: set[str] = set()
            position = np.asarray(tj.position, dtype=float)
            for pair in combinations(tj.grain_ids, 2):
                pair = tuple(sorted(pair))
                candidates = keys_by_pair.get(pair, set())
                if not candidates:
                    continue

                def distance_squared(key: str) -> float:
                    points = boundaries[key].points
                    delta = points - position
                    if self.periodic:
                        box = np.asarray(shape, dtype=float)
                        delta -= np.round(delta / box) * box
                    return float(np.min(np.sum(delta * delta, axis=1)))

                adjoining.add(min(candidates, key=lambda key: (distance_squared(key), key)))
            tj.adjoining_boundaries = adjoining

        self.boundaries = boundaries
        return GeometrySnapshot(base.grains, boundaries, base.triple_junctions)
