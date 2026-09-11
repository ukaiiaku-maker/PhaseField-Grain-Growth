"""Strict periodic polylines from an explicit, conforming boundary edge graph.

This is not an adapter for unordered tracker pixels: adjacency and junction
identities must be supplied by a conforming interface reconstruction.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OrderedBoundary:
    pair: tuple[int, int]
    node_ids: tuple[int, ...]
    points: np.ndarray
    closed: bool
    winding: np.ndarray
    endpoint_junctions: tuple[int | None, int | None]


def order_network(vertices, edges, pairs, *, box=None, junctions=()):
    """Partition all edges exactly once, by grain pair and true connectivity.

    Degree-two vertices continue a boundary; marked junctions terminate it.
    A branch without an explicit junction is rejected. No nearest-pixel jumps,
    inferred history transfer, pixel-count length, or proximity TJ attachment.
    Every periodic edge must be shorter than half a box in each coordinate;
    ambiguous half-box edges are rejected. Input vertices are Cartesian.
    """
    vertices = np.asarray(vertices, float)
    edges, pairs = np.asarray(edges, int), np.asarray(pairs, int)
    if vertices.ndim != 2 or vertices.shape[1] != 2 or np.any(~np.isfinite(vertices)):
        raise ValueError("finite Cartesian vertices required")
    if edges.ndim != 2 or edges.shape[1] != 2 or pairs.shape != edges.shape:
        raise ValueError("edges and pairs must both be (N,2)")
    if np.any(edges < 0) or np.any(edges >= len(vertices)) or np.any(pairs < 0):
        raise ValueError("invalid node or grain index")
    if np.any(edges[:, 0] == edges[:, 1]) or np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError("self edges or identical grain pair")
    junctions = set(junctions)
    if not junctions.issubset(range(len(vertices))):
        raise ValueError("unknown junction vertex")
    if box is not None:
        box = np.asarray(box, float)
        if box.shape != (2,) or np.any(~np.isfinite(box)) or np.any(box <= 0):
            raise ValueError("positive Cartesian box lengths required")
    adjacency = {}
    seen = set()
    for edge, pair in zip(edges, pairs):
        pair = tuple(sorted(map(int, pair)))
        a, b = map(int, edge)
        identity = (min(a, b), max(a, b))
        if identity in seen:
            raise ValueError("duplicated physical edge")
        seen.add(identity)
        graph = adjacency.setdefault(pair, {})
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
    result = []
    for pair, graph in sorted(adjacency.items()):
        if any(len(neighbors) > 2 for node, neighbors in graph.items() if node not in junctions):
            raise ValueError("unmarked boundary branch")
        remaining = {tuple(sorted((a, b))) for a in graph for b in graph[a]}
        while remaining:
            candidates = sorted({v for edge in remaining for v in edge
                                 if len(graph[v]) != 2 or v in junctions})
            start = candidates[0] if candidates else min(min(edge) for edge in remaining)
            ids = [start]
            while True:
                current = ids[-1]
                options = sorted(b for b in graph[current] if tuple(sorted((current, b))) in remaining)
                if not options:
                    break
                chosen = options[0]
                remaining.remove(tuple(sorted((current, chosen))))
                ids.append(chosen)
                if chosen == start or chosen in junctions or len(graph[chosen]) != 2:
                    break
            closed = ids[-1] == start
            raw = vertices[ids]
            delta = np.diff(raw, axis=0)
            if box is not None:
                delta -= np.round(delta/box)*box
                if np.any(np.isclose(np.abs(delta), box/2, rtol=0, atol=1e-12)):
                    raise ValueError("ambiguous periodic edge")
            if np.any(np.linalg.norm(delta, axis=1) < 1e-14):
                raise ValueError("zero-length edge")
            points = np.vstack((raw[0], raw[0]+np.cumsum(delta, axis=0)))
            winding = points[-1]-points[0] if closed else np.zeros(2)
            endpoints = (None, None) if closed else tuple(v if v in junctions else None for v in (ids[0], ids[-1]))
            result.append(OrderedBoundary(pair, tuple(ids[:-1] if closed else ids),
                                          points[:-1] if closed else points,
                                          closed, winding, endpoints))
    return result
