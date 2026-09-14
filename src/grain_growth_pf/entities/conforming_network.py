"""Upper-envelope interfaces of piecewise-linear phase fields on a triangle mesh.

Used for energy/geometry qualification. This reconstruction alone does not
establish a work-conjugate map from network forces to the historical PF kernel.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np

from .ordered_network import order_network


def triangle_interfaces(values, tolerance=1e-12):
    """Return (grain_i, grain_j, barycentric_A, barycentric_B) upper-envelope edges."""
    values = np.asarray(values, float)
    if values.ndim != 2 or values.shape[1] != 3 or np.any(~np.isfinite(values)):
        raise ValueError("finite phase values at three vertices required")
    candidates = np.flatnonzero(np.max(values, axis=1) > tolerance)
    if len(candidates) < 2:
        return []
    result = []
    eye = np.eye(3)
    for i, j in combinations(candidates, 2):
        d = values[i]-values[j]
        if np.all(np.abs(d) < tolerance):
            # Coincident maximal phases make the interface nonunique.
            if np.any(np.abs(values[i]-values.max(axis=0)) < tolerance):
                raise ValueError("coincident maximal phases: unresolved reconstruction")
            continue
        crossings = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if abs(d[a]) <= tolerance:
                crossings.append(eye[a])
            if d[a]*d[b] < 0:
                fraction = d[a]/(d[a]-d[b])
                crossings.append((1-fraction)*eye[a]+fraction*eye[b])
        if len(crossings) < 2:
            continue
        a, b = max(combinations(crossings, 2), key=lambda ab: np.linalg.norm(ab[0]-ab[1]))
        if np.linalg.norm(a-b) < tolerance:
            continue
        low, high = 0., 1.
        for k in candidates:
            if k in (i, j):
                continue
            q = values[i]-values[k]
            qa, slope = q@a, q@(b-a)
            if abs(slope) < tolerance:
                if qa < -tolerance:
                    high = -1.
                    break
            elif slope > 0:
                low = max(low, -qa/slope)
            else:
                high = min(high, -qa/slope)
        if high-low > tolerance:
            result.append((int(i), int(j), a+low*(b-a), a+high*(b-a)))
    return result


def reconstruct_periodic(eta, dx=1.):
    """Conforming periodic graph with exact shared TJ vertices; Cartesian output.

    Fixed SW--NE diagonal per cell is a discretization choice requiring rotation
    and grid qualification. Never interpret this polygonal length as continuum
    perimeter without convergence evidence.
    """
    eta = np.asarray(eta, float)
    if eta.ndim != 3 or min(eta.shape[1:]) < 3 or dx <= 0:
        raise ValueError("periodic phase field (P,H,W), H,W>=3 and positive dx required")
    if np.any(~np.isfinite(eta)) or np.any(eta < 0):
        raise ValueError("finite nonnegative phase field required")
    _, height, width = eta.shape
    box = np.array([width, height])*dx
    node_index, vertices, edge_map = {}, [], {}
    tolerance = 1e-10*dx
    bucket_period = np.rint(box/tolerance).astype(np.int64)
    def node(point):
        point = np.mod(point, box)
        key_array = np.floor(point/tolerance).astype(np.int64) % bucket_period
        # Equality of rounded coordinates is not a tolerance predicate:
        # two roundoff-separated evaluations can straddle a rounding bin.
        # Check neighboring periodic buckets and actual minimum-image distance.
        for offset_x in (-1, 0, 1):
            for offset_y in (-1, 0, 1):
                nearby = tuple((key_array+[offset_x, offset_y]) % bucket_period)
                for index in node_index.get(nearby, ()):
                    delta = point-vertices[index]
                    delta -= np.round(delta/box)*box
                    if np.linalg.norm(delta) <= tolerance:
                        return index
        key = tuple(key_array)
        index = len(vertices)
        node_index.setdefault(key, []).append(index)
        vertices.append(point)
        return index
    for y in range(height):
        for x in range(width):
            coords = np.array([[x, y], [x+1, y], [x+1, y+1], [x, y+1]], float)*dx
            values = eta[:, [y, y, (y+1)%height, (y+1)%height],
                         [x, (x+1)%width, (x+1)%width, x]]
            if len(np.unique(np.argmax(values, axis=0))) == 1:
                continue
            for indices in ([0, 1, 2], [0, 2, 3]):
                active = np.flatnonzero(np.max(values[:, indices], axis=1) > 1e-12)
                for i, j, a, b in triangle_interfaces(values[active][:, indices]):
                    pair = tuple(sorted((int(active[i]), int(active[j]))))
                    u, v = node(a@coords[indices]), node(b@coords[indices])
                    if u == v:
                        continue
                    edge = tuple(sorted((u, v)))
                    if edge in edge_map and edge_map[edge] != pair:
                        raise ValueError("multiple pairs occupy a geometric edge")
                    edge_map[edge] = pair
    if not edge_map:
        return np.empty((0, 2)), np.empty((0, 2), int), np.empty((0, 2), int), []
    incident = {}
    for edge, pair in edge_map.items():
        for v in edge:
            incident.setdefault(v, []).append(pair)
    junctions = {v for v, ps in incident.items() if len(set(ps)) > 1}
    for v, ps in incident.items():
        if v in junctions:
            grains = set(g for pair in ps for g in pair)
            if len(ps) != 3 or len(set(ps)) != 3 or len(grains) != 3:
                raise ValueError("unresolved higher-order junction")
        elif len(ps) != 2:
            raise ValueError(f"unclosed periodic boundary at {vertices[v].tolist()}: incident pairs {ps}")
    vertices = np.asarray(vertices)
    edges = np.asarray(list(edge_map), int)
    pairs = np.asarray(list(edge_map.values()), int)
    paths = order_network(vertices, edges, pairs, box=box, junctions=junctions)
    return vertices, edges, pairs, paths
