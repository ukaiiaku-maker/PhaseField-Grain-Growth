"""First variation of a piecewise-linear anisotropic boundary network.

Points use Cartesian (x,y), NOT the legacy tracker's (row,column) order.
An edge A->B has tangent t and clockwise normal n=(t_y,-t_x).
Interior forces and endpoint forces are disjoint terms of one energy gradient.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .anisotropy import BoundaryLaw


@dataclass(frozen=True)
class PolylineVariation:
    energy: float
    edge_lengths: np.ndarray
    line_vectors: np.ndarray
    integrated_forces: np.ndarray
    dual_lengths: np.ndarray
    normals: np.ndarray
    normal_pressure: np.ndarray
    endpoint_forces: np.ndarray


def variation(points, law: BoundaryLaw, phi_i, phi_j, *, closed=False,
              closing_displacement=None):
    """Exact negative nodal gradient; closed loops may wind around a torus.

    The last edge ends at points[0]+closing_displacement. Closed input does
    not duplicate the first node. Points must already be connected/unwrapped.
    Endpoint forces are returned separately and excluded from interior pressure.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < (3 if closed else 2):
        raise ValueError("ordered 2-D polyline with enough vertices required")
    if np.any(~np.isfinite(points)):
        raise ValueError("nonfinite polyline")
    if closing_displacement is not None and not closed:
        raise ValueError("only a closed periodic loop has a winding displacement")
    edges = np.diff(points, axis=0)
    if closed:
        winding = np.zeros(2) if closing_displacement is None else np.asarray(closing_displacement, float)
        if winding.shape != (2,) or np.any(~np.isfinite(winding)):
            raise ValueError("invalid winding displacement")
        edges = np.vstack((edges, points[0]+winding-points[-1]))
    lengths = np.linalg.norm(edges, axis=1)
    if np.any(lengths <= 1e-14):
        raise ValueError("zero-length edge: resolve topology before taking its variation")
    tangent = edges / lengths[:, None]
    edge_normal = np.column_stack((tangent[:, 1], -tangent[:, 0]))
    theta = np.arctan2(edge_normal[:, 1], edge_normal[:, 0])
    gamma = law.evaluate(theta, phi_i, phi_j)[0]
    line = law.vectors(theta, phi_i, phi_j)[1]
    forces = np.zeros_like(points)
    dual = np.zeros(len(points))
    normals = np.zeros_like(points)
    endpoints = np.zeros_like(points)
    if closed:
        forces = line - np.roll(line, 1, axis=0)
        dual = (lengths + np.roll(lengths, 1))/2
        normals = edge_normal + np.roll(edge_normal, 1, axis=0)
    else:
        forces[1:-1] = line[1:] - line[:-1]
        dual[1:-1] = (lengths[1:] + lengths[:-1])/2
        dual[0], dual[-1] = lengths[0]/2, lengths[-1]/2
        normals[1:-1] = edge_normal[1:] + edge_normal[:-1]
        normals[0], normals[-1] = edge_normal[0], edge_normal[-1]
        endpoints[0], endpoints[-1] = line[0], -line[-1]
    norms = np.linalg.norm(normals, axis=1)
    if np.any(norms < 1e-12):
        raise ValueError("unresolved 180-degree cusp")
    normals /= norms[:, None]
    pressure = np.sum(forces*normals, axis=1)/dual
    return PolylineVariation(float(np.dot(gamma, lengths)), lengths, line,
                             forces, dual, normals, pressure, endpoints)


def herring_force(outward_tangents, pairs, orientations, law):
    tangents = np.asarray(outward_tangents, float)
    pairs = np.asarray(pairs, int)
    if tangents.shape != (3, 2) or pairs.shape != (3, 2):
        raise ValueError("a triple junction requires exactly three incident pairs")
    if np.any(pairs < 0) or np.any(pairs >= len(orientations)) or np.any(pairs[:, 0] == pairs[:, 1]):
        raise ValueError("invalid grain pair")
    if np.any(~np.isfinite(tangents)) or np.any(np.linalg.norm(tangents, axis=1) < 1e-14):
        raise ValueError("nonzero finite outward tangents required")
    if len(set(map(tuple, np.sort(pairs, axis=1)))) != 3 or len(np.unique(pairs)) != 3:
        raise ValueError("three distinct grain pairs of three grains required")
    theta = np.arctan2(-tangents[:, 0], tangents[:, 1])
    return law.vectors(theta, np.asarray(orientations)[pairs[:, 0]],
                       np.asarray(orientations)[pairs[:, 1]])[1].sum(axis=0)
