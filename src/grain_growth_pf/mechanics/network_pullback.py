"""Exact phase-value derivative of a resolved triangular network energy.

This is the work-conjugate geometric pullback, not a qualified PF time update.
Topology-degenerate vertices are rejected rather than regularized silently.
"""
from __future__ import annotations

import numpy as np

from grain_growth_pf.entities.conforming_network import triangle_interfaces
from .anisotropy import BoundaryLaw


def triangle_energy_pullback(values, coordinates, orientations, law: BoundaryLaw):
    values, coordinates = np.asarray(values, float), np.asarray(coordinates, float)
    orientations = np.asarray(orientations, float)
    if values.ndim != 2 or values.shape[1] != 3 or coordinates.shape != (3, 2):
        raise ValueError('phase values at three Cartesian triangle vertices required')
    if (orientations.shape != (len(values),) or np.any(~np.isfinite(coordinates))
            or np.any(~np.isfinite(orientations))):
        raise ValueError('invalid orientations or coordinates')
    basis = np.column_stack((coordinates[1]-coordinates[0], coordinates[2]-coordinates[0]))
    if np.linalg.cond(basis) > 1e10:
        raise ValueError('unresolved triangle')
    spatial_gradients = (values[:, 1:]-values[:, :1]) @ np.linalg.inv(basis)
    gradient = np.zeros_like(values)
    energy = 0.
    for i, j, a, b in triangle_interfaces(values):
        displacement = (b-a)@coordinates
        length = np.linalg.norm(displacement)
        if length < 1e-12:
            raise ValueError('unresolved short edge')
        tangent = displacement/length
        theta = np.arctan2(-tangent[0], tangent[1])
        gamma = law.evaluate(theta, orientations[i], orientations[j])[0]
        line = law.vectors(theta, orientations[i], orientations[j])[1]
        energy += length*float(gamma)
        for barycentric, force in ((a, line), (b, -line)):
            zero = np.flatnonzero(np.abs(barycentric) < 1e-10)
            if len(zero) >= 2:
                raise ValueError('interface at a mesh vertex is not differentiably resolved')
            difference_gradient = spatial_gradients[i]-spatial_gradients[j]
            if len(zero) == 1:
                opposite = int(zero[0])
                direction = coordinates[(opposite+1)%3]-coordinates[(opposite+2)%3]
                denominator = difference_gradient@direction
                if abs(denominator) < 1e-12:
                    raise ValueError('interface tangent to a triangle side')
                coefficient = (force@direction)/denominator
                gradient[i] += coefficient*barycentric
                gradient[j] -= coefficient*barycentric
            else:
                local_values = values@barycentric
                active = np.flatnonzero(np.abs(local_values-local_values.max()) < 1e-10)
                third = [int(k) for k in active if k not in (i, j)]
                if len(third) != 1:
                    raise ValueError('unresolved interior junction')
                k = third[0]
                constraints = np.stack((difference_gradient, spatial_gradients[i]-spatial_gradients[k]))
                if np.linalg.cond(constraints) > 1e10:
                    raise ValueError('singular junction constraints')
                coefficients = np.linalg.solve(constraints.T, force)
                gradient[i] += coefficients.sum()*barycentric
                gradient[j] -= coefficients[0]*barycentric
                gradient[k] -= coefficients[1]*barycentric
    return float(energy), gradient
