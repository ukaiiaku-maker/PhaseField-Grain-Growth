"""Smooth synthetic fourfold energies; Cartesian (x, y), normal-angle convention.

These independent primitives do not change the historical PF equations.
Normalization must be supplied from a verified initial network and then frozen.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np


@dataclass(frozen=True)
class Strength:
    g_min: float
    inclination_weight: float
    support_power: int
    mobility_exponent: float

    def __post_init__(self):
        if not 0 < self.g_min <= 1 or not 0 <= self.inclination_weight < 1:
            raise ValueError("positive misorientation floor and rounded support required")
        if (not isinstance(self.support_power, int) or self.support_power < 2
                or self.support_power % 2):
            raise ValueError("support power must be an even integer >= 2")
        if not np.isfinite(self.mobility_exponent) or self.mobility_exponent < 0:
            raise ValueError("mobility exponent must be finite and nonnegative")


LADDER = {
    "A0_ISOTROPIC": Strength(1., 0., 2, 0.),
    "A1_MODERATE": Strength(.80, .65, 8, 1.),
    "A2_STRONG": Strength(.65, .85, 16, 2.),
    "A3_STRONGER_BOUNDED": Strength(.55, .90, 24, 2.5),
}


def support_derivatives(psi, power):
    """Return h, h', h'' and h+h'' analytically, including crystal axes."""
    if not isinstance(power, int) or power < 2 or power % 2:
        raise ValueError("even support power required")
    psi = np.asarray(psi, dtype=float)
    if np.any(~np.isfinite(psi)):
        raise ValueError("angles must be finite")
    c, s = np.cos(psi), np.sin(psi)
    if power == 2:
        one = np.ones_like(psi)
        zero = np.zeros_like(psi)
        return one, zero, zero, one
    u = c**power + s**power
    h = u**(1. / power)
    first = u**(1. / power - 1.) * (s**(power-1)*c - c**(power-1)*s)
    # This nonnegative form avoids cancellation at axes, where h+h'' = 0.
    stiffness = (power-1) * (c*s)**(power-2) * u**(1. / power - 2.)
    return h, first, stiffness-h, stiffness


@lru_cache(maxsize=16)
def angular_normalization(strength: Strength):
    theta = np.arange(4096) * (2*np.pi/4096)
    h = support_derivatives(theta, strength.support_power)[0]
    return float(1 / ((1-strength.inclination_weight) + strength.inclination_weight*h.mean()))


def disorientation(phi_i, phi_j):
    """Fourfold fundamental zone [0, pi/4]; radians throughout."""
    return np.abs((np.asarray(phi_i)-np.asarray(phi_j)+np.pi/4) % (np.pi/2)-np.pi/4)


@dataclass(frozen=True)
class BoundaryLaw:
    strength: Strength = LADDER["A2_STRONG"]
    gamma0: float = 1.
    mobility0: float = 4.
    energy_normalization: float = 1.
    mobility_normalization: float = 1.

    def __post_init__(self):
        values = [self.gamma0, self.mobility0, self.energy_normalization,
                  self.mobility_normalization]
        if not all(np.isfinite(x) and x > 0 for x in values):
            raise ValueError("reference values and frozen normalizations must be positive")

    def evaluate(self, theta, phi_i, phi_j):
        """Return gamma, gamma_theta, gamma_theta_theta, stiffness and mobility."""
        theta, phi_i, phi_j = np.broadcast_arrays(theta, phi_i, phi_j)
        if any(np.any(~np.isfinite(x)) for x in (theta, phi_i, phi_j)):
            raise ValueError("orientations and inclinations must be finite")
        a = self.strength
        hi = support_derivatives(theta-phi_i, a.support_power)
        hj = support_derivatives(theta-phi_j, a.support_power)
        g = a.g_min + (1-a.g_min)*np.sin(2*disorientation(phi_i, phi_j))**2
        prefactor = self.gamma0*self.energy_normalization*g*angular_normalization(a)
        weight = a.inclination_weight
        gamma = prefactor*((1-weight) + weight*(hi[0]+hj[0])/2)
        first = prefactor*weight*(hi[1]+hj[1])/2
        stiffness = prefactor*((1-weight) + weight*(hi[3]+hj[3])/2)
        second = stiffness-gamma
        mobility = self.mobility0*self.mobility_normalization*(gamma/self.gamma0)**(-a.mobility_exponent)
        return gamma, first, second, stiffness, mobility

    def vectors(self, theta, phi_i, phi_j):
        gamma, first, *_ = self.evaluate(theta, phi_i, phi_j)
        theta = np.broadcast_to(theta, gamma.shape)
        normal = np.stack((np.cos(theta), np.sin(theta)), axis=-1)
        tangent = np.stack((-np.sin(theta), np.cos(theta)), axis=-1)
        xi = gamma[..., None]*normal + first[..., None]*tangent
        line = gamma[..., None]*tangent - first[..., None]*normal
        return xi, line

    def normalize(self, theta, phi_i, phi_j, lengths):
        """Return a new law; callers must persist it, never refit during evolution."""
        theta, phi_i, phi_j, lengths = np.broadcast_arrays(theta, phi_i, phi_j, lengths)
        if not lengths.size or np.any(~np.isfinite(lengths)) or np.any(lengths <= 0):
            raise ValueError("a nonempty network of positive finite edge lengths is required")
        raw = replace(self, energy_normalization=1., mobility_normalization=1.)
        gamma = raw.evaluate(theta, phi_i, phi_j)[0]
        normalized = replace(raw, energy_normalization=float(self.gamma0/np.average(gamma, weights=lengths)))
        mobility = normalized.evaluate(theta, phi_i, phi_j)[4]
        return replace(normalized, mobility_normalization=float(self.mobility0/np.average(mobility, weights=lengths)))
