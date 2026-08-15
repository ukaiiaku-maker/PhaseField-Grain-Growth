import numpy as np
import pytest

from grain_growth_pf.pf.free_energy import laplacian
from grain_growth_pf.pf.kernels import pairwise_free_energy, pairwise_obstacle_step


def _vectorized_reference(eta, active, scale, external, dt, periodic):
    indices = np.flatnonzero(active)
    values = eta[indices]
    boundary = "periodic" if periodic else "neumann"
    lap = laplacian(values, 1.0, boundary)
    present = values > 1e-14
    if periodic:
        support = (
            present
            | np.roll(present, 1, -2)
            | np.roll(present, -1, -2)
            | np.roll(present, 1, -1)
            | np.roll(present, -1, -1)
        )
    else:
        padded = np.pad(present, ((0, 0), (1, 1), (1, 1)), mode="constant")
        support = (
            present
            | padded[:, 2:, 1:-1]
            | padded[:, :-2, 1:-1]
            | padded[:, 1:-1, 2:]
            | padded[:, 1:-1, :-2]
        )
    count = np.maximum(support.sum(axis=0, keepdims=True), 1)
    phase_sum = (values * support).sum(axis=0, keepdims=True)
    lap_sum = (lap * support).sum(axis=0, keepdims=True)
    obstacle = 2.0 * np.sin(np.pi / 8.0) ** 2
    rate = 4.0 * (
        lap * phase_sum - values * lap_sum
        + obstacle * (count * values - phase_sum)
    ) * support
    ext = external[indices].copy()
    ext -= (ext * support).sum(axis=0, keepdims=True) / count
    rate += 4.0 * ext * support
    trial = np.clip(values + dt * rate * scale[None], 0.0, 1.0)
    result = eta.copy()
    result[indices] = trial / trial.sum(axis=0, keepdims=True)
    return result


@pytest.mark.parametrize("periodic", [False, True])
def test_compiled_update_matches_vectorized_equation(periodic):
    rng = np.random.default_rng(741)
    labels = rng.integers(0, 5, size=(13, 11))
    eta = np.eye(5)[labels].transpose(2, 0, 1).astype(float)
    # Add compact diffuse pixels while retaining the filling constraint.
    eta = 0.8 * eta + 0.2 * np.roll(eta, 1, axis=2)
    active = np.array([True, True, False, True, True])
    eta[2] = 0.0
    eta[0] += np.maximum(1.0 - eta.sum(axis=0), 0.0)
    scale = rng.uniform(0.2, 1.0, labels.shape)
    external = rng.normal(0.0, 0.03, eta.shape)
    expected = _vectorized_reference(eta, active, scale, external, 0.01, periodic)
    actual = pairwise_obstacle_step(
        eta, active, scale, external, True, 0.01, 4.0, 1.0, 4.0, 1.0, periodic
    )
    assert np.allclose(actual, expected, rtol=2e-15, atol=2e-15)


@pytest.mark.parametrize("periodic", [False, True])
def test_compiled_energy_matches_vectorized_definition(periodic):
    rng = np.random.default_rng(992)
    eta = rng.random((7, 12, 9))
    eta /= eta.sum(axis=0, keepdims=True)
    if periodic:
        gx = np.roll(eta, -1, axis=1) - eta
        gy = np.roll(eta, -1, axis=2) - eta
    else:
        gx = np.zeros_like(eta)
        gy = np.zeros_like(eta)
        gx[:, :-1] = eta[:, 1:] - eta[:, :-1]
        gy[:, :, :-1] = eta[:, :, 1:] - eta[:, :, :-1]
    phase_sum = eta.sum(axis=0)
    pair_potential = 0.5 * (phase_sum**2 - np.sum(eta**2, axis=0))
    pair_gradient = 0.5 * (
        gx.sum(axis=0) ** 2 + gy.sum(axis=0) ** 2
        - np.sum(gx**2 + gy**2, axis=0)
    )
    expected = float((pair_potential - 16.0 / np.pi**2 * pair_gradient).sum())
    actual = pairwise_free_energy(eta, 1.0, 4.0, 1.0, periodic)
    assert np.isclose(actual, expected, rtol=2e-15, atol=2e-13)
