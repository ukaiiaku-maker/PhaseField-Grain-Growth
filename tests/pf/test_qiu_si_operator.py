import numpy as np
import pytest

from grain_growth_pf.pf.qiu_si import (
    QiuSI4RefParameters,
    QiuSIControl,
    qiu_native_beta,
    qiu_native_capillary,
    qiu_native_laplacian,
    qiu_pair_anisotropy_correction,
    qiu_si_pairwise_rate,
)


def _native_ordered_loop(phi, elastic, eij, parameters):
    phases = len(phi)
    lap = np.stack([qiu_native_laplacian(field, parameters.dx) for field in phi])
    present = phi > 0.0
    local = present | np.roll(present, 1, axis=1) | np.roll(present, -1, axis=1) | np.roll(present, 1, axis=2) | np.roll(present, -1, axis=2)
    rate = np.zeros_like(phi)
    for i in range(phases):
        for j in range(phases):
            if i == j:
                continue
            ppp = qiu_native_capillary(
                phi[i], phi[j], lap[i], lap[j], parameters.interface_width
            )
            barrier = (
                np.pi / parameters.interface_width
                * np.sqrt(phi[i] * phi[j]) * eij[i, j]
            )
            rate[i] += (local[i] & local[j]) * parameters.native_mobility * (
                ppp + elastic[i, j] + barrier
            )
    return parameters.dt * rate


def test_archived_parameter_identity_and_beta_values():
    p = QiuSI4RefParameters()
    assert (p.shape, p.phases, p.dx, p.dy, p.dt, p.steps) == (
        (500, 500), 17, 1.0, 1.0, 0.1, 200_000
    )
    assert p.reference_angle == pytest.approx(np.pi / 4)
    expected = 2 * np.tan(np.deg2rad(22.6) / 2)
    assert qiu_native_beta(0.0, np.deg2rad(-22.6)) == pytest.approx(
        (expected, -expected)
    )
    assert qiu_native_beta(0.0, np.deg2rad(45.0)) == (0.0, -0.0)


def test_a0_port_matches_archived_ordered_pair_algebra_exactly():
    rng = np.random.default_rng(31)
    phi = rng.uniform(0.1, 1.0, (3, 5, 6))
    phi /= phi.sum(axis=0)
    orientations = np.array([0.0, -0.2, 0.4])
    elastic = np.zeros((3, 3, 5, 6))
    elastic[0, 1] = rng.normal(scale=0.01, size=(5, 6))
    elastic[1, 0] = -elastic[0, 1]
    elastic[0, 2] = rng.normal(scale=0.01, size=(5, 6))
    elastic[2, 0] = -elastic[0, 2]
    elastic[1, 2] = rng.normal(scale=0.01, size=(5, 6))
    elastic[2, 1] = -elastic[1, 2]
    eij = np.array([[0.0, -20.0, -20.0], [20.0, 0.0, 0.0], [20.0, 0.0, 0.0]])
    p = QiuSI4RefParameters(shape=(5, 6), phases=3)
    actual, audit = qiu_si_pairwise_rate(
        phi, orientations, elastic, eij, QiuSIControl.A0_PORT, p
    )
    expected = _native_ordered_loop(phi, elastic, eij, p)
    np.testing.assert_array_equal(actual, expected)
    assert audit["phase_sum_rate_max_abs"] <= 1e-14


def test_a0_port_matches_native_local_phase_support():
    phi = np.zeros((3, 6, 7))
    phi[0] = 1.0
    phi[1, 1:3, 1:4] = 0.4
    phi[0, 1:3, 1:4] -= 0.4
    phi[2, 4:6, 4:7] = 0.3
    phi[0, 4:6, 4:7] -= 0.3
    elastic = np.zeros((3, 3, 6, 7))
    eij = np.zeros((3, 3))
    p = QiuSI4RefParameters(shape=(6, 7), phases=3)
    actual, _ = qiu_si_pairwise_rate(
        phi, [0.0, -0.2, 0.4], elastic, eij, QiuSIControl.A0_PORT, p
    )
    np.testing.assert_array_equal(actual, _native_ordered_loop(phi, elastic, eij, p))


def test_anisotropic_correction_force_is_its_discrete_energy_derivative():
    rng = np.random.default_rng(32)
    first = rng.uniform(0.1, 0.8, (4, 5))
    second = rng.uniform(0.1, 0.8, (4, 5))
    direction = rng.normal(size=(4, 5))
    arguments = (0.1, 0.7, 1.0, 0.65, 0.85, 16, 1.0, 1.0)
    energy, drive = qiu_pair_anisotropy_correction(first, second, *arguments)
    epsilon = 1e-7
    plus = qiu_pair_anisotropy_correction(first + epsilon * direction, second, *arguments)[0]
    minus = qiu_pair_anisotropy_correction(first - epsilon * direction, second, *arguments)[0]
    measured = (plus - minus) / (2 * epsilon)
    assert measured == pytest.approx(-np.sum(drive * direction), rel=2e-7, abs=2e-7)
    assert np.isfinite(energy)


def test_pair_exchange_is_antisymmetric_and_mobility_scales_complete_drive():
    phi = np.empty((2, 4, 5))
    phi[0] = 0.4
    phi[1] = 0.6
    orientations = np.array([0.0, 0.3])
    elastic = np.zeros((2, 2, 4, 5))
    elastic[0, 1] = 0.7
    elastic[1, 0] = -0.7
    eij = np.array([[0.0, -20.0], [20.0, 0.0]])
    p = QiuSI4RefParameters(shape=(4, 5), phases=2, dt=1e-5)
    a0, _ = qiu_si_pairwise_rate(
        phi, orientations, elastic, eij, QiuSIControl.A0_PORT, p
    )
    mobile, audit = qiu_si_pairwise_rate(
        phi, orientations, elastic, eij, QiuSIControl.ANISO_M, p
    )
    np.testing.assert_allclose(np.sum(mobile, axis=0), 0.0, atol=1e-15)
    ratio = audit["mobility_min"] / p.native_mobility
    np.testing.assert_allclose(mobile, a0 * ratio, rtol=2e-15, atol=1e-15)


def test_rejects_nonantisymmetric_native_pair_inputs():
    phi = np.full((2, 2, 2), 0.5)
    elastic = np.zeros((2, 2, 2, 2))
    elastic[0, 1] = 1.0
    with pytest.raises(ValueError, match="elastic"):
        qiu_si_pairwise_rate(
            phi, np.zeros(2), elastic, np.zeros((2, 2)),
            QiuSIControl.A0_PORT, QiuSI4RefParameters(shape=(2, 2), phases=2),
        )
