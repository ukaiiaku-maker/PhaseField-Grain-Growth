import numpy as np
import pytest

from grain_growth_pf.mechanics.qiu_full_field import (
    FFTEigenstrainV2,
    QiuFullFieldLegacy,
    isotropic_lame_lambda,
)


def _symmetric_random(shape, seed=7):
    raw = np.random.default_rng(seed).normal(size=(2, 2, *shape))
    return 0.5 * (raw + raw.swapaxes(0, 1))


def _spectral_gradient(field, dx=1.0):
    ny, nx = field.shape[-2:]
    waves = (
        2 * np.pi * np.fft.fftfreq(ny, d=dx),
        2 * np.pi * np.fft.fftfreq(nx, d=dx),
    )
    transformed = np.fft.fftn(field, axes=(-2, -1))
    gradient = np.empty((2, *field.shape), dtype=float)
    for axis, wave in enumerate(waves):
        shape = [1] * field.ndim
        shape[-2 + axis] = len(wave)
        gradient[axis] = np.fft.ifftn(
            1j * wave.reshape(shape) * transformed, axes=(-2, -1)
        ).real
    return gradient


@pytest.mark.parametrize("state", ["plane_stress", "plane_strain"])
def test_fft_v2_mechanical_equilibrium_symmetry_and_energy(state):
    backend = FFTEigenstrainV2((17, 19), shear_modulus=2.3,
                               poisson_ratio=0.27, constitutive_state=state)
    backend.set_eigenstrain(_symmetric_random(backend.shape))
    backend.solve()
    assert backend.last_equilibrium_residual < 1e-10
    assert np.max(np.abs(backend.stress - backend.stress.swapaxes(0, 1))) < 1e-13
    assert backend.elastic_energy() >= 0.0


def test_fft_v2_is_linear_sign_reversing_and_translation_invariant():
    source = _symmetric_random((15, 17))
    first = FFTEigenstrainV2((15, 17)); first.set_eigenstrain(source)
    stress = first.solve().copy()
    opposite = FFTEigenstrainV2((15, 17)); opposite.set_eigenstrain(-source)
    assert np.allclose(opposite.solve(), -stress, rtol=2e-14, atol=2e-14)
    shifted = FFTEigenstrainV2((15, 17))
    shifted.set_eigenstrain(np.roll(source, (3, -4), axis=(-2, -1)))
    assert np.allclose(
        shifted.solve(), np.roll(stress, (3, -4), axis=(-2, -1)),
        rtol=2e-13, atol=2e-13,
    )


def test_fft_v2_traction_free_uniform_and_compatible_strain_are_stress_free():
    shape = (15, 17)
    uniform = np.zeros((2, 2, *shape))
    uniform[0, 0] = 0.2; uniform[1, 1] = -0.1
    uniform[0, 1] = uniform[1, 0] = 0.07
    backend = FFTEigenstrainV2(shape)
    backend.set_eigenstrain(uniform)
    assert np.max(np.abs(backend.solve())) < 1e-13

    displacement = np.random.default_rng(11).normal(size=(2, *shape))
    gradient = _spectral_gradient(displacement)
    compatible = 0.5 * (gradient + gradient.swapaxes(0, 1))
    backend.set_eigenstrain(compatible)
    assert np.max(np.abs(backend.solve())) < 2e-12


@pytest.mark.parametrize("seed", [19, 23, 29])
def test_fft_v2_energy_directional_derivative(seed):
    shape = (13, 15)
    source = 0.1 * _symmetric_random(shape, seed)
    direction = _symmetric_random(shape, seed + 1)
    backend = FFTEigenstrainV2(shape, shear_modulus=1.7, poisson_ratio=0.21)
    backend.set_eigenstrain(source); backend.solve()
    predicted = -float(np.sum(backend.stress * direction))
    increment = 2e-7
    energies = []
    for sign in (-1.0, 1.0):
        trial = FFTEigenstrainV2(shape, shear_modulus=1.7, poisson_ratio=0.21)
        trial.set_eigenstrain(source + sign * increment * direction)
        trial.solve(); energies.append(trial.elastic_energy())
    observed = (energies[1] - energies[0]) / (2 * increment)
    relative = abs(observed - predicted) / max(abs(predicted), 1e-12)
    assert relative < 1e-6


def _dense_spectral_oracle(eigenstrain, mu, nu, state):
    """Independent real-space minimization using a dense spectral B matrix."""
    ny, nx = eigenstrain.shape[-2:]
    points = ny * nx
    derivative = []
    for axis, length in enumerate((ny, nx)):
        matrix = np.empty((points, points))
        wave = 2 * np.pi * np.fft.fftfreq(length)
        for column in range(points):
            basis = np.zeros((ny, nx)); basis.ravel()[column] = 1.0
            transformed = np.fft.fftn(basis)
            shape = (length, 1) if axis == 0 else (1, length)
            matrix[:, column] = np.fft.ifftn(
                1j * wave.reshape(shape) * transformed
            ).real.ravel()
        derivative.append(matrix)
    dy, dx = derivative
    b_matrix = np.zeros((3 * points, 2 * points))
    b_matrix[0:points, 0:points] = dy
    b_matrix[points:2 * points, points:2 * points] = dx
    b_matrix[2 * points:, 0:points] = 0.5 * dx
    b_matrix[2 * points:, points:2 * points] = 0.5 * dy
    lam = isotropic_lame_lambda(mu, nu, state)
    point_c = np.asarray([
        [lam + 2 * mu, lam, 0.0],
        [lam, lam + 2 * mu, 0.0],
        [0.0, 0.0, 4 * mu],
    ])
    c_matrix = np.kron(point_c, np.eye(points))
    star = np.concatenate((
        eigenstrain[0, 0].ravel(), eigenstrain[1, 1].ravel(),
        eigenstrain[0, 1].ravel(),
    ))
    stiffness = b_matrix.T @ c_matrix @ b_matrix
    rhs = b_matrix.T @ c_matrix @ star
    displacement = np.linalg.lstsq(stiffness, rhs, rcond=1e-12)[0]
    elastic = b_matrix @ displacement - star
    generalized_stress = c_matrix @ elastic
    result = np.empty_like(eigenstrain)
    result[0, 0] = generalized_stress[:points].reshape(ny, nx)
    result[1, 1] = generalized_stress[points:2 * points].reshape(ny, nx)
    result[0, 1] = result[1, 0] = (
        0.5 * generalized_stress[2 * points:].reshape(ny, nx)
    )
    return result


@pytest.mark.parametrize("state", ["plane_stress", "plane_strain"])
def test_fft_v2_matches_independent_dense_periodic_equilibrium(state):
    shape = (3, 5)
    source = _symmetric_random(shape, 37)
    source -= source.mean(axis=(-2, -1), keepdims=True)
    mu, nu = 1.4, 0.22
    backend = FFTEigenstrainV2(shape, shear_modulus=mu, poisson_ratio=nu,
                               constitutive_state=state)
    backend.set_eigenstrain(source)
    observed = backend.solve()
    expected = _dense_spectral_oracle(source, mu, nu, state)
    assert np.allclose(observed, expected, rtol=2e-11, atol=2e-11)


def test_fft_v2_single_fourier_mode_satisfies_equilibrium_mode_by_mode():
    shape = (17, 19)
    yy, xx = np.indices(shape)
    phase = 2 * np.pi * (3 * yy / shape[0] - 4 * xx / shape[1])
    source = np.zeros((2, 2, *shape))
    source[0, 0] = np.cos(phase)
    source[0, 1] = source[1, 0] = 0.3 * np.cos(phase)
    backend = FFTEigenstrainV2(shape); backend.set_eigenstrain(source)
    stress = backend.solve()
    stress_hat = np.fft.fftn(stress, axes=(-2, -1))
    ky = 2 * np.pi * np.fft.fftfreq(shape[0])
    kx = 2 * np.pi * np.fft.fftfreq(shape[1])
    wave = np.stack(np.meshgrid(ky, kx, indexing="ij"))
    divergence = np.einsum("jYX,ijYX->iYX", wave, stress_hat)
    assert np.max(np.abs(divergence)) < 2e-11


def test_legacy_projection_is_preserved_but_fails_equilibrium_gate():
    backend = QiuFullFieldLegacy((17, 19))
    backend.eigenstrain = _symmetric_random(backend.shape, 41)
    backend.solve()
    assert backend.last_equilibrium_residual > 0.1
