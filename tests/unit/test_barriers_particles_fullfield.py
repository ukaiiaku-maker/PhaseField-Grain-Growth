import numpy as np

from grain_growth_pf.disconnections.barriers import assign_barriers, renew_barrier
from grain_growth_pf.disconnections.spectrum import isotropic_surrogate_library
from grain_growth_pf.mechanics.qiu_full_field import QiuFullField
from grain_growth_pf.obstacles.particles import ParticleField


def test_barrier_disorder_is_quenched_and_renewal_is_explicit():
    modes = isotropic_surrogate_library(b_shells=(0.2,), directions=4, step_heights=(0.2,))
    a = assign_barriers(modes, "truncated_gaussian", 4, 0.5, 0.2, (0.3, 0.7))
    b = assign_barriers(modes, "truncated_gaussian", 4, 0.5, 0.2, (0.3, 0.7))
    assert [m.barrier_ev for m in a] == [m.barrier_ev for m in b]
    assert all(0.3 <= m.barrier_ev <= 0.7 for m in a)
    renewed = renew_barrier(a[0], np.random.default_rng(8), 0.1, (0.3, 0.7))
    assert 0.3 <= renewed.barrier_ev <= 0.7


def test_full_field_is_nonlocal_zero_mean_and_sign_reversing():
    backend = QiuFullField((24, 24))
    increment = np.array([[0.0, 0.5], [0.5, 0.0]])
    backend.add_event((12, 12), increment)
    stress = backend.solve().copy()
    assert np.max(np.abs(stress.mean(axis=(-2, -1)))) < 1e-12
    assert np.count_nonzero(np.abs(stress[0, 1]) > 1e-10) > 10
    opposite = QiuFullField((24, 24)); opposite.add_event((12, 12), -increment)
    assert np.allclose(opposite.solve(), -stress)


def test_full_field_self_stress_opposes_its_source_eigenstrain():
    backend = QiuFullField((24, 24))
    increment = np.array([[0.2, 0.3], [0.3, -0.1]])
    backend.add_event((12, 12), increment)
    backend.solve()
    self_work = float(np.sum(backend.stress * backend.eigenstrain))
    assert self_work < 0.0
    assert -0.5 * self_work > 0.0


def test_zener_pressure_scaling():
    field1 = ParticleField.random(10, 1.0, (20, 20), 1)
    field2 = ParticleField.random(10, 2.0, (20, 20), 1)
    assert np.isclose(field1.zener_pressure_3d(0.1, 1.0), 2 * field2.zener_pressure_3d(0.1, 1.0))
    assert np.isclose(field1.zener_pressure_3d(0.2, 1.0), 2 * field1.zener_pressure_3d(0.1, 1.0))
