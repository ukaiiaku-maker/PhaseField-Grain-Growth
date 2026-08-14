import numpy as np

from grain_growth_pf.climb.exchange import butler_volmer_flux, linear_onsager_coefficient
from grain_growth_pf.climb.free_volume import FreeVolumeState
from grain_growth_pf.climb.serial_cycle import SerialClimbCycle
from grain_growth_pf.climb.transport import diffusivity, transport_time
from grain_growth_pf.disconnections.admissibility import combination_rates, feasible_combinations
from grain_growth_pf.disconnections.mode import DisconnectionMode, ModeDriving
from grain_growth_pf.disconnections.spectrum import isotropic_surrogate_library
from grain_growth_pf.mechanics.local_shear_memory import LocalShearMemory


def test_mode_attempt_limit_and_direction_selection():
    plus = DisconnectionMode("plus", (1, 0), 1, 0, 0.5, 1e8, 2,
                             activation_volume_normal=1, activation_volume_shear=1)
    minus = DisconnectionMode("minus", (-1, 0), -1, 0, 0.5, 1e8, 2,
                              activation_volume_normal=-1, activation_volume_shear=-1)
    assert plus.rate(900, ModeDriving(normal_pressure=100)) == 2e8
    assert plus.rate(900, ModeDriving(resolved_shear=0.2)) > minus.rate(900, ModeDriving(resolved_shear=0.2))
    assert plus.rate(900, ModeDriving(normal_pressure=0.2)) > minus.rate(900, ModeDriving(normal_pressure=0.2))


def test_discrete_spectrum_is_isotropic_and_has_minimum_burgers():
    modes = isotropic_surrogate_library(b_shells=(0.2,), directions=8, step_heights=(0.2,), disorder_std_ev=0)
    magnitudes = [np.linalg.norm(m.burgers) for m in modes]
    rates = [m.rate(900) for m in modes]
    assert np.isclose(min(magnitudes), 0.2)
    assert np.ptp(rates) / np.mean(rates) < 1e-14
    assert all(m > 0 for m in magnitudes)


def test_feasible_combinations_are_serial():
    a = DisconnectionMode("a", (1, 0), 1, 0, 0, 2, 1)
    b = DisconnectionMode("b", (-1, 0), 1, 0, 0, 3, 1)
    combos = feasible_combinations([a, b], (0, 0), 2, max_events=2)
    assert len(combos) == 1
    assert np.isclose(combination_rates(combos, 900, ModeDriving())[0], 1 / (1 / 2 + 1 / 3))


def test_free_volume_conservation_and_scaling():
    state = FreeVolumeState(0.12, 0.03, 2.0)
    q1 = state.require_for_area_change(-10)
    q2 = state.require_for_area_change(-20)
    assert q2 == 2 * q1
    assert state.accommodate(q1 / 2) == q1 / 2
    assert state.check_balance()
    zero = FreeVolumeState(0, 0.03, 1)
    assert zero.require_for_area_change(20) == 0


def test_butler_volmer_onsager_and_transport_limits():
    temp, j0 = 1000.0, 2.0
    delta = 1e-8
    assert np.isclose(butler_volmer_flux(delta, temp, j0) / delta,
                      linear_onsager_coefficient(temp, j0), rtol=1e-6)
    d1 = diffusivity(800, 1e-6, 0.8)
    d2 = diffusivity(1000, 1e-6, 0.8)
    assert d2 > d1
    assert transport_time(2, d1) == 4 * transport_time(1, d1)


def test_serial_mean_not_parallel_and_shear_sign():
    assert SerialClimbCycle.mean_completion_time(2, 3, 4) == 1 / 2 + 1 / 3 + 1 / 4
    memory = LocalShearMemory(stiffness=2)
    memory.migrate(beta=0.5, normal_displacement=2)
    assert memory.state == 1
    assert memory.internal_shear_stress == -2
    assert memory.normal_velocity(1, capillary_pressure=0.5, beta=1) < 0
    assert memory.release(0.5) > 0


def test_tj_closed_burgers_sequence_and_combined_driving():
    from grain_growth_pf.entities.triple_junction import TripleJunction
    tj = TripleJunction((1, 2, 3), (0.0, 0.0))
    increments = [np.array([1.0, 0.0]), np.array([-0.5, 0.5]), np.array([-0.5, -0.5])]
    for increment in increments:
        tj.add_burgers(increment)
    assert np.linalg.norm(tj.residual_burgers) < 1e-15

    shear_mode = DisconnectionMode("shear", (1, 0), -1, 0, 0.4, 1e6, 1,
                                    activation_volume_normal=-1, activation_volume_shear=2)
    step_mode = DisconnectionMode("step", (-1, 0), 1, 0, 0.4, 1e6, 1,
                                   activation_volume_normal=2, activation_volume_shear=-1)
    shear_favored = ModeDriving(normal_pressure=0.05, resolved_shear=0.2)
    curvature_favored = ModeDriving(normal_pressure=0.2, resolved_shear=0.05)
    assert shear_mode.rate(900, shear_favored) > step_mode.rate(900, shear_favored)
    assert step_mode.rate(900, curvature_favored) > shear_mode.rate(900, curvature_favored)
