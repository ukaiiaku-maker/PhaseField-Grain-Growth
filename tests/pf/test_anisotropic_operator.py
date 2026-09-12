import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.mechanics.anisotropy import (
    BoundaryLaw, LADDER, angular_normalization,
)
from grain_growth_pf.pf.anisotropic import anisotropic_energy_gradient
from grain_growth_pf.pf.geometry import circular_grain
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


def _gradient_arguments(eta):
    strength = LADDER["A2_STRONG"]
    return (
        eta,
        np.ones(len(eta), dtype=bool),
        np.array([0.1, 0.5, 1.0]),
        1.0,
        4.0,
        4.0,
        1.0,
        True,
        strength.g_min,
        strength.inclination_weight,
        strength.support_power,
        strength.mobility_exponent,
        angular_normalization(strength),
        1.2060485247581734,
        0.924620319961451,
        True,
    )


def test_discrete_anisotropic_force_is_energy_derivative():
    rng = np.random.default_rng(4)
    eta = rng.uniform(0.1, 1.0, (3, 4, 5))
    eta /= eta.sum(axis=0)
    energy, derivative = anisotropic_energy_gradient(*_gradient_arguments(eta))
    index = (1, 2, 3)
    epsilon = 1e-7
    plus = eta.copy()
    minus = eta.copy()
    plus[index] += epsilon
    minus[index] -= epsilon
    plus_energy = anisotropic_energy_gradient(*_gradient_arguments(plus))[0]
    minus_energy = anisotropic_energy_gradient(*_gradient_arguments(minus))[0]
    finite_difference = (plus_energy - minus_energy) / (2.0 * epsilon)
    assert np.isfinite(energy)
    assert abs(derivative[index] - finite_difference) <= 1e-7


def test_pair_energy_is_continuous_when_third_phase_leaves_local_support():
    profile = np.array([0.05, 0.3, 0.7, 0.95])
    two_phase = np.empty((3, 3, 4))
    two_phase[0] = profile
    two_phase[1] = 1.0 - profile
    two_phase[2] = 0.0
    epsilon = 1e-12
    three_phase = two_phase.copy()
    three_phase[:2] *= 1.0 - epsilon
    three_phase[2] = epsilon
    energy_two = anisotropic_energy_gradient(*_gradient_arguments(two_phase))[0]
    energy_three = anisotropic_energy_gradient(*_gradient_arguments(three_phase))[0]
    assert abs(energy_three - energy_two) <= 1e-8


def test_a0_configuration_is_exact_historical_kernel_nesting():
    eta = circular_grain((16, 16), 4, 3)
    common = dict(
        shape=(16, 16), interface_width=3, time_step=0.005,
        intrinsic_mobility=0.2,
    )
    historical = MultiphaseFieldSolver(eta.copy(), PFConfig(**common))
    nested = MultiphaseFieldSolver(
        eta.copy(), PFConfig(**common, anisotropy_strength="A0_ISOTROPIC")
    )
    historical_record = historical.step()
    nested_record = nested.step()
    assert np.array_equal(historical.eta, nested.eta)
    assert historical_record == nested_record


def test_anisotropic_pair_exchange_preserves_phase_sum_before_projection_effects():
    rng = np.random.default_rng(8)
    eta = rng.uniform(0.2, 1.0, (3, 5, 6))
    eta /= eta.sum(axis=0)
    cfg = PFConfig(
        shape=(5, 6), interface_width=4, time_step=1e-8,
        intrinsic_mobility=0.2, anisotropy_strength="A2_STRONG",
        anisotropy_energy_normalization=1.2060485247581734,
        anisotropy_mobility_normalization=0.924620319961451,
    )
    solver = MultiphaseFieldSolver(
        eta.copy(), cfg, orientations=np.array([0.1, 0.5, 1.0])
    )
    initial_energy = solver._anisotropic_energy()
    record = solver.step()
    assert np.max(np.abs(solver.eta.sum(axis=0) - 1.0)) <= 2e-15
    assert record.interfacial_energy < initial_energy


def test_pair_mobility_multiplies_noncapillary_drive():
    eta = np.full((2, 3, 4), 0.5)
    orientations = np.array([0.1, 0.5])
    energy_normalization = 1.2060485247581734
    mobility_normalization = 0.924620319961451
    cfg = PFConfig(
        shape=(3, 4), interface_width=4, time_step=1e-8,
        intrinsic_mobility=4.0, anisotropy_strength="A2_STRONG",
        anisotropic_energy=False, anisotropic_mobility=True,
        anisotropy_energy_normalization=energy_normalization,
        anisotropy_mobility_normalization=mobility_normalization,
        grain_extinction_threshold=0.1,
    )
    external = np.empty_like(eta)
    external[0] = 1.0
    external[1] = -1.0
    solver = MultiphaseFieldSolver(
        eta.copy(), cfg, driving=lambda _eta, _time: external,
        orientations=orientations,
    )
    solver.step(compute_energy=False)
    law = BoundaryLaw(
        strength=LADDER["A2_STRONG"], gamma0=1.0, mobility0=4.0,
        energy_normalization=energy_normalization,
        mobility_normalization=mobility_normalization,
    )
    expected_mobility = law.evaluate(0.0, *orientations)[4]
    measured_rate = (solver.eta[0, 0, 0] - 0.5) / cfg.time_step
    assert np.isclose(measured_rate, expected_mobility, rtol=1e-8)


def test_pair_flux_obstacle_constraint_is_conservative_without_projection():
    eta = np.empty((3, 3, 4))
    eta[0] = 0.1
    eta[1] = 0.2
    eta[2] = 0.7
    external = np.empty_like(eta)
    external[0] = 100.0
    external[1] = 0.0
    external[2] = -100.0
    cfg = PFConfig(
        shape=(3, 4), interface_width=4, time_step=0.1,
        intrinsic_mobility=4.0, adaptive_stepping=False,
        anisotropy_strength="A2_STRONG", grain_extinction_threshold=0.01,
    )
    solver = MultiphaseFieldSolver(
        eta, cfg, driving=lambda _eta, _time: external,
        orientations=np.array([0.1, 0.5, 1.0]),
    )
    solver.step(compute_energy=False)
    assert np.min(solver.eta) >= -1e-14
    assert np.max(np.abs(solver.eta.sum(axis=0) - 1.0)) <= 2e-15
    assert np.any(solver.eta == 0.0)


def test_anisotropic_activity_does_not_delete_positive_subthreshold_phase():
    eta = np.empty((3, 4, 5))
    eta[0] = 0.495
    eta[1] = 0.495
    eta[2] = 0.01
    cfg = PFConfig(
        shape=(4, 5), interface_width=4, time_step=1e-8,
        intrinsic_mobility=0.2, anisotropy_strength="A2_STRONG",
        grain_extinction_threshold=0.05,
    )
    solver = MultiphaseFieldSolver(
        eta, cfg, orientations=np.array([0.1, 0.5, 1.0])
    )
    assert np.array_equal(solver.active_phases, np.ones(3, dtype=bool))
    solver.step()
    assert solver.active_phases[2]
    assert np.max(solver.eta[2]) > 0.0


def test_anisotropic_restart_preserves_orientation_and_path():
    rng = np.random.default_rng(11)
    eta = rng.uniform(0.2, 1.0, (3, 5, 5))
    eta /= eta.sum(axis=0)
    orientations = np.array([0.1, 0.5, 1.0])
    cfg = PFConfig(
        shape=(5, 5), interface_width=4, time_step=1e-6,
        intrinsic_mobility=0.2, anisotropy_strength="A2_STRONG",
    )
    continuous = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    continuous.run(2)
    interrupted = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    interrupted.run(1)
    restored = MultiphaseFieldSolver(eta.copy(), cfg, orientations=orientations)
    restored.load_state_dict(interrupted.state_dict())
    restored.run(1)
    assert np.array_equal(continuous.eta, restored.eta)
    assert np.array_equal(continuous.orientations, restored.orientations)
    assert continuous.time == restored.time
