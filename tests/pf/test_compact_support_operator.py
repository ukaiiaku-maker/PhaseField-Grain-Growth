import numpy as np
import pytest

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.compact_support import exact_candidate_graph
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


def _three_phase_solver(*, tolerance: float = 1e-10) -> MultiphaseFieldSolver:
    height = width = 12
    eta = np.zeros((3, height, width), dtype=float)
    eta[0, :, : width // 2] = 1.0
    eta[1, :, width // 2 :] = 1.0
    eta[2] = 1e-12
    eta /= eta.sum(axis=0, keepdims=True)
    config = PFConfig(
        shape=(height, width),
        time_step=0.002,
        interface_width=4.0,
        anisotropy_strength="A2_STRONG",
        anisotropic_support_mode="compact_active_set",
        anisotropic_kkt_tolerance=tolerance,
    )
    return MultiphaseFieldSolver(
        eta, config, orientations=np.array([0.0, 0.3, 0.7])
    )


def test_compact_support_retires_uniform_tiny_tail_at_exact_zero():
    solver = _three_phase_solver()
    diagnostics = solver.step().compact_support
    assert diagnostics is not None
    assert np.count_nonzero(solver.eta[2]) == 0
    assert diagnostics.retirements >= solver.config.shape[0] * solver.config.shape[1]
    assert diagnostics.active_max <= 2
    assert diagnostics.kkt_residual <= solver.config.anisotropic_kkt_tolerance
    np.testing.assert_array_equal(solver.eta.sum(axis=0), np.ones(solver.config.shape))


def test_candidate_graph_has_no_remote_phase_nucleation():
    eta = np.zeros((3, 9, 9), dtype=float)
    eta[0] = 1.0
    eta[0, 4, 4] = 0.0
    eta[1, 4, 4] = 1.0
    eta[0, 0, 0] = 0.0
    eta[2, 0, 0] = 1.0
    graph, counts = exact_candidate_graph(eta, periodic=False)
    center = set(graph[4, 4, : counts[4, 4]])
    assert center == {0, 1}
    assert 2 not in center


@pytest.mark.parametrize("tolerance", [1e-8, 1e-10, 1e-12])
def test_compact_support_energy_descends_and_tolerance_converges(tolerance):
    solver = _three_phase_solver(tolerance=tolerance)
    energies = [solver._anisotropic_energy()]
    for _ in range(12):
        energies.append(solver.step().interfacial_energy)
    assert np.max(np.diff(energies)) <= 128 * np.finfo(float).eps * max(energies)
    assert np.min(solver.eta) >= 0.0
    assert np.max(np.abs(solver.eta.sum(axis=0) - 1.0)) <= 1e-14


def test_compact_support_kkt_tolerance_ladder_has_same_short_trajectory():
    fields = []
    for tolerance in (1e-8, 1e-10, 1e-12):
        solver = _three_phase_solver(tolerance=tolerance)
        solver.run(12)
        fields.append(solver.eta)
    np.testing.assert_array_equal(fields[0], fields[1])
    np.testing.assert_array_equal(fields[1], fields[2])


def test_compact_support_restart_is_exact_and_energy_cadence_invariant():
    continuous = _three_phase_solver()
    continuous.run(10)

    split = _three_phase_solver()
    for _ in range(5):
        split.step(compute_energy=False)
    state = split.state_dict()
    resumed = _three_phase_solver()
    resumed.load_state_dict(state)
    for _ in range(5):
        resumed.step(compute_energy=True)

    np.testing.assert_array_equal(resumed.eta, continuous.eta)
    assert resumed.time == continuous.time
    assert resumed.step_number == continuous.step_number
