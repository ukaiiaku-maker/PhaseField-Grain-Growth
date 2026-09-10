import numpy as np
import pytest

from grain_growth_pf.mechanics.fft_eigenstrain_coupling import (
    LocalInterfaceSweepCoupling,
)
from grain_growth_pf.mechanics.qiu_full_field import FFTEigenstrainV2


def _one_hot(labels, phases):
    return np.eye(phases, dtype=float)[labels].transpose(2, 0, 1)


def _stripe(shape=(24, 28)):
    labels = np.zeros(shape, dtype=int)
    labels[:, shape[1] // 4:3 * shape[1] // 4] = 1
    return _one_hot(labels, 2)


def _circle(shape=(41, 41), radius=10, phases=2):
    yy, xx = np.indices(shape)
    cy, cx = (np.asarray(shape) - 1) / 2
    labels = ((yy - cy) ** 2 + (xx - cx) ** 2 <= radius**2).astype(int)
    return _one_hot(labels, phases)


def test_planar_translation_uses_each_local_transfer_once():
    before = _stripe()
    after = np.roll(before, 1, axis=2)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    coupling.source_increment(before, after, np.asarray([0.1, 0.8]))
    expected = float(np.maximum(after - before, 0.0).sum())
    assert coupling.last_diagnostics.absolute_swept_area == expected
    assert expected == 2 * before.shape[1]


def test_shrinking_circle_and_multineighbor_transfer_are_not_duplicated():
    before = _circle(radius=10)
    after = _circle(radius=9)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    coupling.source_increment(before, after, np.asarray([0.2, 1.0]))
    lost_area = float(np.count_nonzero(before[1] - after[1] > 0.5))
    assert coupling.last_diagnostics.absolute_swept_area == lost_area

    labels_before = np.asarray([
        [1, 1, 0, 2, 2],
        [1, 1, 0, 2, 2],
        [3, 3, 0, 4, 4],
        [3, 3, 0, 4, 4],
    ])
    labels_after = labels_before.copy()
    labels_after[0, 2] = 1; labels_after[1, 2] = 2
    labels_after[2, 2] = 3; labels_after[3, 2] = 4
    before_multi = _one_hot(labels_before, 5)
    after_multi = _one_hot(labels_after, 5)
    coupling.source_increment(
        before_multi, after_multi, np.asarray([0.0, 0.2, 0.5, 0.9, 1.2])
    )
    assert coupling.last_diagnostics.absolute_swept_area == 4.0


def test_disconnected_same_pair_segments_are_counted_by_pixels_not_pair_key():
    labels = np.zeros((20, 24), dtype=int)
    labels[3:7, 3:7] = 1; labels[12:16, 15:19] = 1
    before = _one_hot(labels, 2)
    labels_after = labels.copy()
    labels_after[3:7, 3] = 0; labels_after[12:16, 15] = 0
    after = _one_hot(labels_after, 2)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    coupling.source_increment(before, after, np.asarray([0.0, 0.6]))
    assert coupling.last_diagnostics.absolute_swept_area == 8.0


def test_exact_rigid_translation_and_stationary_boundary_create_no_source():
    before = _stripe()
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=2)
    rigid = coupling.source_increment(
        before, np.roll(before, 1, axis=2), np.asarray([0.1, 0.8])
    )
    assert coupling.last_diagnostics.rigid_translation_removed
    assert np.count_nonzero(rigid) == 0
    stationary = coupling.source_increment(before, before.copy(), np.asarray([0.1, 0.8]))
    assert coupling.last_diagnostics.absolute_swept_area == 0.0
    assert np.count_nonzero(stationary) == 0


def test_reversing_migration_reverses_the_source_exactly():
    before = _circle(radius=9)
    after = _circle(radius=8)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    forward = coupling.source_increment(before, after, np.asarray([0.2, 0.9]))
    reverse = coupling.source_increment(after, before, np.asarray([0.2, 0.9]))
    assert np.array_equal(reverse, -forward)


def test_symmetric_opposing_sweeps_cancel_integrated_source():
    labels = np.zeros((31, 31), dtype=int)
    labels[7:14, 5:12] = 1
    labels[17:24, 19:26] = 1
    before = _one_hot(labels, 2)
    after_labels = labels.copy()
    after_labels[7:14, 5] = 0; after_labels[7:14, 12] = 1
    after_labels[17:24, 25] = 0; after_labels[17:24, 18] = 1
    after = _one_hot(after_labels, 2)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    increment = coupling.source_increment(before, after, np.asarray([0.1, 0.8]))
    integrated = increment.sum(axis=(-2, -1))
    assert np.max(np.abs(integrated)) < 1e-14


@pytest.mark.parametrize("dx", [1.0, 0.5, 0.25])
def test_planar_swept_area_is_exact_under_grid_refinement(dx):
    physical_height, physical_width = 12.0, 16.0
    shape = (round(physical_height / dx), round(physical_width / dx))
    labels = np.zeros(shape, dtype=int)
    left, right = round(4.0 / dx), round(12.0 / dx)
    labels[:, left:right] = 1
    before = _one_hot(labels, 2)
    pixels = round(1.0 / dx)
    after = np.roll(before, pixels, axis=2)
    coupling = LocalInterfaceSweepCoupling(dx, rigid_shift_search=0)
    coupling.source_increment(before, after, np.asarray([0.1, 0.8]))
    expected = 2.0 * physical_height * 1.0
    relative = abs(coupling.last_diagnostics.absolute_swept_area - expected) / expected
    assert relative < 1e-3


def test_source_and_feedback_are_discrete_work_conjugates():
    eta = _stripe((17, 19))
    orientations = np.asarray([0.15, 0.85])
    elasticity = FFTEigenstrainV2((17, 19), shear_modulus=1.6, poisson_ratio=0.24)
    rng = np.random.default_rng(71)
    base = rng.normal(scale=0.02, size=elasticity.eigenstrain.shape)
    base = 0.5 * (base + base.swapaxes(0, 1))
    elasticity.set_eigenstrain(base); elasticity.solve()

    y, x, amount = 5, 3, 1e-7
    after = eta.copy(); after[0, y, x] -= amount; after[1, y, x] += amount
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    increment = coupling.source_increment(eta, after, orientations)
    predicted_derivative = -float(np.sum(elasticity.stress * increment)) / amount
    driving = coupling.driving_field(eta, orientations, elasticity.stress)
    force_derivative = -float(driving[1, y, x] - driving[0, y, x])
    assert abs(predicted_derivative) > 1e-5
    assert np.isclose(predicted_derivative, force_derivative, rtol=1e-6, atol=1e-10)

    alpha = 2e-5
    energies = []
    direction = increment / amount
    for sign in (-1.0, 1.0):
        trial = FFTEigenstrainV2((17, 19), shear_modulus=1.6, poisson_ratio=0.24)
        trial.set_eigenstrain(base + sign * alpha * direction)
        trial.solve(); energies.append(trial.elastic_energy())
    observed = (energies[1] - energies[0]) / (2 * alpha)
    relative = abs(observed - predicted_derivative) / max(abs(predicted_derivative), 1e-12)
    assert relative < 1e-4


def test_accumulated_source_persists_across_topology_loss():
    before = _circle(radius=5)
    after = _circle(radius=4)
    coupling = LocalInterfaceSweepCoupling(rigid_shift_search=0)
    increment = coupling.source_increment(before, after, np.asarray([0.2, 1.0]))
    elasticity = FFTEigenstrainV2(before.shape[1:])
    elasticity.add_field(increment); stored = elasticity.eigenstrain.copy()
    # No mapper operation silently edits already accumulated mechanical state.
    coupling.source_increment(after, _circle(radius=0), np.asarray([0.2, 1.0]))
    assert np.array_equal(elasticity.eigenstrain, stored)
