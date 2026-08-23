from __future__ import annotations

import numpy as np
import pytest

from grain_growth_pf.climb.area_loss import (
    ConservedDefectInventory,
    best_compatible_tj_velocity,
    gb_excess_volume_density,
    released_excess_volume,
    released_point_defect_quota,
    validate_tj_sink_candidate,
)


def test_lower_density_gb_has_positive_excess_volume():
    assert np.isclose(gb_excess_volume_density(
        delta_gb=2.0, rho_lattice=10.0, rho_gb=8.0,
    ), 0.4)


def test_equal_density_gb_has_zero_excess_volume():
    assert gb_excess_volume_density(
        delta_gb=2.0, rho_lattice=10.0, rho_gb=10.0,
    ) == 0.0


def test_decreasing_gb_density_monotonically_increases_excess_volume():
    values = [
        gb_excess_volume_density(delta_gb=1.0, rho_lattice=10.0, rho_gb=rho)
        for rho in (9.0, 8.0, 7.0)
    ]
    assert values[0] < values[1] < values[2]


def test_direct_excess_volume_override_preserves_production_value():
    assert gb_excess_volume_density(
        excess_volume_per_area=0.01,
        delta_gb=-1.0,
        rho_lattice=0.0,
        rho_gb=20.0,
    ) == 0.01


@pytest.mark.parametrize(
    "kwargs",
    [
        {"delta_gb": 1.0, "rho_lattice": 10.0, "rho_gb": 12.0},
        {"delta_gb": 1.0, "rho_lattice": 0.0, "rho_gb": 0.0},
        {"delta_gb": -1.0, "rho_lattice": 10.0, "rho_gb": 8.0},
    ],
)
def test_nonphysical_derived_excess_volume_is_rejected(kwargs):
    with pytest.raises(ValueError):
        gb_excess_volume_density(**kwargs)


def test_straight_gb_translation_has_zero_area_loss_demand():
    assert released_point_defect_quota(12.0, 12.0, 0.01, 0.02) == 0.0


def test_controlled_gb_shortening_matches_analytic_demand():
    # ΔV_ex = 0.01 * (12 - 9) = 0.03 and ΔN = 0.03 / 0.02 = 1.5.
    assert np.isclose(released_excess_volume(12.0, 9.0, 0.01), 0.03)
    assert np.isclose(released_point_defect_quota(12.0, 9.0, 0.01, 0.02), 1.5)


def test_gb_length_increase_releases_no_excess_volume():
    assert released_excess_volume(9.0, 12.0, 0.01) == 0.0
    assert released_point_defect_quota(9.0, 12.0, 0.01, 0.02) == 0.0


def test_domain_split_and_merge_preserve_inventory():
    ledger = ConservedDefectInventory()
    ledger.require(4.0, "parent")
    ledger.split("parent", ["left", "right"], [1.0, 3.0])
    assert np.isclose(ledger.active_vacancy["left"], 1.0)
    assert np.isclose(ledger.active_vacancy["right"], 3.0)
    ledger.merge(["left", "right"], "merged")
    assert np.isclose(ledger.active_vacancy["merged"], 4.0)
    assert ledger.conservation_residual == 0.0


def test_entity_retirement_moves_inventory_without_destroying_it():
    ledger = ConservedDefectInventory()
    ledger.require(2.25, "vanishing-gb")
    ledger.retire_missing({"surviving-gb"})
    assert "vanishing-gb" not in ledger.active_vacancy
    assert np.isclose(ledger.retired_vacancy, 2.25)
    assert np.isclose(ledger.stored_total, 2.25)
    assert ledger.conservation_residual == 0.0


def test_gb_sink_consumes_available_quota_once_without_overshoot():
    ledger = ConservedDefectInventory()
    ledger.require(0.4, "gb:1-2:0")
    accepted = ledger.accommodate(1.0, "gb", "gb:1-2:0")
    assert np.isclose(accepted, 0.4)
    assert np.isclose(ledger.accommodated_gb, 0.4)
    assert ledger.stored_total == 0.0
    assert ledger.accommodate(1.0, "gb", "gb:1-2:0") == 0.0
    assert np.isclose(ledger.accommodated_gb, 0.4)


def test_tj_sink_allows_correct_sign_flux_and_finite_residual_energy():
    decision = validate_tj_sink_candidate(
        available_signed_quota=2.0,
        requested_sign=1,
        compatibility_norm=1e-5,
        compatibility_tolerance=1e-3,
        residual_burgers=np.zeros(2),
        burgers_increment=np.asarray([0.2, 0.0]),
        residual_stiffness=2.0,
    )
    assert decision.allowed
    assert decision.reason == "allowed_finite_residual"
    assert np.isclose(decision.residual_energy_change, 0.04)


def test_tj_sink_forbids_wrong_sign_flux():
    decision = validate_tj_sink_candidate(
        available_signed_quota=2.0,
        requested_sign=-1,
        compatibility_norm=0.0,
        compatibility_tolerance=1e-3,
        residual_burgers=np.zeros(2),
        burgers_increment=np.zeros(2),
        residual_stiffness=1.0,
    )
    assert not decision.allowed
    assert decision.reason == "wrong_sign_or_no_flux"


def test_tj_sink_forbids_incompatible_geometry():
    decision = validate_tj_sink_candidate(
        available_signed_quota=2.0,
        requested_sign=1,
        compatibility_norm=0.2,
        compatibility_tolerance=0.1,
        residual_burgers=np.zeros(2),
        burgers_increment=np.zeros(2),
        residual_stiffness=1.0,
    )
    assert not decision.allowed
    assert decision.reason == "incompatible_geometry"


def test_burgers_incompatibility_can_be_strict_or_finite_penalty():
    arguments = dict(
        available_signed_quota=1.0,
        requested_sign=1,
        compatibility_norm=0.0,
        compatibility_tolerance=0.1,
        residual_burgers=np.asarray([0.1, 0.0]),
        burgers_increment=np.asarray([0.2, 0.0]),
        residual_stiffness=1.0,
    )
    finite = validate_tj_sink_candidate(**arguments)
    strict = validate_tj_sink_candidate(**arguments, strict_burgers_tolerance=0.2)
    assert finite.allowed and np.isfinite(finite.residual_energy_change)
    assert not strict.allowed and strict.reason == "burgers_incompatible"


def test_compatible_tj_velocity_has_zero_residual():
    normals = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    normals[2] /= np.linalg.norm(normals[2])
    velocity = np.asarray([0.3, -0.2])
    desired = normals @ velocity
    fitted, residual, norm = best_compatible_tj_velocity(normals, desired)
    assert np.allclose(fitted, velocity)
    assert np.allclose(residual, 0.0, atol=1e-14)
    assert norm < 1e-14


def test_incompatible_tj_velocity_has_finite_residual():
    normals = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    desired = np.asarray([0.0, 0.0, 1.0])
    _, _, norm = best_compatible_tj_velocity(normals, desired)
    assert norm > 0.0


@pytest.mark.parametrize("dx,pixels", [(1.0, 12), (0.5, 24), (0.25, 48)])
def test_grid_refinement_preserves_physical_area_loss_quota(dx, pixels):
    initial_length = pixels * dx
    final_length = initial_length - 2.0
    quota = released_point_defect_quota(initial_length, final_length, 0.01, 0.02)
    assert np.isclose(quota, 1.0)
