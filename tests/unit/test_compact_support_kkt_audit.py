import numpy as np

from scripts.audit_compact_support_kkt_tolerance_v2 import T_STAR, TARGETS, comparison


def metrics() -> dict[str, object]:
    return {
        "energy": 10.0, "boundary_density": 0.1, "grain_area_cv": 0.2,
        "area_weighted_mean_radius": 3.0, "fourth_harmonic": 0.01,
        "eighth_harmonic": 0.02, "support_mean": 2.0, "support_p95": 3.0,
        "support_max": 4, "grain_count": 20,
        "exact_active_set_sha256": "a", "candidate_graph_sha256": "b",
        "orientations_sha256": "c", "mobility_scale_sha256": "d",
    }


def test_registered_endpoint_alignment_is_exact():
    # Solver time is accumulated stepwise, so direct multiplication may differ
    # by a few rounding units even though the registered accepted steps match.
    assert abs(TARGETS[6.36493341483619e-05] * 6.36493341483619e-05 - T_STAR) < 5e-14
    assert abs(TARGETS[3.182466707418095e-05] * 3.182466707418095e-05 - T_STAR) < 5e-14


def test_identical_matched_fields_are_bitwise_converged():
    eta = np.zeros((2, 2, 2)); eta[0] = 1.0
    state = {"eta": eta}
    result = comparison(metrics(), metrics(), state, state)
    assert result["exact_eta_equal"]
    assert result["field_rms"] == 0.0
    assert result["non_tie_label_disagreement"] == 0.0
