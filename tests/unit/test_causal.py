from __future__ import annotations

from grain_growth_pf.analysis.causal import trace_center_causal_null


def test_trace_center_null_detects_release_centered_response():
    profiles = []
    for _ in range(20):
        profile = {step: 1.0 for step in range(-25, 51)}
        for step in range(1, 6):
            profile[step] = 3.0
        profiles.append(profile)
    result = trace_center_causal_null(profiles, shuffles=100, seed=8)
    assert result["events"] == 20
    assert result["actual_mean_delta"] == 2.0
    assert result["causal_excess"] > 1.5
    assert result["p_one_sided"] < 0.02
