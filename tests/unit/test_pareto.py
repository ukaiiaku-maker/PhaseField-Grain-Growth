import numpy as np
import pandas as pd

from grain_growth_pf.analysis.pareto import pareto_summary


def test_pareto_ranking_respects_scaling_constraint_and_dominance():
    summary = pd.DataFrame({
        "regime": ["balanced", "jerky", "dominated", "bad-scaling"],
        "temperature": [900.0] * 4,
        "n": [2.0, 2.3, 2.4, 3.0],
        "n_ci_low": [1.9, 2.1, 2.2, 2.9],
        "n_ci_high": [2.1, 2.5, 2.6, 3.1],
        "K": [1.0] * 4,
        "jerkiness_CV": [2.0, 4.0, 3.0, 9.0],
        "Fano": [2.0, 4.0, 3.0, 9.0],
        "burstiness": [0.1, 0.4, 0.2, 0.9],
    })
    result = pareto_summary(summary, exponent_tolerance=0.5).set_index("regime")
    assert result.loc["balanced", "pareto_front"]
    assert result.loc["jerky", "pareto_front"]
    assert result.loc["dominated", "pareto_rank"] == 2
    assert not result.loc["bad-scaling", "physically_admissible"]
    assert pd.isna(result.loc["bad-scaling", "pareto_rank"])


def test_pareto_treats_missing_event_metrics_as_worst_observed():
    summary = pd.DataFrame({
        "regime": ["with-events", "without-events"],
        "temperature": [900.0, 900.0], "n": [2.0, 2.0],
        "n_ci_low": [1.8, 1.8], "n_ci_high": [2.2, 2.2], "K": [1.0, 1.0],
        "jerkiness_CV": [2.0, 2.0], "Fano": [3.0, np.nan],
        "burstiness": [0.2, 0.2],
    })
    result = pareto_summary(summary).set_index("regime")
    assert result.loc["with-events", "pareto_rank"] == 1
    assert result.loc["without-events", "pareto_rank"] == 2
