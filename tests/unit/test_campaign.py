from grain_growth_pf.campaign import enumerate_campaign


def test_campaign_regime_overrides_merge_parameters_and_pf():
    spec = {
        "base": {
            "pf": {"shape": [24, 24], "temperature": 800},
            "parameters": {"initial_grains": 8, "shared": 1},
            "max_steps": 10,
        },
        "regimes": {
            "A": {},
            "B": {"parameters": {"shared": 2, "barrier_ev": 0.5},
                  "pf": {"intrinsic_mobility": 3}, "max_steps": 20},
        },
        "temperatures": [700, 900],
        "seeds": [4],
    }
    configs = enumerate_campaign(spec)
    assert len(configs) == 4
    b = next(config for config in configs if config.regime == "B" and config.pf.temperature == 900)
    assert b.parameters == {"initial_grains": 8, "shared": 2, "barrier_ev": 0.5}
    assert b.pf.intrinsic_mobility == 3
    assert b.max_steps == 20
