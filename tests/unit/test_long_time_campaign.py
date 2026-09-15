from __future__ import annotations

from pathlib import Path

import yaml

from grain_growth_pf.campaign import enumerate_campaign


def test_phase1_long_time_matrix_and_scaling_are_frozen():
    path = Path(__file__).parents[2] / "configs/production/long_time_kinetics_900K.yaml"
    spec = yaml.safe_load(path.read_text())
    configs = enumerate_campaign(spec)
    assert [config.regime for config in configs] == [
        "B0", "G", "T", "GT", "GTC_GB", "GSC_GBTJ_Ks025",
        "GTSC_GBTJ_Ks025", "QIU", "GTSC_GBTJ_Ks030_LONG",
    ]
    assert {config.seed for config in configs} == {5101}
    assert {config.pf.temperature for config in configs} == {900.0}
    assert {config.pf.shape for config in configs} == {(384, 384)}
    assert {config.termination_grains for config in configs} == {100}
    assert {config.max_steps for config in configs} == {100000}
    assert {config.output_cadence for config in configs} == {200}
    base_area = 192**2 / 200
    large_area = 384**2 / 800
    assert abs(large_area / base_area - 1.0) < 0.02
    initial_controls = {
        (
            config.parameters["initial_grains"],
            config.parameters["equilibrate_to_grains"],
            config.parameters["event_domain_length"],
        )
        for config in configs
    }
    assert initial_controls == {(1068, 800, 12.0)}
    assert all(config.parameters["event_trace_sample_fraction"] == 0.10 for config in configs)
    assert all(config.parameters["checkpoint_cadence"] == 500 for config in configs)
    assert all(config.parameters["energy_diagnostic_cadence"] == 100 for config in configs)


def test_phase1_shear_values_and_minimum_model_are_explicit():
    path = Path(__file__).parents[2] / "configs/production/long_time_kinetics_900K.yaml"
    configs = {config.regime: config for config in enumerate_campaign(yaml.safe_load(path.read_text()))}
    assert configs["GSC_GBTJ_Ks025"].parameters["shear_stiffness"] == 0.25
    assert configs["GTSC_GBTJ_Ks025"].parameters["shear_stiffness"] == 0.25
    assert configs["GTSC_GBTJ_Ks030_LONG"].parameters["shear_stiffness"] == 0.30
    assert "tj_compatibility" not in configs["GSC_GBTJ_Ks025"].active_modules
    assert "tj_compatibility" in configs["GTSC_GBTJ_Ks025"].active_modules
