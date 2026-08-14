#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import kstest

from grain_growth_pf.analysis.activation_energy import fit_activation_energy
from grain_growth_pf.disconnections.mode import DisconnectionMode, K_B_EV, ModeDriving
from grain_growth_pf.io.provenance import git_sha, software_versions
from grain_growth_pf.stochastic.hazard import CumulativeHazardClock
from grain_growth_pf.stochastic.multihit import poisson_completion_probability


def sampled_by_clock(rate: float, dt: float, seeds: np.ndarray) -> np.ndarray:
    samples = []
    for seed in seeds:
        clock = CumulativeHazardClock(np.random.default_rng(int(seed)))
        t = 0.0
        while True:
            events = clock.advance(rate, dt, t)
            if events:
                samples.append(events[0].event_time)
                break
            t += dt
    return np.asarray(samples)


def main() -> None:
    rng = np.random.default_rng(20260814)
    rate = 2.5; seeds = rng.integers(0, 2**32 - 1, 6000)
    dt_samples = {str(dt): sampled_by_clock(rate, dt, seeds) for dt in (0.08, 0.01)}
    single = dt_samples["0.08"]
    single_stats = {
        "input_rate": rate, "sample_mean": float(single.mean()), "expected_mean": 1 / rate,
        "sample_cv": float(single.std() / single.mean()),
        "ks_pvalue_exponential": float(kstest(single, "expon", args=(0, 1 / rate)).pvalue),
        "timestep_mean_relative_difference": float(abs(dt_samples["0.08"].mean() / dt_samples["0.01"].mean() - 1)),
        "paired_max_difference": float(np.max(np.abs(dt_samples["0.08"] - dt_samples["0.01"]))),
    }

    hits = 5
    erlang = rng.gamma(hits, 1 / rate, 100000)
    multihit_stats = {
        "K": hits, "sample_mean": float(erlang.mean()), "expected_mean": hits / rate,
        "sample_cv": float(erlang.std() / erlang.mean()), "expected_cv": 1 / np.sqrt(hits),
    }
    lam = 2.2
    packet_fraction = float(np.mean(rng.poisson(lam, 200000) >= 3))
    packet_stats = {"Lambda": lam, "K": 3, "sample_fraction": packet_fraction,
                    "exact_poisson_tail": poisson_completion_probability(3, lam)}

    r1, r2 = 2.0, 5.0
    parallel = np.minimum(rng.exponential(1 / r1, 100000), rng.exponential(1 / r2, 100000))
    serial = rng.exponential(1 / r1, 100000) + rng.exponential(1 / r2, 100000)
    composition = {
        "parallel_mean": float(parallel.mean()), "parallel_expected": 1 / (r1 + r2),
        "serial_mean": float(serial.mean()), "serial_expected": 1 / r1 + 1 / r2,
    }

    temperatures = np.array([700.0, 800.0, 900.0, 1050.0])
    barrier = 0.65
    mode = DisconnectionMode("activation", (0.2, 0), 0.2, 0, barrier, 1e7, 1)
    exact_rates = np.array([mode.rate(t, ModeDriving()) for t in temperatures])
    event_samples = [rng.exponential(1 / r, 20000) for r in exact_rates]
    measured_single_rates = np.array([1 / sample.mean() for sample in event_samples])
    multihit_samples = [rng.gamma(5, 1 / r, 20000) for r in exact_rates]
    measured_multihit_rates = np.array([5 / sample.mean() for sample in multihit_samples])
    single_fit = fit_activation_energy(temperatures, measured_single_rates)
    multi_fit = fit_activation_energy(temperatures, measured_multihit_rates)
    activation = {
        "temperatures": temperatures.tolist(), "input_barrier_ev": barrier,
        "exact_rates": exact_rates.tolist(), "measured_single_rates": measured_single_rates.tolist(),
        "measured_multihit_elementary_rates": measured_multihit_rates.tolist(),
        "single_hit_Q_ev": single_fit.activation_energy_ev,
        "single_hit_Q_standard_error_ev": single_fit.standard_error_ev,
        "persistent_K5_Q_ev": multi_fit.activation_energy_ev,
        "persistent_K5_Q_standard_error_ev": multi_fit.standard_error_ev,
    }
    result = {
        "git_sha": git_sha(), "software": software_versions(), "seed": 20260814,
        "single_hit": single_stats, "persistent_multihit": multihit_stats,
        "packet_reset": packet_stats, "parallel_serial": composition,
        "mechanism_isolation_activation": activation,
    }
    result["passed"] = bool(
        abs(single_stats["sample_mean"] * rate - 1) < 0.03
        and single_stats["ks_pvalue_exponential"] > 0.01
        and single_stats["paired_max_difference"] < 1e-10
        and abs(multihit_stats["sample_mean"] / (hits / rate) - 1) < 0.015
        and abs(multihit_stats["sample_cv"] - 1 / np.sqrt(hits)) < 0.01
        and abs(packet_fraction - packet_stats["exact_poisson_tail"]) < 0.005
        and abs(composition["parallel_mean"] / composition["parallel_expected"] - 1) < 0.015
        and abs(composition["serial_mean"] / composition["serial_expected"] - 1) < 0.015
        and abs(single_fit.activation_energy_ev - barrier) < 0.01
        and abs(multi_fit.activation_energy_ev - barrier) < 0.01
    )
    target = Path("results/validation/stochastic_validation.json")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["passed"]:
        raise SystemExit("stochastic validation failed")


if __name__ == "__main__":
    main()

