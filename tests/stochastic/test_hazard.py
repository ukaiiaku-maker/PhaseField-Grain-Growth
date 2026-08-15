import numpy as np
from scipy.stats import kstest

from grain_growth_pf.stochastic.hazard import CumulativeHazardClock, ParallelHazardClock
from grain_growth_pf.stochastic.multihit import MultiHitProcess, poisson_completion_probability


def sample_constant(rate, dt, count, seed=1):
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(count):
        clock = CumulativeHazardClock(rng)
        t = 0.0
        while True:
            events = clock.advance(rate, dt, t)
            if events:
                samples.append(events[0].event_time)
                break
            t += dt
    return np.asarray(samples)


def test_constant_rate_is_exponential_and_timestep_invariant():
    rate = 2.5
    coarse = sample_constant(rate, 0.08, 2500, 31)
    fine = sample_constant(rate, 0.01, 2500, 31)
    assert abs(coarse.mean() - 1 / rate) < 0.02
    assert abs(coarse.mean() - fine.mean()) < 0.02
    assert kstest(coarse, fine).pvalue > 0.01
    assert kstest(coarse, "expon", args=(0, 1 / rate)).pvalue > 0.01


def test_persistent_multihit_is_erlang():
    rng = np.random.default_rng(8)
    rate, hits = 3.0, 4
    values = []
    for _ in range(2500):
        process = MultiHitProcess(hits, rng)
        t = 0.0
        while True:
            complete = process.advance(rate, 0.02, t)
            if complete:
                values.append(complete[0].time)
                break
            t += 0.02
    values = np.asarray(values)
    assert abs(values.mean() - hits / rate) < 0.04
    assert abs(values.std() / values.mean() - 1 / np.sqrt(hits)) < 0.035


def test_packet_reset_poisson_tail():
    rng = np.random.default_rng(4)
    lam, hits, samples = 2.2, 3, 40000
    measured = np.mean(rng.poisson(lam, samples) >= hits)
    assert abs(measured - poisson_completion_probability(hits, lam)) < 0.01


def test_packet_and_persistent_windows_have_distinct_memory():
    persistent = MultiHitProcess(3, np.random.default_rng(10), "persistent_hits")
    packet = MultiHitProcess(3, np.random.default_rng(10), "packet_reset")
    persistent.hit_count = packet.hit_count = 2
    persistent.begin_window()
    packet.begin_window()
    assert persistent.hit_count == 2
    assert packet.hit_count == 0


def test_time_dependent_linear_hazard():
    rng = np.random.default_rng(88)
    a, dt = 0.7, 0.005
    sampled, exact = [], []
    for _ in range(500):
        clock = CumulativeHazardClock(rng)
        threshold = clock.threshold
        exact.append(np.sqrt(2 * threshold / a))
        t = 0.0
        previous = 0.0
        while True:
            rate = a * (t + dt)
            events = clock.advance(rate, dt, t, previous)
            if events:
                sampled.append(events[0].event_time)
                break
            previous, t = rate, t + dt
    assert np.max(np.abs(np.asarray(sampled) - exact)) < 0.006


def test_parallel_channel_probabilities():
    rng = np.random.default_rng(3)
    selected = []
    for _ in range(2500):
        clock = ParallelHazardClock(rng)
        t = 0.0
        while not (events := clock.advance([1.0, 3.0], 0.1, t)):
            t += 0.1
        event = events[0]
        selected.append(event.channel)
    assert abs(np.mean(np.asarray(selected) == 1) - 0.75) < 0.025
