import csv
import json

import numpy as np

from grain_growth_pf.config import ModelConfig, PFConfig
from grain_growth_pf.disconnections.mode import ModeDriving
from grain_growth_pf.pf.initial_conditions import prepare_initial_condition
from grain_growth_pf.simulation import EventResolvedSimulation
from grain_growth_pf.stochastic.multihit import MultiHitProcess, poisson_completion_probability


def test_event_resolved_smoke_writes_reproducible_schema(tmp_path):
    config = ModelConfig(
        regime="G1", seed=17,
        pf=PFConfig(shape=(24, 24), interface_width=3, time_step=0.02,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "single_hit_poisson", "shear_memory", "free_volume"),
        output_cadence=2, max_steps=6, termination_grains=1,
        parameters={"initial_grains": 6, "equilibration_steps": 0, "encounter_density": 5.0,
                    "attempt_frequency": 100.0, "event_domain_length": 100.0},
    )
    output = tmp_path / "run"
    EventResolvedSimulation(config, output, code_sha="captured-launch-sha").run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["config"]["seed"] == 17
    assert manifest["git_sha"] == "captured-launch-sha"
    assert {item["path"].split("/")[-1] for item in manifest["restart_artifacts"]} == {
        "checkpoint.npz", "checkpoint.json"
    }
    assert all(len(item["sha256"]) == 64 for item in manifest["restart_artifacts"])
    with (output / "grain_tracks.csv").open() as handle:
        tracks = list(csv.DictReader(handle))
    with (output / "events.csv").open() as handle:
        fields = csv.DictReader(handle).fieldnames
    assert tracks
    assert "random_hazard_threshold" in fields
    assert "burgers_vector_b" in fields


def test_qiu_full_field_backend_smoke(tmp_path):
    config = ModelConfig(
        regime="E1", seed=9,
        pf=PFConfig(shape=(20, 20), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="qiu_full_field", compatibility_model="explicit_modes",
        active_modules=("event_modes", "shear_feedback"), output_cadence=1,
        max_steps=2, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 0,
                    "attempt_frequency": 100.0, "barrier_core_ev": 0.0,
                    "b_coefficient_ev": 0.0, "h_coefficient_ev": 0.0,
                    "event_domain_length": 100.0},
    )
    output = tmp_path / "qiu"
    simulation = EventResolvedSimulation(config, output)
    simulation.run()
    assert simulation.full_field is not None
    assert simulation.full_field.stress.shape == (2, 2, 20, 20)
    assert np.any(simulation.full_field.eigenstrain != 0)
    assert np.any(simulation.full_field.stress != 0)
    assert np.any(simulation.driving_field != 0)
    with (tmp_path / "qiu" / "events.csv").open() as handle:
        strain_sum = sum(float(row["shear_strain_increment"]) for row in csv.DictReader(handle))
    assert np.isclose(strain_sum, simulation.accumulated_shear_strain)


def test_simulation_wires_quenched_barrier_distribution(tmp_path):
    config = ModelConfig(
        regime="barrier-disorder", seed=91,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={
            "initial_grains": 5,
            "barrier_distribution": "truncated_gaussian",
            "barrier_mean_ev": 0.55,
            "barrier_std_ev": 0.2,
            "barrier_bounds_ev": [0.4, 0.7],
        },
    )
    simulation = EventResolvedSimulation(config, tmp_path / "barrier-disorder")
    barriers = np.asarray([mode.barrier_ev for mode in simulation.modes])
    assert np.all((barriers >= 0.4) & (barriers <= 0.7))
    assert np.std(barriers) > 0
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_vectorized_mode_rates_match_individual_mode_equations(tmp_path):
    config = ModelConfig(
        regime="rate-equivalence", seed=191,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    temperature=875.0),
        compatibility_model="explicit_modes", active_modules=("event_modes",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "event_domain_length": 100.0},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "rate-equivalence")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    domain.shear.state = 0.17
    domain.free_volume.required_total = 0.31
    candidates, rates, normal, shear, vacancy = simulation._activation_rates(domain, segment)
    expected = np.asarray([
        mode.rate(config.pf.temperature, ModeDriving(normal, shear[index], vacancy))
        for index, mode in enumerate(candidates)
    ])
    assert np.allclose(rates, expected, rtol=3e-14, atol=0.0)
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_production_climb_rates_use_butler_volmer_and_l_squared_transport(tmp_path):
    config = ModelConfig(
        regime="climb-rates", seed=192,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    temperature=900.0),
        active_modules=("free_volume", "serial_climb"),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={
            "initial_grains": 5, "event_domain_length": 100.0,
            "exchange_prefactor": 1e4, "exchange_barrier_ev": 0.4,
            "transport_prefactor": 1e4, "transport_barrier_ev": 0.5,
            "free_volume_stiffness": 0.1,
        },
    )
    simulation = EventResolvedSimulation(config, tmp_path / "climb-rates")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    domain.free_volume.required_total = 0.2
    _, exchange_1, transport_1 = simulation._stage_rates(domain, segment)
    domain.free_volume.required_total = 0.4
    _, exchange_2, _ = simulation._stage_rates(domain, segment)
    segment.length *= 2.0
    _, _, transport_2 = simulation._stage_rates(domain, segment)
    assert exchange_2 > exchange_1
    assert np.isclose(transport_2, transport_1 / 4.0)
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_packet_renewal_window_matches_poisson_tail_and_checkpoints_age(tmp_path):
    config = ModelConfig(
        regime="packet-renewal", seed=193,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "multihit_packet_reset"),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "required_hits": 3,
                    "packet_window_time": 1.0},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "packet-tail")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    rng = np.random.default_rng(194)
    completed = 0
    samples, window_hazard = 12000, 2.2
    for _ in range(samples):
        domain.activation = MultiHitProcess(3, rng, "packet_reset")
        domain.packet_window_elapsed = 0.0
        completions, _ = simulation._advance_activation(domain, window_hazard, 1.0, 0.0)
        completed += bool(completions)
    expected = poisson_completion_probability(3, window_hazard)
    assert abs(completed / samples - expected) < 0.015
    domain.packet_window_elapsed = 0.37
    restored = simulation._new_domain(segment)
    restored.load_state_dict(domain.state_dict())
    assert restored.packet_window_elapsed == 0.37
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_packet_reset_cannot_accumulate_hits_across_renewal_windows(tmp_path):
    config = ModelConfig(
        regime="packet-memory", seed=195,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        active_modules=("multihit_packet_reset",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "required_hits": 2,
                    "packet_window_time": 0.5},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "packet-memory")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    packet = simulation.domains[segment.entity_id]
    persistent = simulation._new_domain(segment)
    persistent.activation = MultiHitProcess(2, np.random.default_rng(196), "persistent_hits")
    for domain in (packet, persistent):
        domain.activation.clock.cumulative_hazard = 0.0
        domain.activation.clock.threshold = 0.25
        domain.activation.clock.last_rate = None
        first, _ = simulation._advance_activation(domain, 1.0, 0.251, 0.0)
        assert not first and domain.activation.hit_count == 1
        domain.activation.clock.threshold = 100.0
        simulation._advance_activation(domain, 1.0, 0.249, 0.251)
    assert packet.activation.hit_count == 0
    assert persistent.activation.hit_count == 1
    for domain in (packet, persistent):
        domain.activation.clock.threshold = domain.activation.clock.cumulative_hazard + 0.25
    packet_completion, _ = simulation._advance_activation(packet, 1.0, 0.251, 0.5)
    persistent_completion, _ = simulation._advance_activation(persistent, 1.0, 0.251, 0.5)
    assert not packet_completion
    assert persistent_completion
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_packet_release_stops_window_at_completion_time(tmp_path):
    config = ModelConfig(
        regime="packet-stop", seed=197,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        active_modules=("multihit_packet_reset",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "required_hits": 1,
                    "packet_window_time": 1.0},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "packet-stop")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    domain.activation.clock.cumulative_hazard = 0.0
    domain.activation.clock.threshold = 0.25
    domain.activation.clock.last_rate = None
    completions, hits = simulation._advance_activation(
        domain, rate=10.0, dt=1.0, start_time=4.0,
        stop_after_completion=True,
    )
    assert len(completions) == len(hits) == 1
    assert np.isclose(completions[0].time, 4.025)
    assert np.isclose(domain.packet_window_elapsed, 0.025)
    assert domain.activation.clock.cumulative_hazard == 0.25
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_single_hit_blocked_gb_ledger_has_one_hit_per_release(tmp_path):
    config = ModelConfig(
        regime="single-hit-ledger", seed=198,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "single_hit_poisson"),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={
            "initial_grains": 5, "required_hits": 1,
            "attempt_frequency": 1e6, "barrier_core_ev": 0.0,
            "b_coefficient_ev": 0.0, "h_coefficient_ev": 0.0,
        },
    )
    output = tmp_path / "single-hit-ledger"
    simulation = EventResolvedSimulation(config, output)
    simulation.solver.step_number = 1
    simulation.solver.time = config.pf.time_step
    for domain in simulation.domains.values():
        domain.blocked = True
        simulation._begin_activation_window(domain)
        domain.activation.clock.threshold = 1e-9
    simulation._update_physics()
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()
    with (output / "events.csv").open() as handle:
        event_types = [row["event_type"] for row in csv.DictReader(handle)]
    hits = event_types.count("activation_hit")
    releases = event_types.count("compatibility_release")
    assert hits == releases
    assert releases > 0


def test_accumulated_atomic_step_triggers_finite_pf_release(tmp_path):
    config = ModelConfig(
        regime="subgrid-release", seed=92,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01),
        compatibility_model="explicit_modes", active_modules=("event_modes",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={
            "initial_grains": 5, "easy_beta": 0.0,
            "pf_release_displacement": 0.25, "event_normal_pressure": 1.0,
        },
    )
    simulation = EventResolvedSimulation(config, tmp_path / "subgrid-release")
    segment = next(iter(simulation.snapshot.boundaries.values()))
    domain = simulation.domains[segment.entity_id]
    mode = next(mode for mode in simulation.modes if mode.step_height > 0)
    simulation._record_event(
        domain, segment, mode, 1.0, ModeDriving(), "test-release", 0.0
    )
    assert np.isclose(domain.normal_displacement_ledger, mode.step_height)
    simulation._update_physics()
    assert domain.normal_displacement_ledger == 0.0
    assert np.isclose(domain.normal_release_remaining, 0.25)
    assert np.any(simulation.driving_field != 0.0)
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_continuous_qiu_reference_converts_boundary_motion_to_eigenstrain(tmp_path):
    config = ModelConfig(
        regime="Q1", seed=19,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="qiu_full_field", active_modules=("qiu_reference_shear",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "easy_beta": 0.5},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "qiu-continuous")
    for segment in simulation.snapshot.boundaries.values():
        domain = simulation.domains[segment.entity_id]
        domain.previous_area_i -= 1.0
        domain.previous_area_j += 1.0
        domain.previous_time = -config.pf.time_step
    simulation.solver.step_number = 1
    simulation._update_physics()
    assert simulation.full_field is not None
    assert np.any(simulation.full_field.eigenstrain != 0)
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()


def test_mixed_event_releases_shear_and_free_volume_together(tmp_path):
    config = ModelConfig(
        regime="SC2", seed=23,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        mechanics_backend="local_memory", compatibility_model="explicit_modes",
        active_modules=("mixed_shear_climb_event",), output_cadence=1,
        max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "attempt_frequency": 100.0,
                    "barrier_core_ev": 0.0, "b_coefficient_ev": 0.0,
                    "h_coefficient_ev": 0.0, "shear_trigger": 0.1,
                    "climb_trigger_quota": 0.1},
    )
    simulation = EventResolvedSimulation(config, tmp_path / "mixed")
    for domain in simulation.domains.values():
        domain.shear.state = 1.0
        domain.free_volume.required_total = 1.0
    simulation.solver.step_number = 1
    simulation.solver.time = config.pf.time_step
    simulation._update_physics()
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()
    with (tmp_path / "mixed" / "events.csv").open() as handle:
        events = list(csv.DictReader(handle))
    assert events
    assert all(row["barrier_type"] != "easy" for row in events)
    assert any(float(row["release_Delta_s"]) > 0 and float(row["release_Delta_q"]) > 0
               for row in events)


def test_event_simulation_checkpoint_restart_is_exact(tmp_path):
    config = ModelConfig(
        regime="G2", seed=44,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "multihit_persistent", "shear_memory", "free_volume"),
        output_cadence=1, max_steps=4, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 0,
                    "encounter_density": 2.0, "attempt_frequency": 2.0,
                    "event_domain_length": 100.0},
    )
    continuous = EventResolvedSimulation(config, tmp_path / "continuous")
    continuous.run()

    interrupted = EventResolvedSimulation(config, tmp_path / "interrupted")
    for _ in range(2):
        interrupted.solver.step()
        interrupted.snapshot = interrupted.tracker.update(interrupted.solver.labels)
        interrupted._update_physics()
    interrupted._save_checkpoint()
    interrupted.ledger.close(); interrupted.track_handle.close(); interrupted.boundary_handle.close()

    resumed = EventResolvedSimulation(config, tmp_path / "interrupted", resume=True)
    resumed.run()
    assert resumed.solver.step_number == continuous.solver.step_number
    assert (resumed.solver.eta == continuous.solver.eta).all()
    assert {k: v.state_dict() for k, v in resumed.domains.items()} == {
        k: v.state_dict() for k, v in continuous.domains.items()
    }


def test_checkpoint_archive_is_authoritative_if_metadata_replacement_is_interrupted(tmp_path):
    config = ModelConfig(
        regime="G1", seed=45,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="geometric_surrogate",
        active_modules=("gb_compatibility", "single_hit_poisson"),
        output_cadence=1, max_steps=3, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 0,
                    "encounter_density": 2.0, "attempt_frequency": 2.0,
                    "event_domain_length": 100.0},
    )
    output = tmp_path / "interrupted-metadata"
    simulation = EventResolvedSimulation(config, output)
    simulation.solver.step()
    simulation.snapshot = simulation.tracker.update(simulation.solver.labels)
    simulation._update_physics()
    simulation._save_checkpoint()
    saved_eta = simulation.solver.eta.copy()
    saved_step = simulation.solver.step_number
    simulation.ledger.close(); simulation.track_handle.close(); simulation.boundary_handle.close()

    # Model a process interruption after the archive replacement but before the
    # companion human-readable metadata replacement.
    (output / "checkpoint.json").write_text('{"step_number": -1}\n')
    resumed = EventResolvedSimulation(config, output, resume=True)
    assert resumed.solver.step_number == saved_step
    assert np.array_equal(resumed.solver.eta, saved_eta)
    resumed.ledger.close(); resumed.track_handle.close(); resumed.boundary_handle.close()


def test_named_temperature_and_tj_particle_regimes_are_distinct(tmp_path):
    common = dict(shape=(18, 18), interface_width=3, time_step=0.01,
                  intrinsic_mobility=1.0, adaptive_stepping=True)
    low = ModelConfig(regime="B1", seed=3, pf=PFConfig(**common, temperature=700),
                      active_modules=("arrhenius_intrinsic",), max_steps=1,
                      parameters={"initial_grains": 5, "intrinsic_barrier_ev": 0.4})
    high = ModelConfig(regime="B1", seed=3, pf=PFConfig(**common, temperature=1050),
                       active_modules=("arrhenius_intrinsic",), max_steps=1,
                       parameters={"initial_grains": 5, "intrinsic_barrier_ev": 0.4})
    low_sim = EventResolvedSimulation(low, tmp_path / "low")
    high_sim = EventResolvedSimulation(high, tmp_path / "high")
    assert high_sim.solver.config.intrinsic_mobility > low_sim.solver.config.intrinsic_mobility
    low_sim.ledger.close(); low_sim.track_handle.close()
    low_sim.boundary_handle.close()
    high_sim.ledger.close(); high_sim.track_handle.close(); high_sim.boundary_handle.close()

    tj_config = ModelConfig(
        regime="P5", seed=12, pf=PFConfig(**common, temperature=900),
        compatibility_model="off", active_modules=("gb_pinning", "tj_pinning", "random_spatial_pinning"),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 7, "equilibration_steps": 0, "particle_count": 100,
                    "particle_radius": 3.0, "encounter_density": 10.0},
    )
    simulation = EventResolvedSimulation(tj_config, tmp_path / "pins")
    simulation.run()
    assert simulation.particles is not None
    assert simulation.tj_domains
    assert (simulation.solver.mobility_scale == 0).any()


def test_equilibration_precedes_physical_time_and_hazard(tmp_path):
    config = ModelConfig(
        regime="G1", seed=5,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01, intrinsic_mobility=0.2),
        compatibility_model="geometric_surrogate", active_modules=("gb_compatibility",),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 5, "equilibration_steps": 3, "encounter_density": 100.0},
    )
    output = tmp_path / "equilibrated"
    simulation = EventResolvedSimulation(config, output)
    assert simulation.solver.time == 0
    assert simulation.solver.step_number == 0
    simulation.run()
    with (output / "grain_tracks.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert any(float(row["time"]) == 0 for row in rows)
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["equilibration_steps_completed_before_time_zero"] == 3


def test_target_equilibration_compacts_phases_before_time_zero(tmp_path):
    config = ModelConfig(
        regime="B0", seed=1,
        pf=PFConfig(shape=(24, 24), interface_width=3, time_step=0.04,
                    intrinsic_mobility=2.0, adaptive_stepping=True),
        output_cadence=1, max_steps=1, termination_grains=1,
        parameters={"initial_grains": 8, "equilibration_steps": 0,
                    "equilibrate_to_grains": 7, "equilibration_max_steps": 500},
    )
    output = tmp_path / "target-equilibrated"
    simulation = EventResolvedSimulation(config, output)
    assert simulation.solver.time == 0
    assert simulation.solver.eta.shape[0] == 7
    assert len(simulation.orientations) == 7
    simulation.run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["grains_after_equilibration"] == 7
    assert manifest["equilibration_steps_completed_before_time_zero"] > 0
    resumed = EventResolvedSimulation(config, output, resume=True)
    resumed.run()
    assert resumed.solver.eta.shape[0] == 7
    assert len(resumed.orientations) == 7


def test_cached_initial_condition_is_loaded_without_reaging(tmp_path):
    pf = PFConfig(shape=(24, 24), interface_width=3, time_step=0.04,
                  intrinsic_mobility=2.0, adaptive_stepping=True)
    parameters = {"initial_grains": 8, "equilibrate_to_grains": 7,
                  "equilibration_max_steps": 500}
    state_path = prepare_initial_condition(pf, 1, parameters, tmp_path / "initial.npz", "test-sha")
    config = ModelConfig(
        regime="B0", seed=1, pf=pf, output_cadence=1, max_steps=1, termination_grains=1,
        parameters={**parameters, "initial_state_file": str(state_path)},
    )
    output = tmp_path / "cached-run"
    simulation = EventResolvedSimulation(config, output)
    assert simulation.solver.step_number == 0
    assert simulation.solver.eta.shape[0] == 7
    simulation.run()
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["initial_condition_source"] == str(state_path)
    assert manifest["grains_after_equilibration"] == 7


def test_output_cadence_does_not_change_stochastic_trajectory(tmp_path):
    common = dict(
        regime="E0", seed=88,
        pf=PFConfig(shape=(18, 18), interface_width=3, time_step=0.01,
                    intrinsic_mobility=0.1, adaptive_stepping=True),
        compatibility_model="explicit_modes", active_modules=("event_modes",),
        max_steps=5, termination_grains=1,
        parameters={"initial_grains": 5, "barrier_core_ev": 0.0,
                    "b_coefficient_ev": 0.0, "h_coefficient_ev": 0.0,
                    "attempt_frequency": 10.0},
    )
    frequent = EventResolvedSimulation(ModelConfig(**common, output_cadence=1), tmp_path / "frequent")
    sparse = EventResolvedSimulation(ModelConfig(**common, output_cadence=5), tmp_path / "sparse")
    frequent.run(); sparse.run()
    assert np.array_equal(frequent.solver.eta, sparse.solver.eta)
    with (tmp_path / "frequent" / "events.csv").open() as handle:
        frequent_events = [{k: v for k, v in row.items() if k != "run_id"} for row in csv.DictReader(handle)]
    with (tmp_path / "sparse" / "events.csv").open() as handle:
        sparse_events = [{k: v for k, v in row.items() if k != "run_id"} for row in csv.DictReader(handle)]
    assert frequent_events == sparse_events
