import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.free_energy import free_energy
from grain_growth_pf.pf.geometry import circular_grain, equivalent_radius, planar_interface
from grain_growth_pf.pf.geometry import voronoi_polycrystal
from grain_growth_pf.pf.solver import MultiphaseFieldSolver


def test_planar_interface_stationary_and_constraint():
    cfg = PFConfig(shape=(48, 48), interface_width=4, time_step=0.05,
                   intrinsic_mobility=0.2, adaptive_stepping=True,
                   boundary_conditions="neumann")
    initial = planar_interface(cfg.shape, cfg.interface_width, angle=0)
    solver = MultiphaseFieldSolver(initial.copy(), cfg)
    solver.run(300)
    assert np.max(np.abs(solver.eta - initial)) < 2e-3
    assert np.max(np.abs(solver.eta.sum(axis=0) - 1)) < 1e-14


def test_circle_shrinks_parabolically():
    cfg = PFConfig(shape=(64, 64), interface_width=4, time_step=0.05,
                   intrinsic_mobility=0.2, gb_energy=1.0, adaptive_stepping=True)
    solver = MultiphaseFieldSolver(circular_grain(cfg.shape, 14, cfg.interface_width), cfg)
    times, radius2 = [], []
    for step in range(1800):
        solver.step()
        if step % 30 == 0:
            times.append(solver.time)
            radius2.append(equivalent_radius(solver.eta[1]) ** 2)
    slope, intercept = np.polyfit(times[8:-3], radius2[8:-3], 1)
    expected = -2 * cfg.intrinsic_mobility * cfg.gb_energy
    fitted = slope * np.asarray(times[8:-3]) + intercept
    r2 = 1 - np.sum((np.asarray(radius2[8:-3]) - fitted) ** 2) / np.sum((np.asarray(radius2[8:-3]) - np.mean(radius2[8:-3])) ** 2)
    # Four cells across a compact obstacle interface is the coarsest supported
    # resolution; the mesh-validation campaign requires convergence below it.
    assert abs(slope / expected - 1) < 0.12
    assert r2 > 0.999


def test_mobility_scaling_and_energy_decrease():
    slopes = []
    for mobility in (0.1, 0.2):
        cfg = PFConfig(shape=(48, 48), interface_width=4, time_step=0.05,
                       intrinsic_mobility=mobility, adaptive_stepping=True)
        solver = MultiphaseFieldSolver(circular_grain(cfg.shape, 11, 4), cfg)
        ts, r2s, energies = [], [], []
        for step in range(600):
            diag = solver.step()
            if step % 20 == 0:
                ts.append(solver.time)
                r2s.append(equivalent_radius(solver.eta[1]) ** 2)
                energies.append(diag.interfacial_energy)
        slopes.append(np.polyfit(ts[5:], r2s[5:], 1)[0])
        assert np.all(np.diff(energies) <= 1e-9)
    assert abs(slopes[1] / slopes[0] - 2) < 0.08


def test_high_mobility_time_rescaling_remains_quantitative():
    cfg = PFConfig(shape=(64, 64), interface_width=4, time_step=0.04,
                   intrinsic_mobility=4.0, adaptive_stepping=True)
    solver = MultiphaseFieldSolver(circular_grain(cfg.shape, 15, 4), cfg)
    times, radius2 = [], []
    for step in range(250):
        solver.step()
        if step % 10 == 0:
            times.append(solver.time)
            radius2.append(equivalent_radius(solver.eta[1]) ** 2)
    slope, intercept = np.polyfit(times[4:-3], radius2[4:-3], 1)
    fitted = slope * np.asarray(times[4:-3]) + intercept
    r_squared = 1 - np.sum((np.asarray(radius2[4:-3]) - fitted) ** 2) / np.sum(
        (np.asarray(radius2[4:-3]) - np.mean(radius2[4:-3])) ** 2
    )
    assert abs(slope / (-2 * cfg.intrinsic_mobility * cfg.gb_energy) - 1) < 0.12
    assert r_squared > 0.9999


def test_restart_is_exact():
    cfg = PFConfig(shape=(32, 32), interface_width=4, time_step=0.04, intrinsic_mobility=0.2)
    eta = circular_grain(cfg.shape, 8, 4)
    continuous = MultiphaseFieldSolver(eta.copy(), cfg)
    continuous.run(30)
    restarted = MultiphaseFieldSolver(eta.copy(), cfg)
    restarted.run(12)
    state = restarted.state_dict()
    restored = MultiphaseFieldSolver(eta.copy(), cfg)
    restored.load_state_dict(state)
    restored.run(18)
    assert np.array_equal(continuous.eta, restored.eta)
    assert continuous.time == restored.time


def test_extinct_grains_cannot_resurrect():
    eta, _, _ = voronoi_polycrystal((32, 32), 18, seed=10, width=2)
    cfg = PFConfig(shape=(32, 32), interface_width=3, time_step=0.04,
                   intrinsic_mobility=4.0, adaptive_stepping=True)
    solver = MultiphaseFieldSolver(eta, cfg)
    counts = []
    for step in range(300):
        solver.step()
        if step % 10 == 0:
            counts.append(np.count_nonzero(solver.active_phases))
            assert set(np.unique(solver.labels)).issubset(set(np.flatnonzero(solver.active_phases)))
    assert np.all(np.diff(counts) <= 0)
    assert counts[-1] < counts[0]


def test_active_phase_can_advance_one_cell_but_not_nucleate_remotely():
    eta = np.zeros((3, 8, 8))
    eta[0, :, :4] = 1.0
    eta[1, :, 4:6] = 1.0
    eta[2, :, 6:] = 1.0
    eta[0, :, 5] = 0.5
    eta[1, :, 5] = 0.5
    cfg = PFConfig(shape=(8, 8), interface_width=3, time_step=0.005,
                   intrinsic_mobility=0.2, boundary_conditions="neumann")
    solver = MultiphaseFieldSolver(eta, cfg)
    solver.step()
    assert solver.eta[2, 4, 5] > 0
    assert solver.eta[2, 4, 0] == 0
