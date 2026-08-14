import numpy as np

from grain_growth_pf.config import PFConfig
from grain_growth_pf.pf.free_energy import free_energy
from grain_growth_pf.pf.geometry import circular_grain, equivalent_radius, planar_interface
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
    assert abs(slope / expected - 1) < 0.08
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
    assert abs(slopes[1] / slopes[0] - 2) < 0.06


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
