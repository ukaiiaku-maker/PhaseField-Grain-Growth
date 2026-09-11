#!/usr/bin/env python3
"""HPC3 manufactured SHARP NETWORK gradient flows, not PF qualification."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from grain_growth_pf.mechanics.anisotropy import BoundaryLaw, LADDER
from grain_growth_pf.mechanics.cahn_hoffman import variation, herring_force


def evolve_loop(points, law, horizon, dt, orientations=(.1, .4)):
    points = np.asarray(points, float).copy()
    time = 0.
    records = []
    rejections = 0
    while time < horizon-1e-14:
        old = variation(points, law, *orientations, closed=True)
        theta = np.arctan2(old.normals[:, 1], old.normals[:, 0])
        mobility = law.evaluate(theta, *orientations)[4]
        velocity = mobility*old.normal_pressure
        dissipation = np.sum(old.dual_lengths*velocity**2/mobility)
        used_dt = min(dt, horizon-time)
        while True:
            proposed = points + used_dt*velocity[:, None]*old.normals
            new = variation(proposed, law, *orientations, closed=True)
            if new.energy <= old.energy + 1e-12:
                break
            used_dt /= 2
            rejections += 1
            if used_dt < 1e-10:
                raise RuntimeError('manufactured loop timestep exhausted')
        derivative = (new.energy-old.energy)/used_dt
        records.append([time+used_dt, new.energy, used_dt, dissipation,
                        abs(derivative+dissipation)/max(dissipation, 1e-14)])
        points, time = proposed, time+used_dt
    return points, np.asarray(records), rejections


def relax_junction(law, dt):
    orientations = [.1, .4, .7]
    angles = np.arange(3)*2*np.pi/3
    anchors = 3*np.column_stack((np.cos(angles), np.sin(angles)))
    pairs = [(0, 1), (1, 2), (0, 2)]
    junction = np.array([.3, -.2])
    def energy(x):
        return sum(variation([x, endpoint], law, orientations[i], orientations[j]).energy
                   for endpoint, (i, j) in zip(anchors, pairs))
    initial = energy(junction)
    previous = initial
    largest_increase = 0.
    for step in range(20000):
        force = herring_force(anchors-junction, pairs, orientations, law)
        if np.linalg.norm(force) < 1e-7:
            break
        candidate = junction+dt*force
        value = energy(candidate)
        largest_increase = max(largest_increase, value-previous)
        if value > previous+1e-10:
            raise RuntimeError('TJ relaxation energy increased')
        junction, previous = candidate, value
    else:
        raise RuntimeError('TJ residual did not converge')
    return {'junction': junction.tolist(), 'steps': step, 'dt': dt,
            'initial_energy': initial, 'final_energy': previous,
            'force_residual': float(np.linalg.norm(force)),
            'largest_energy_increase': largest_increase}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path)
    args = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURMD_NODENAME'):
        raise RuntimeError('Scientific evolution requires an HPC3 compute allocation')
    args.output.mkdir(parents=True, exist_ok=False)
    report = {'job_id': os.environ['SLURM_JOB_ID'], 'model': 'sharp polygon network',
              'pf_coupling_qualified': False, 'production_release': False, 'cases': {}}
    try:
        for name in ('A0_ISOTROPIC', 'A2_STRONG'):
            law = BoundaryLaw(LADDER[name])
            n = 128
            theta = np.arange(n)*2*np.pi/n
            points = 8*np.column_stack((np.cos(theta), np.sin(theta)))
            results = []
            for dt in (.002, .001):
                final, records, rejected = evolve_loop(points, law, .2, dt)
                np.savez_compressed(args.output/f'{name}-{dt}.npz', initial=points, final=final,
                                    time_energy_dt_dissipation_relative_balance_error=records)
                results.append({'dt': dt, 'steps': len(records), 'rejections': rejected,
                                'initial_energy': variation(points, law, .1, .4, closed=True).energy,
                                'final_energy': float(records[-1, 1]),
                                'maximum_relative_dissipation_error': float(records[:, 4].max()),
                                'mean_final_radius': float(np.linalg.norm(final, axis=1).mean())})
                if name == 'A0_ISOTROPIC':
                    expected = np.sqrt(8**2-2*4*.2)
                    results[-1]['circle_exact_radius_relative_error'] = abs(results[-1]['mean_final_radius']-expected)/expected
            report['cases'][name] = results
        report['junction'] = relax_junction(BoundaryLaw(), .01)
        report['junction_half_dt'] = relax_junction(BoundaryLaw(), .005)
        report['junction_refinement_distance'] = float(np.linalg.norm(
            np.array(report['junction']['junction'])-report['junction_half_dt']['junction']))
        report['sharp_network_tests_passed'] = (
            all(c['final_energy'] < c['initial_energy'] and c['maximum_relative_dissipation_error'] < .01
                for cases in report['cases'].values() for c in cases)
            and report['junction_refinement_distance'] < 1e-5)
        if not report['sharp_network_tests_passed']:
            raise RuntimeError('sharp-network energy or timestep gate failed')
    except Exception as exc:
        report['failure'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        (args.output/'manufactured.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
