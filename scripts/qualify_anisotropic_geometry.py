#!/usr/bin/env python3
"""HPC3-only initial-network and constitutive qualification (no PF evolution)."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.stats import rankdata

from grain_growth_pf.entities.conforming_network import reconstruct_periodic
from grain_growth_pf.mechanics.anisotropy import BoundaryLaw, LADDER
from grain_growth_pf.mechanics.cahn_hoffman import variation


NPZ_SHA = '106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6'
META_SHA = 'ff66cb5f938b5b70fe9e611674956f1e45f30f0979714507a853750f61e1a9bb'


def weighted_quantile(values, weights, probabilities):
    order = np.argsort(values, kind='stable')
    values, weights = np.asarray(values)[order], np.asarray(weights)[order]
    cdf = (np.cumsum(weights)-weights/2)/weights.sum()
    return np.interp(probabilities, cdf, values)


def describe(values, weights):
    q = weighted_quantile(values, weights, [.05, .5, .95])
    return dict(min=float(np.min(values)), p05=float(q[0]), median=float(q[1]),
                p95=float(q[2]), max=float(np.max(values)),
                mean=float(np.average(values, weights=weights)), ratio95_05=float(q[2]/q[0]))


def correlation(a, b, weights):
    a = a-np.average(a, weights=weights)
    b = b-np.average(b, weights=weights)
    denominator = np.sqrt(np.sum(weights*a*a)*np.sum(weights*b*b))
    return float(np.sum(weights*a*b)/denominator) if denominator > 1e-24 else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('initial', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not os.environ.get('SLURMD_NODENAME'):
        raise RuntimeError('Scientific qualification requires an HPC3 compute allocation')
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    for path, expected in [(args.initial, NPZ_SHA), (args.initial.with_suffix('.json'), META_SHA)]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'initial input checksum mismatch: {path}')
    report = {'classification': 'ANISOTROPIC_IMPLEMENTATION_UNRESOLVED',
              'job_id': os.environ['SLURM_JOB_ID'], 'initial_sha256': NPZ_SHA,
              'initial_metadata_sha256': META_SHA, 'matched_isotropic_initial_state': True,
              'pf_coupling_qualified': False, 'production_release': False,
              'scientific_evolution_performed': False}
    try:
        with np.load(args.initial, allow_pickle=False) as data:
            eta, orientations = data['eta'], data['orientations']
        assert eta.shape == (800, 384, 384)
        assert orientations.shape == (800,)
        assert len(np.unique(np.argmax(eta, axis=0))) == 800
        vertices, edges, pairs, paths = reconstruct_periodic(eta)
        del eta
        box = np.array([384., 384.])
        delta = vertices[edges[:, 1]]-vertices[edges[:, 0]]
        delta -= np.round(delta/box)*box
        lengths = np.linalg.norm(delta, axis=1)
        theta = np.arctan2(-delta[:, 0], delta[:, 1])
        phi_i, phi_j = orientations[pairs[:, 0]], orientations[pairs[:, 1]]
        np.savez_compressed(output/'initial_network.npz', vertices=vertices, edges=edges,
                            pairs=pairs, lengths=lengths, theta=theta, orientations=orientations)
        report['network'] = {'vertices': len(vertices), 'edges': len(edges),
                             'polylines': len(paths), 'closed_loops': sum(p.closed for p in paths),
                             'total_length': float(lengths.sum()),
                             'unique_pairs': len(np.unique(pairs, axis=0)),
                             'reconstruction': 'piecewise-linear upper envelope; fixed cell diagonal',
                             'continuum_geometry_converged': False}
        report['ladder'] = {}
        for name in ('A0_ISOTROPIC', 'A1_MODERATE', 'A2_STRONG'):
            law = BoundaryLaw(LADDER[name]).normalize(theta, phi_i, phi_j, lengths)
            gamma, first, second, stiffness, mobility = law.evaluate(theta, phi_i, phi_j)
            stats = {key: describe(value, lengths) for key, value in
                     [('gamma', gamma), ('mobility', mobility), ('stiffness', stiffness),
                      ('reduced_mobility', mobility*stiffness)]}
            dense_min = float('inf')
            dense_max = 0.
            angles = np.arange(8192)*(np.pi/2/8192)
            for i, j in np.unique(pairs, axis=0):
                dense = law.evaluate(angles, orientations[i], orientations[j])[3]
                dense_min, dense_max = min(dense_min, float(dense.min())), max(dense_max, float(dense.max()))
            # Fundamental-zone cross-product; orientation objectivity removes one angle.
            for misorientation in np.linspace(0, np.pi/4, 257):
                dense = law.evaluate(angles, 0., misorientation)[3]
                dense_min = min(dense_min, float(dense.min()))
            gates = {'energy_range': stats['gamma']['ratio95_05'] >= 1.25,
                     'mobility_range': stats['mobility']['ratio95_05'] >= 2.,
                     'reduced_mobility_range': stats['reduced_mobility']['ratio95_05'] >= 1.75,
                     'stiffness_margin': dense_min >= .08}
            stats.update(law=asdict(law), dense_min_stiffness=dense_min,
                         dense_actual_pair_max_stiffness=dense_max,
                         pearson=correlation(gamma, mobility, lengths),
                         rank_correlation=correlation(rankdata(gamma), rankdata(mobility), lengths),
                         constitutive_gates=gates, all_constitutive_gates=all(gates.values()))
            # Sum energy from ordered polylines independently of the raw edge inventory.
            energy = sum(variation(p.points, law, orientations[p.pair[0]], orientations[p.pair[1]],
                                   closed=p.closed, closing_displacement=p.winding if p.closed else None).energy
                         for p in paths)
            stats['edge_polyline_energy_relative_error'] = abs(energy-np.dot(gamma, lengths))/energy
            if stats['edge_polyline_energy_relative_error'] > 1e-10:
                raise RuntimeError('edge partition changed total energy')
            report['ladder'][name] = stats
            (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        report['A3_status'] = 'NOT_RELEASED: A2 numerical qualification must pass before escalation'
        report['normalization_status'] = 'PROVISIONAL_POLYGONAL: reconstruction convergence pending'
    except Exception as exc:
        report['failure'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        report['elapsed_seconds'] = time.monotonic()-start
        (output/'qualification.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
