#!/usr/bin/env python3
"""Assemble compact retrieved qualification evidence; never launches a solver."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('retrieved_output', type=Path)
    p.add_argument('destination', type=Path)
    p.add_argument('--earlier-output', type=Path, action='append', default=[])
    args = p.parse_args()
    src, dst = args.retrieved_output, args.destination
    dst.mkdir(parents=True, exist_ok=False)
    report = json.loads((src/'geometry/qualification.json').read_text()) if (src/'geometry/qualification.json').exists() else {}
    manufactured = json.loads((src/'manufactured/manufactured.json').read_text()) if (src/'manufactured/manufactured.json').exists() else {}
    suites = ET.parse(src/'tests.xml').getroot().iter('testsuite') if (src/'tests.xml').exists() else []
    tests = {k: 0 for k in ('tests', 'failures', 'errors', 'skipped')}
    for suite in suites:
        for key in tests:
            tests[key] += int(suite.get(key, '0'))
    test_exit = int((src/'test_exit_code.txt').read_text()) if (src/'test_exit_code.txt').exists() else None
    dump(dst/'qualification_decision.json', {
        'classification': 'ANISOTROPIC_IMPLEMENTATION_UNRESOLVED',
        'scope': 'geometry foundations and mathematical network checks',
        'production_release': False, 'pf_coupling_implemented': False,
        'junit_counts': tests, 'test_process_exit_code': test_exit,
        'full_suite_passed': test_exit == 0 and tests['tests'] > 0 and tests['failures'] == 0 and tests['errors'] == 0,
        'geometry_status': ('failed' if report.get('failure') else 'completed') if report else 'not_run',
        'manufactured_status': ('failed' if manufactured.get('failure') else 'completed') if manufactured else 'not_run',
        'geometry_failure': report.get('failure'),
        'manufactured_failure': manufactured.get('failure'),
        'sharp_network_tests_passed': manufactured.get('sharp_network_tests_passed', False),
        'effect_size_requirement': 'not evaluated',
        'minimum_model_change': 'undetermined',
        'production_trajectories_executed': 0,
    })
    dump(dst/'initial_network_normalization.json', {
        'status': 'PROVISIONAL_POLYGONAL_NOT_PRODUCTION_FROZEN' if report.get('ladder') else 'NOT_COMPUTED',
        'network': report.get('network'),
        'laws': {k: v['law'] for k, v in report.get('ladder', {}).items()},
        'initial_sha256': report.get('initial_sha256'),
    })
    dump(dst/'anisotropy_parameters.json', {
        'selected_production_strength': None, 'intended_strength': 'A2_STRONG',
        'A3_released': False, 'laws': {k: v['law'] for k, v in report.get('ladder', {}).items()},
        'hypothesis': 'synthetic inverse energy-mobility correlation',
    })
    with (dst/'constitutive_distributions.csv').open('w') as out:
        writer = csv.writer(out)
        writer.writerow(['strength', 'quantity', 'min', 'p05', 'median', 'p95', 'max', 'mean', 'ratio95_05'])
        for strength, data in report.get('ladder', {}).items():
            for quantity in ('gamma', 'mobility', 'stiffness', 'reduced_mobility'):
                row = data[quantity]
                writer.writerow([strength, quantity]+[row[k] for k in ('min', 'p05', 'median', 'p95', 'max', 'mean', 'ratio95_05')])
    prereg = json.loads((Path(__file__).resolve().parents[1]/'configs/production/anisotropic_phase1_preregistration.json').read_text())
    with (dst/'run_matrix.csv').open('w') as out:
        writer = csv.writer(out)
        writer.writerow(['trajectory_id', 'isotropic_counterpart', 'state', 'slurm_job_id'])
        for row in prereg['trajectories']:
            writer.writerow([row['trajectory_id'], row['isotropic_counterpart'], 'PREPARED_NOT_RELEASED', ''])
    with (dst/'isotropic_anisotropic_comparison.csv').open('w') as out:
        writer = csv.writer(out)
        writer.writerow(['trajectory_id', 'comparison_status'])
        writer.writerows([row['trajectory_id'], 'NOT_EVALUATED'] for row in prereg['trajectories'])
    for relative in ('tests.xml', 'tests.log', 'test_exit_code.txt', 'environment-fingerprint.json',
                     'pip-freeze.txt', 'numpy-build.txt', 'cpu.txt', 'runtime-libraries.txt',
                     'source_commit.txt', 'inputs.sha256', 'geometry/qualification.json',
                     'geometry/geometry.log', 'geometry.log', 'manufactured/manufactured.json',
                     'manufactured.log', 'geometry_exit_code.txt', 'manufactured_exit_code.txt'):
        original = src/relative
        if original.is_file():
            target = dst/'evidence'/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, target)
    figures = []
    if report.get('ladder') or manufactured.get('cases'):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    if report.get('ladder'):
        fig, axes = plt.subplots(1, 4, figsize=(12, 3.5), layout='constrained')
        strengths = list(report['ladder'])
        for ax, quantity in zip(axes, ('gamma', 'mobility', 'stiffness', 'reduced_mobility')):
            for i, strength in enumerate(strengths):
                stats = report['ladder'][strength][quantity]
                ax.plot([i, i], [stats['p05'], stats['p95']], lw=3)
                ax.plot(i, stats['median'], 'ko', ms=4)
            ax.set_xticks(range(len(strengths)), [s.split('_')[0] for s in strengths])
            ax.set_title(quantity.replace('_', ' '))
            ax.set_ylabel('Initial length-weighted p05–p95')
            ax.grid(axis='y', alpha=.2)
        fig.suptitle('Provisional polygonal initial network — no PF qualification')
        fig.savefig(dst/'constitutive_ranges.png', dpi=160)
        plt.close(fig)
        figures.append({'path': 'constitutive_ranges.png', 'scope': 'provisional initial-network constitutive summaries'})
    if manufactured.get('cases'):
        studies = []
        for previous in args.earlier_output:
            path = previous/'manufactured/manufactured.json'
            if path.exists():
                studies.append(json.loads(path.read_text()))
        studies.append(manufactured)
        dump(dst/'sharp_network_timestep_evidence.json', studies)
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), layout='constrained')
        for name in ('A0_ISOTROPIC', 'A2_STRONG'):
            points = sorted((c['dt'], c['maximum_relative_dissipation_error'], c['final_energy'])
                            for study in studies for c in study.get('cases', {}).get(name, []))
            if points:
                axes[0].loglog([v[0] for v in points], [v[1] for v in points], 'o-', label=name)
                axes[1].semilogx([v[0] for v in points], [v[2] for v in points], 'o-', label=name)
        axes[0].axhline(.01, color='black', ls='--', label='1% balance gate')
        axes[0].set_ylabel('Maximum relative dissipation-balance error')
        axes[1].set_ylabel('Final polygon energy at t=0.2')
        for ax in axes:
            ax.set_xlabel('Requested timestep')
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
        fig.suptitle('Sharp-network time refinement — no PF qualification')
        fig.savefig(dst/'sharp_network_timestep.png', dpi=160)
        plt.close(fig)
        figures.append({'path': 'sharp_network_timestep.png', 'scope': 'retrieved sharp-network time evidence; failed coarse tests retained'})
    dump(dst/'figure_manifest.json', {'figures': figures, 'production_figures': [], 'status': 'production not run'})
    (dst/'movie_index.csv').write_text('trajectory_id,movie_path,status\n'+''.join(
        f"{row['trajectory_id']},,NOT_PRODUCED_NO_TRAJECTORY\n" for row in prereg['trajectories']))
    dump(dst/'checksum_manifest.json', {str(path.relative_to(dst)): hashlib.sha256(path.read_bytes()).hexdigest()
                                      for path in sorted(dst.rglob('*')) if path.is_file()})


if __name__ == '__main__':
    main()
