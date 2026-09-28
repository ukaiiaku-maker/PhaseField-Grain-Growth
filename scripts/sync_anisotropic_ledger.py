#!/usr/bin/env python3
"""Idempotently reconcile dedicated local runner records; no SSH or submission."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path


def reconcile(root):
    root = Path(root).resolve()
    path = root/'ownership.json'
    with (root/'ownership.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        original = path.read_text()
        ledger = json.loads(original)
        jobs = ledger['jobs']
        by_id = {j['runner_run_id']: j for j in jobs}
        if len(by_id) != len(jobs):
            raise RuntimeError('duplicate runner ownership records')
        known_job_ids = [j['slurm_job_id'] for j in jobs if j.get('slurm_job_id')]
        if len(set(known_job_ids)) != len(known_job_ids):
            raise RuntimeError('duplicate Slurm ownership records')
        for item in sorted((root/'.hpc3/runs').glob('*.json')):
            record = json.loads(item.read_text())
            if record['project'] != 'pfgg-anisotropic-cahn-hoffman-v1':
                raise RuntimeError('refusing a different project')
            run_id = record['run_id']
            if run_id not in by_id:
                raise RuntimeError(f'unclaimed runner record: {run_id}')
            job = by_id[run_id]
            if record['input_checksums']['source.tar.gz'] != job['source_archive_sha256']:
                raise RuntimeError('source identity changed')
            previous_id = job.get('slurm_job_id')
            if previous_id and previous_id != record.get('job_id'):
                raise RuntimeError('one prepared run acquired two Slurm job IDs')
            job['slurm_job_id'] = record.get('job_id')
            job['scheduler_state'] = record['state']
            job['runner_fetch_status'] = record['fetch_status']
            job['runner_state_history'] = record['state_history']
            job['submission_timestamp'] = record.get('submission_timestamp')
            job['runner_manifest_sha256'] = hashlib.sha256(json.dumps(
                record['manifest'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        new_job_ids = [j['slurm_job_id'] for j in jobs if j.get('slurm_job_id')]
        if len(set(new_job_ids)) != len(new_job_ids):
            raise RuntimeError('two calculations claim one Slurm job')
        updated = json.dumps(ledger, indent=2)+'\n'
        if updated == original:
            return False
        temporary = path.with_name(path.name+f'.tmp.{os.getpid()}')
        with temporary.open('x') as handle:
            handle.write(updated)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        descriptor = os.open(root, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('state_root', type=Path)
    args = parser.parse_args()
    print('updated' if reconcile(args.state_root) else 'unchanged')
