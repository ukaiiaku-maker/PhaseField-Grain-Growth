# Anisotropic HPC3 execution evidence

Classification: **ANISOTROPIC_IMPLEMENTATION_UNRESOLVED**. These are environment, regression and mathematical network jobs, not production PF trajectories. No reduced arrays or array task IDs exist; each job is one task. No full production or continuation IDs exist.

All submitted jobs use **SDILLON1 / standard**, one node and one task, no GPU, no high QOS, and one numerical-library thread. Requests respect the 6 GiB/core partition limit. The latest preflight recorded 547 SU available, remote usage 6.42 GiB/1 TiB and 43.75k/8m inodes, and local free space 197 GiB. The other account SDILLON1_LAB was not used by this session.

Total maximum requested exposure: **7 CPU-hours**. Recorded actual CPU use: **0.134444 CPU-hours**. Recorded allocated core-walltime: **1.065556 hours**. All five submitted jobs are terminal and have accounting records. Allocated core-walltime is distinct from measured CPU use.

| Slurm job | CPUs | Memory | Limit | Actual CPU (s) | Elapsed (s) | Peak RSS | seff CPU efficiency | Outcome |
| --- | ---: | --- | --- | ---: | ---: | --- | ---: | --- |
| 55932113 | 2 | 12G | 01:00:00 | 1 | 10 | 209.25 MB | 5.0% | environment import failed before science |
| 55932211 | 2 | 12G | 01:00:00 | 19 | 304 | 473.55 MB | 3.12% | test collection error before science |
| 55932614 | 1 | 6G | 01:00:00 | 72 | 600 | 2.14 GB | 12.0% | initial network topology failed |
| 55932951 | 1 | 3G | 01:00:00 | 181 | 866 | 2.14 GB | 20.9% | geometry passed sharp network time tests failed |
| 55933300 | 1 | 3G | 01:00:00 | 211 | 1742 | 2.14 GB | 12.11% | mathematical checks passed pf campaign unresolved |

## Exact immutable identities

Audited base: `9f66c8d7a5a266687284d8da35aefbc6062808f7`. Historical production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`.

Initial NPZ SHA-256: `106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`.

Initial metadata SHA-256: `ff66cb5f938b5b70fe9e611674956f1e45f30f0979714507a853750f61e1a9bb`.

Both inputs were verified locally and by the compute wrappers. The matched initial state was not anisotropically equilibrated. Runner IDs contain `nogit` because the dedicated orchestration state is outside Git; the separate source commit and archive hash below identify the scientific source unambiguously.

### Job 55932113

Run ID: `20260911T002118Z-nogit-d279e0`.

Scientific Git commit: `db89a772ca01076abf177b0ce0c31cfedb8c2834`.

Source archive SHA-256: `a99aebe9f1d59d8a1b11fe194aa51645bf729d89cb7ffcf9305399f8df4e25e3`.

Result archive SHA-256: `4e1f446f06dbd0248485a0c0721f23b1f5a04ef1f59c83b2ee6b9b5465b85b73`.

Retrieval: `fetched_twice`; audit: `archive_checksum_verified_application_failed`.

### Job 55932211

Run ID: `20260911T002932Z-nogit-2371f7`.

Scientific Git commit: `f4aa53f2725c7f9bb3523e1165fe80203d5f487d`.

Source archive SHA-256: `6359b2660309bb3244706f9b16843b36bfc4a5b166346715401a865de076eb0c`.

Result archive SHA-256: `a0d61c5cf56a4406bd1ddbc6a49a1fac722924523cfafe37aced80db262bd01e`.

Retrieval: `fetched_twice`; audit: `archive_checksum_verified_application_failed`.

### Job 55932614

Run ID: `20260911T004050Z-nogit-3e2ea8`.

Scientific Git commit: `f4aa53f2725c7f9bb3523e1165fe80203d5f487d`.

Source archive SHA-256: `6359b2660309bb3244706f9b16843b36bfc4a5b166346715401a865de076eb0c`.

Result archive SHA-256: `901bd3d767793a577db1c31769162ee8c9e9764960c0f9e3471d57be3fc1d677`.

Retrieval: `fetched_twice`; audit: `archive_checksum_verified_application_failed`.

### Job 55932951

Run ID: `20260911T010258Z-nogit-ae4ebc`.

Scientific Git commit: `81c8785e516ae1ab480ca8a2916479f519d7192e`.

Source archive SHA-256: `96424071a82b2675bc136fa4efbfd4f4dacbf3c024cd3908845c090ebad63501`.

Result archive SHA-256: `86e69f6d877bef320b26e82dcf5d7e5f3ca6212e475fc93e49e2a22789ac8168`.

Retrieval: `fetched_twice`; audit: `archive_checksum_verified_application_failed`.

### Job 55933300

Run ID: `20260911T013122Z-nogit-e7f36e`.

Scientific Git commit: `54189b6fc7767535018afc5112bbffb3f47b8831`.

Source archive SHA-256: `7eb06d15256487ef35df2a3b6cb85babbc60c25a2f9ddf1c50718b62b5a27b54`.

Result archive SHA-256: `de36bf39e1d2fc6862f0a3166cb9187ab3c9c322938bb80b703cc64e103b0923`.

Retrieval: `fetched_twice`; audit: `finalized_archive_checksum_verified`.

## Environment and retained paths

Initial environment verification job **55932113** failed on missing pandas. Environment setup job **55932211** created the separate overlay and verified package versions, then failed test collection because the wrapper omitted the source root from PYTHONPATH. Later jobs reuse that overlay unchanged:

`/pub/sdillon1/sw/pfgg-anisotropy/f4aa53f-714dd16b2f3624b8`

Python 3.13.5; NumPy 2.1.3; SciPy 1.15.3; Numba 0.61.0; pandas 2.2.3; PyArrow 19.0.0; PyYAML 6.0.2; pytest 8.3.5; Matplotlib 3.10.1. The original Conda environment and QIU environments were not modified. Each retrieved output includes the environment fingerprint, package inventory, CPU details, NumPy build information and runtime libraries.

All local retrieved jobs are under:

`/Users/sdillon/HPC3/anisotropic-phase1-state/hpc3-results/pfgg-anisotropic-cahn-hoffman-v1/<run-id>/`

All remote jobs remain under:

`/pub/sdillon1/codex-runs/pfgg-anisotropic-cahn-hoffman-v1/<run-id>/`

Each terminal root retains Slurm stdout/stderr, `sacct-seff.txt`, finalization JSON, result archive and its checksum. Small result archives are decoded into `decoded/output/` locally after hash verification. No source archive was extracted locally. No remote result was deleted.

## Failure and ownership audit

Jobs 55932113 and 55932211 failed before scientific execution. Job 55932614 passed 190 tests but failed periodic topology; job 55932951 passed 193 tests and fixed topology but failed sharp-network time tests. These failures remain visible. Their archives were independently checksum verified despite the runner correctly retaining `partial_unverified` for failed applications.

Remote directory and status-marker operations intermittently timed out. Before any submission retry, local records and Slurm were reconciled; retries reused the same immutable bundle. Job 55933300 had two directory timeouts before upload/sbatch and was subsequently submitted exactly once. The first status timeout caused no resubmission. A failed exploratory `lfs quota` command was replaced by the existing runner quota tool, which returned the BeeGFS quota successfully.

Superseded prepared run `20260911T005754Z-nogit-669634` has no Slurm ID and incurred no allocation. The locked, atomic ownership ledger enforces one prepared run to one Slurm ID and unique trajectory ownership. No anisotropic production trajectory was launched, so no duplicate or restarted production chain exists.

Only approved seconds-scale static unit/derivative checks, syntax/collection checks, file hashing and lightweight report plotting ran locally. All full tests, initial-network scans and evolving network calculations ran in compute allocations. No PF evolution or large scientific analysis ran locally or on login nodes. Login access was limited to scheduler/account/filesystem control, and transfers used `uci-hpc3-transfer`.

QIU jobs 55930486 (216 CPU-hour maximum on SDILLON1) and 55932457 (720 CPU-hour maximum on its separate SDILLON1_LAB authorization) were read-only observations. This session issued no QIU cancellation, reprioritization, task assignment, artifact write or source edit. Its worktree and status were inspected read-only for shared-core reconciliation. All 75 pre-existing source files in this anisotropic branch are byte-identical to the audited base.

## Final reconciliation and diagnostic steps

Job 55933300 passed all 193 tests in 1219.61 s and completed the tested
sharp-network time gates. Its plotting regression temporarily waited in
BeeGFS I/O, then recovered before any cancellation. Diagnostic Slurm steps
55933300.0, 55933300.1 and 55933300.2 ran inside the existing allocation, each
with a one-minute limit and zero elapsed seconds at displayed accounting
precision. They read progress/process state and staged a diagnostic archive;
they were not additional scientific trajectories or allocations.

The separately staged diagnostic archive verifies as
`f2e66cdeab1d23a600f4ca964fee650e4ee93a3a5684d77729d68f753c5b3a65`.
Despite its `diagnostic-before-cancel.tar.gz` filename, no cancellation occurred.
The locally staged geometry-v7 fallback was never assigned a saved runner plan
or Slurm ID. It is marked NOT_SUBMITTED and was not charged.

The final job archive verifies as
`de36bf39e1d2fc6862f0a3166cb9187ab3c9c322938bb80b703cc64e103b0923`.
Application and finalization exit codes are zero and its complete marker is true.
An idempotent second fetch succeeded, and reconciliation reports no unretrieved
jobs. All older failed archives remain retained and explicitly failed.
The last live query showed only the two independently owned QIU jobs and an
SDILLON1 balance of 545 SU. No anisotropy job remains active.
