# Anisotropic campaign handoff

This is an incomplete campaign with implemented and tested geometry foundations.
Classification: `ANISOTROPIC_IMPLEMENTATION_UNRESOLVED`. Do not release production.

## Locations and ownership

- Worktree: `/private/tmp/pfgg-anisotropic-cahn-hoffman-v1`.
- Branch: `codex/anisotropic-cahn-hoffman-phase1-v1` (pushed).
- Local state: `/Users/sdillon/HPC3/anisotropic-phase1-state`.
- Ledger: `ownership.json` there, with lock/fsync/atomic replacement transitions.
- Existing runner: `/Users/sdillon/HPC3/hpc3-runner/.venv/bin/hpc3`.
- Remote root: `/pub/sdillon1/codex-runs/pfgg-anisotropic-cahn-hoffman-v1`.
- Historical comparison root: `/Users/sdillon/PF-graingrowth/results/long_time_kinetics_900K_20260824/20260825T012949Z-6d83ee5c82` (read-only).

## Immutable jobs

Job 55932113, runner `20260911T002118Z-nogit-d279e0`, used source
`db89a77` and source archive SHA-256
`a99aebe9f1d59d8a1b11fe194aa51645bf729d89cb7ffcf9305399f8df4e25e3`.
It failed before science due to missing pandas. Failed diagnostics were retrieved
twice, the archive hash verified, and remote artifacts were retained.

Job 55932211, runner `20260911T002932Z-nogit-2371f7`, uses source
`f4aa53f2725c7f9bb3523e1165fe80203d5f487d` and archive SHA-256
`6359b2660309bb3244706f9b16843b36bfc4a5b166346715401a865de076eb0c`.
Its immutable bundle is in `geometry-v2/`. The separate environment overlay
succeeded, but test collection failed because the wrapper omitted the source
root from PYTHONPATH. Its failed archive was verified and fetched twice.
The numerical environment is now reused read-only at
`/pub/sdillon1/sw/pfgg-anisotropy/f4aa53f-714dd16b2f3624b8`.

Job 55932614, runner `20260911T004050Z-nogit-3e2ea8`, uses the same scientific
source and corrected wrapper. All 190 tests passed; exact initial reconstruction
then failed with an unclosed periodic boundary. Its archive was verified and
an idempotent second fetch completed after one monitoring timeout.

Job 55932951, runner `20260911T010258Z-nogit-ae4ebc`, uses scientific source
`81c8785e516ae1ab480ca8a2916479f519d7192e` and archive SHA-256
`96424071a82b2675bc136fa4efbfd4f4dacbf3c024cd3908845c090ebad63501`.
It tests minimum-image vertex matching and phase-value pullback, followed by
independent initial geometry and manufactured sharp-network evolution.
All 193 tests and A2 constitutive gates passed. Sharp time checks failed: A2
balance error above 1% and an unequal-physical-time TJ iteration budget.
Its verified archive is
`86e69f6d877bef320b26e82dcf5d7e5f3ca6212e475fc93e49e2a22789ac8168`.

Job **55933300**, runner `20260911T013122Z-nogit-e7f36e`, uses source
`54189b6fc7767535018afc5112bbffb3f47b8831`, archive SHA-256
`7eb06d15256487ef35df2a3b6cb85babbc60c25a2f9ddf1c50718b62b5a27b54`.
It refines loop time and matches TJ physical horizons without changing A2 or
acceptance thresholds. Latest state: SUBMITTED. Request: 1 CPU, 3 GiB, 1 h on
SDILLON1 / standard. Two pre-upload directory timeouts caused no Slurm submission;
the same immutable bundle was reconciled and then submitted exactly once.
The scripts refuse scientific execution outside a Slurm compute allocation.

The regression-v4 plan `20260911T005754Z-nogit-669634` was superseded before
submission. It has no Slurm ID; do not submit it.

For bounded status and retrieval, run from the dedicated state directory:

```sh
/Users/sdillon/HPC3/hpc3-runner/.venv/bin/hpc3 status 20260911T013122Z-nogit-e7f36e
/Users/sdillon/HPC3/hpc3-runner/.venv/bin/hpc3 fetch 20260911T013122Z-nogit-e7f36e
/Users/sdillon/HPC3/hpc3-runner/.venv/bin/hpc3 reconcile --fetch-completed
```

Do not submit another copy on an SSH/status timeout. Use the known job ID and
preserved runner records. Scientific completion requires terminal Slurm state,
application success, complete finalization and verified output hashes. The
manufactured flow is a polygon model; it does not qualify PF coupling.

## Remaining implementation

The historical PF kernel has scalar spatial mobility and a discrete capillary
rate that is not simply geometric pressure. Establish a variational native
kernel or a work-conjugate reference correction before connecting the new
geometry. Audit pair mobility, activation pressure, TJ endpoint terms, physical
time cadence, checkpoint state and stochastic clocks together. Preserve exact
A0 nesting and do not invent a new constitutive family or use A3 prematurely.

QIU shared-core changes were reviewed at clean commit `847eb06e`. A later
snapshot at `4ef9390` contains additional QIU-specific scripts/diagnostics and a
dirty spatial-convergence script, without further shared solver/kernel changes.
See `anisotropic_qiu_shared_core_review.md`. No changes were imported. Recheck
its durable status before production; never write to its worktree, jobs or
artifacts. Necessary selective integration requires repeated qualification.

The non-executable preregistration holds exact historical configs for all ten
trajectories. There is deliberately no runnable anisotropic production YAML:
the implementation and normalization have not qualified. Never run the old
isotropic solver with anisotropic case names as a substitute.
