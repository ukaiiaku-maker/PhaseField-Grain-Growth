# Anisotropic campaign handoff

This is an incomplete campaign with implemented and tested geometry foundations.
Classification: `ANISOTROPIC_PF_TJ_ENERGY_GATE_FAILED`. Do not release reduced
or production trajectories with the current diffuse operator.

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
It refined loop time and matched TJ physical horizons without changing A2 or
acceptance thresholds. Final state: **COMPLETED, verified, fetched twice**.
All 193 tests and the tested sharp-network time gates passed. Request: 1 CPU,
3 GiB, 1 h on SDILLON1 / standard; elapsed 1742 s, CPU 211 s, peak RSS 2.14 GB,
CPU efficiency 12.11%. Archive SHA-256:
`de36bf39e1d2fc6862f0a3166cb9187ab3c9c322938bb80b703cc64e103b0923`.
Two pre-upload directory timeouts caused no Slurm submission; the same immutable
bundle was reconciled and then submitted exactly once. A later BeeGFS plotting
stall recovered before any cancellation. Diagnostic steps 55933300.0/.1/.2
read progress/process status and staged a partial archive inside its allocation.
A locally staged geometry-v7 fallback has no saved plan or Slurm ID and must
not be submitted as a duplicate.
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
snapshot at `4ef9390` contained additional QIU-specific scripts/diagnostics and a
dirty spatial-convergence script. The final read-only snapshot is clean at
`c68312799b4e41ba404cc0dd1ec3865bdcc9f9b8`, with no further source changes since
4ef9390.
See `anisotropic_qiu_shared_core_review.md`. No changes were imported. Recheck
its durable status before production; never write to its worktree, jobs or
artifacts. Necessary selective integration requires repeated qualification.

The non-executable preregistration holds exact historical configs for all ten
trajectories. There is deliberately no runnable anisotropic production YAML:
the implementation and normalization have not qualified. Never run the old
isotropic solver with anisotropic case names as a substitute.

## Final evidence and remaining scientific work

All five submitted anisotropy qualification jobs are terminal and retrieved.
Reconciliation reports no unretrieved jobs. The source tested on HPC3 is 54189b6;
later commits add reports and compact evidence. A2 passes the four provisional
initial-network constitutive gates. Refined sharp-loop dissipation errors are
0.517%/0.258%; TJ residuals are below 1e-7 and the refinement distance is 3.10e-10.
No production strength or continuum normalization is frozen.

Primary compact result root:
`results/long_time_kinetics_900K_anisotropic_20260910/20260911T013122Z-mathematical-qualification/`.
It contains the decision, parameters, provisional normalizations, constitutive
CSV, ten-row run/comparison matrices, two actual mathematical summary figures,
movie nonproduction index, checksum inventory and compact HPC evidence.
Earlier coarse failures remain in the 20260911T010258Z companion directory.

## Diffuse PF implementation outcome

The resumed branch implemented an opt-in pairwise anisotropic energy derivative,
pair mobility, diffuse capillary-pressure diagnostics and accepted-physical-time
propagation. HPC3 job 55949185 passed all 200 repository tests and most focused
diffuse checks, but failed A2 triple-junction energy descent. The focused
timestep job 55949331 showed maximum energy increases of 2.64325, 2.76135,
3.66392 and 3.30394 when accepted dt was reduced by 8x. This excludes the
stability cap as the explanation and fails the current active-set/projection
formulation structurally.

Compact evidence is in
`results/long_time_kinetics_900K_anisotropic_20260910/20260912T001247Z-pf-tj-gate-failed/`.
Full immutable output remains under runner ID
`20260912T001247Z-nogit-ff8282`, Slurm job 55949331, result archive SHA-256
`b70bb2188ac817112794eade58b98b0ef4c8b311bec7d86913e9f8957c9ab01e`.
All ten production rows remain `PREPARED_NOT_RELEASED`; reconstruction
normalization remains provisional. A future continuation must replace or
derive a continuous work-conjugate multiphase constraint treatment and repeat
the full diffuse/convergence gates before any reduced or production release.
