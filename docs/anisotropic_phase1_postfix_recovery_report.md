# Phase-1 anisotropic postfix recovery

## Recovery identity and scope

This recovery uses branch `codex/aniso-phase1-postfix-recovery-v1`, created in
an isolated worktree directly from `1316cc89dbabfb41cb883b0d4a4c74738cc2bef6`.
The executable scientific equations are unchanged from that commit. Recovery
edits are limited to checkpoint/finalization code, focused tests, and reports.
The later Qiu/SI source is absent and outside this recovery.

## Job 55968528 archive audit

The canonical local record is
`/Users/sdillon/HPC3/anisotropic-phase1-state/hpc3-results/pfgg-anisotropic-cahn-hoffman-v1/20260912T212636Z-nogit-6d3d72`.
The partial archive was fetched twice with matching SHA-256
`f4b97f553ec5e6923bb205eb95700ecec63a7fb60b74dd6c3cd94edf749980f9`.
Slurm recorded `TIMEOUT`, 04:00:25 elapsed against 04:00:00 requested, two
allocated CPUs, 03:58:45 CPU time, 49.65% CPU efficiency, and 1.19 GiB peak
memory. The application wrapper had not returned; emergency finalization
created a partial archive. The final marker therefore records incomplete
finalization and cannot confer scientific completion.

The archived identities are:

- scientific source commit: `1316cc89dbabfb41cb883b0d4a4c74738cc2bef6`
- submitted source archive: `57197ed07f845a4823a745403dcb227466864ba32b04db21091667db59e42f6e`
- runner input bundle: `8d922df923de1bf2d57a305e16e309830837a73c7fbad9e495604e4adf3c6c25`
- captured transition state: `75f189987c95a94645b628548078ecb9d37c460ada921b81c8ad233658558a0c`
- environment manifest: `714dd16b2f3624b877352abbbe3e75a576883507e8b286dfd0691eb70d30438a`
- resolved job manifest: `54f5fe3d606a37f9c1a3252625a5212e41e306e9691b6113f91ba70c9b3bff65`
- job wrapper: `81271fe933acc28b72b0f79d6ce1c18d81a17990309a0732caadfad8d8ac179c`

The archive contains a complete `pytest.log` (`206 passed in 85.67s`) and a
passing manufactured `qualification.json`. Those results are reusable upstream
evidence because the scientific source is identical. It contains only progress
text for the polycrystal calculation and an empty `output/postfix/` directory.
There are no phase fields, active masks, orientations, physical times, accepted
step records, anisotropy/configuration identities, pair-support state,
timestep histories, checkpoint schemas, energy traces, final states, case
summaries, or restart comparisons.

| Requested record | Preserved progress | Classification | Recovery action |
|---|---:|---|---|
| A0 continuous | 64/64 | `LOG_ONLY` | rerun from initial state |
| A0 midpoint restart | replay not independently serialized | `LOG_ONLY` | rerun from new midpoint checkpoint |
| A2 energy-only continuous | 64/64 | `LOG_ONLY` | rerun from initial state |
| A2 midpoint restart | replay not independently serialized | `LOG_ONLY` | rerun from new midpoint checkpoint |
| half-dt A2 continuous | 96/128 | `LOG_ONLY` | rerun all 128 steps |
| half-dt A2 midpoint restart | no state | `LOG_ONLY` | rerun from new midpoint checkpoint |

No exact scientific state from job 55968528 is reusable or restartable.

## Finalization repair and preregistered decision

Each recovery case writes atomic, fsynced checkpoints at steps divisible by 16
and at stage boundaries. The checkpoint embeds the complete field, active mask,
orientations, mobility scale, exact step/time, solver configuration, energy and
support histories, timestep records, schema, source, input, configuration, and
environment identities. Pair support is a pure derivation of the field and
active mask and has no mutable solver cache. A separate exact midpoint and
continuous final state support the restart replay. `CASE_INCOMPLETE` remains
until the case summary, energy CSV, timestep CSV, and `CASE_COMPLETE` marker are
durable. `USR1` requests an exact checkpoint and never creates a completion
marker.

The timestep result is preregistered as converged when the matched-horizon
coarse/half-dt final-energy relative difference is at most 1% and each of
boundary density, grain-area CV, and area-weighted mean-radius differs by at
most 10%. Classification priority is energy failure, restart failure, timestep
nonconvergence, complete qualification, then operational incompleteness.

## Sizing and preflight

Job 55968528 did not timestamp its individual progress records, so it cannot
support separate defensible measurements for coarse A2, pre/post-saturation
half-dt A2, or restart replay. The only exact rate is the aggregate bound:
03:58:45 CPU time and 04:00:25 elapsed while completing the tests,
manufactured/actual-state work, A0 continuous/restart, A2 continuous/restart,
and 96 half-dt continuous steps. Assigning that full elapsed time to the 224
logged continuous accepted steps gives a conservative upper bound of 64.4
elapsed seconds per logged step; it is not represented as a case-specific
measurement. The recovery records monotonic per-stage timings so this omission
will not recur.

The recovery repeats 256 continuous steps plus 128 restart steps, omits the
already passed full test/manufactured work and actual-state audit, and adds
per-step diagnostics plus atomic I/O. Its conservative predicted runtime is
eight hours; twelve hours supplies the required 50% margin. One CPU and 3 GiB
limit exposure to 12 CPU-hours. On 2026-09-14 the `standard` partition was UP
with a 14-day maximum, the SDILLON1 balance was 477 SUs, the remote project root
existed, and no job with the recovery name/comment appeared in `squeue`. One
unrelated FFT job, 55950433, was running. No Qiu recovery or production job was
submitted.

## Completion result

Pending the single bounded recovery job.
