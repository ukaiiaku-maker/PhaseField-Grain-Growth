# Matched anisotropic Phase-1 event-production status

Status timestamp: 2026-10-26 16:31 PDT

This campaign is the corrected execution of the original matched-isotropic
Phase-1 experiment. It is separate from the native-Qiu campaign. The Qiu
outputs do not substitute for these trajectories.

## Fixed scientific identity

- Branch: `codex/aniso-phase1-event-production-v1`
- Original production source commit: `8310877207fbee26e5eea6ae14617d10d4d0a10a`
- GTC continuation source commit: `8754e382f38236c8514c34d3c356d2dbd8ccfa15`
- Coupled-run stability-repair commit: `2f5ad882b6102dacfe0904fe6b3e27461a4ed696`
- Configuration: `configs/production/anisotropic_phase1_event_900K.yaml`
- Initial state: exact 384 x 384, 800-grain, seed-5101 matched state
- Initial-state SHA-256:
  `106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`
- Temperature: 900 K
- Anisotropy: A2 strong, energy and mobility enabled
- Support constraint: compact active set, KKT tolerance `1e-10`
- Stop: grain count at most 100. No scientific physical-time cap is present in the continuation configuration
- Proposed timestep: 0.04; the physical clock advances by the timestep accepted
  by energy backtracking.

The mechanism matrix matches the original isotropic factorial. Controls omit
mechanisms by design; the three fully coupled regimes activate shear memory,
area-loss climb, GB sinks, TJ sinks, and the stochastic event engine.

| Regime | Stochastic events | Climb | Shear | GB sink | TJ sink |
|---|---:|---:|---:|---:|---:|
| `B0` | no | no | no | no | no |
| `G` | yes | no | no | no | no |
| `T` | yes | no | no | no | no |
| `GT` | yes | no | no | no | no |
| `GTC_GB` | yes | yes | no | yes | no |
| `GSC_GBTJ_Ks025` | yes | yes | yes | yes | yes |
| `GTSC_GBTJ_Ks025` | yes | yes | yes | yes | yes |
| `GTSC_GBTJ_Ks030_LONG` | yes | yes | yes | yes | yes |

Here, stochastic events are supplied by `multihit_persistent` with the
condition-specific GB/TJ compatibility modules. Identity verification also
requires checkpointed domain RNG state and the applicable event ledgers.

## Qualification

All eight original exact-state HPC3 preflights completed under runner ID
`20261008T162623Z-nogit-a56629`, Slurm array `57937613_[0-7]`. Results were
retrieved twice and checksum-verified. All eight identity reports passed. Peak
RSS was 7.36-7.80 GiB, supporting the 16 GiB production request. The focused
compact-support, anisotropy, restart, and accepted-clock test set has 22 passing
tests.

Those short preflights did not detect the forced-step stability defect described
below. After repair, a new exact 384 x 384 production-scale gate advanced all
three affected conditions for two steps. Both steps accepted `dt = 0.005`, all
three retained 800 grains, and the focused repair suite passed 18 tests.

The `GTSC_GBTJ_Ks030_LONG` preflight is much slower than the other conditions:
about 148 seconds per step versus 8.7-14.0 seconds per step. It remains a full
production condition and is not being replaced by a short run.

## Production outcomes and continuation

Four trajectories have valid terminal results. `GTC_GB` is advancing from its
exact timeout checkpoint. The three shear/climb trajectories previously marked
complete are scientifically invalid and have been replaced by repaired reruns.

| Regime | State | Grain count | Physical time | Disposition |
|---|---|---:|---:|---|
| `B0` | complete | 99 | 216.0025 | terminal |
| `G` | complete | 100 | 419.9675 | terminal |
| `T` | complete | 100 | 350.6650 | terminal |
| `GT` | complete | 100 | 535.30875 | terminal; original archive verified |
| `GTC_GB` | running | 108 at latest track output | 528.0050 | timestep stream has advanced to 531.505 |
| `GSC_GBTJ_Ks025` | invalid; repaired rerun running | 800 restart | 0 | Slurm array task 0 |
| `GTSC_GBTJ_Ks025` | invalid; repaired rerun running | 800 restart | 0 | Slurm array task 1 |
| `GTSC_GBTJ_Ks030_LONG` | invalid; repaired rerun running | 800 restart | 0 | Slurm array task 2 |

The original `GTC_GB` task `57937728` exhausted the standard partition's
14-day scheduler allocation. Its step-105921 checkpoint is intact, has SHA-256
`cc551e9ef2b3907976736c24baa6d8d927e80c174d177e3b13e50298bf16bc01`,
and contains the exact field state, active-phase mask, stochastic RNG states,
hazards, climb states, and event counters. The scientific continuation has
`termination_grains = 100`, `max_steps = 100000000`, and no
`maximum_physical_time` setting. HPC3 still requires a finite Slurm request;
the continuation uses the partition maximum of 14 days and will be chained
from its latest exact checkpoint if another scheduler segment is needed.

The timeout archive omitted early Parquet fragments from the detailed event
ledgers (`events` parts 0-3, `event_trace_events` parts 0-3, and
`event_traces` parts 0-65). The node-local scratch directory no longer exists.
This creates an explicit early-ledger gap for final mechanism statistics but
does not affect the checkpointed physical or stochastic state. Zero-row schema
placeholders allow the append-only ledgers to continue, and
`RECOVERED_LEDGER_GAPS.json` records the missing ranges.

Two continuation launches stopped before accepting a new timestep. Runner
`20261023T225217Z-nogit-17c561` exposed the non-contiguous Parquet parts;
runner `20261023T225631Z-nogit-d78f67` then exposed a legacy restart bug in
which active serial-climb stage history had not been checkpointed. Neither
attempt modified the step-105921 checkpoint. Commit `8754e38` persists climb
history for new checkpoints and reconstructs the diagnostic stage-history
lower bound for legacy checkpoints without changing the saved stochastic
clock. The relevant unit and resume/integration suite passes 22 tests.

The repaired exact continuation is runner
`20261023T230538Z-nogit-47b1c2`, Slurm job `58449253`. It was submitted with
16 GiB memory and is running on `hpc3-15-28`. Its
immutable lineage record binds the original source, continuation source,
checkpoint step, physical time, grain count, and checkpoint SHA-256. Once it
starts, acceptance requires an observed step greater than 105921 and physical
time greater than 520.005 before this restart is classified as advancing. The
process was found in Linux stopped state on 2026-10-26 after a prior diagnostic
attach had not resumed it. It had not changed the exact checkpoint, timestep
history, or 114-grain state. `SIGCONT` returned the same process to CPU-running
state on `hpc3-15-28`; stderr remains empty and memory remains below the 16 GiB
request. It is now confirmed advancing: the durable timestep stream reached
step 108221 and physical time 531.505, while the latest grain-track output at
time 528.005 contains 108 grains.

## Coupled-run instability repair

The prior `GSC_GBTJ_Ks025`, `GTSC_GBTJ_Ks025`, and
`GTSC_GBTJ_Ks030_LONG` results are invalid. Their compact-support steps accepted
the full `0.04` proposal whenever the mechanism layer supplied a nonzero
driving field, because forced steps bypassed energy backtracking. The controls
mostly accepted `0.00125-0.005`. The affected runs consequently collapsed from
800 grains to 220-247 grains by physical time 8 and to about 100 grains by time
13-14.4.

Commit `2f5ad88` replaces the bypass with deterministic backtracking of the
discrete interfacial energy minus the frozen external driving potential. It
also enforces a 0.1 phase-fraction trust region so a strong external drive
cannot project a cell directly across a large fraction of the simplex in one
step. The production configuration has no physical-time cap.

No valid restart-complete checkpoint exists between time zero and the unstable
final checkpoints. The intermediate NPZ files are visualization frames, not
restart states. All three replacements therefore use the exact matched
800-grain initial state at time zero, SHA-256
`106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`.
The invalid outputs and videos are quarantined under
`completed-results/invalidated-unstable` with an explicit invalidation record.

The repaired production array is runner `20261026T232929Z-nogit-a6542b`,
Slurm job `58521780_[0-2]`. All three 16 GiB tasks are running concurrently on
`hpc3-20-25`, have clean stderr, and verified the exact repair source. The only
scientific termination criterion is at most 100 grains. Because HPC3 enforces a
14-day scheduler maximum, the wrapper copies each new exact checkpoint and its
complete run history to shared storage for subsequent continuation segments.

A redundant replacement `GT` job (`58245823`) was cancelled after the original
`GT` result was proven complete at 100 grains.
