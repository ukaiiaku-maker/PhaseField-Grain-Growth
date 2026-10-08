# Matched anisotropic Phase-1 event-production status

Status timestamp: 2026-10-08 09:53 PDT

This campaign is the corrected execution of the original matched-isotropic
Phase-1 experiment. It is separate from the native-Qiu campaign. The Qiu
outputs do not substitute for these trajectories.

## Fixed scientific identity

- Branch: `codex/aniso-phase1-event-production-v1`
- Source commit: `8310877207fbee26e5eea6ae14617d10d4d0a10a`
- Configuration: `configs/production/anisotropic_phase1_event_900K.yaml`
- Initial state: exact 384 x 384, 800-grain, seed-5101 matched state
- Initial-state SHA-256:
  `106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`
- Temperature: 900 K
- Anisotropy: A2 strong, energy and mobility enabled
- Support constraint: compact active set, KKT tolerance `1e-10`
- Stop: physical time 4000 or grain count at most 100
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

All eight exact-state HPC3 preflights completed under runner ID
`20261008T162623Z-nogit-a56629`, Slurm array `57937613_[0-7]`. Results were
retrieved twice and checksum-verified. All eight identity reports passed. Peak
RSS was 7.36-7.80 GiB, so each production task requests 32 GiB. The focused
compact-support, anisotropy, restart, and accepted-clock test set has 22 passing
tests.

The `GTSC_GBTJ_Ks030_LONG` preflight is much slower than the other conditions:
about 148 seconds per step versus 8.7-14.0 seconds per step. It remains a full
production condition and is not being replaced by a short run.

## Production submissions

| Conditions | Runner ID | Slurm identity | State at timestamp |
|---|---|---|---|
| `B0`, `G`, `T`, `GT`, `GTC_GB`, `GSC_GBTJ_Ks025`, `GTSC_GBTJ_Ks025` | `20261008T165132Z-nogit-9fa9b5` | `57937717_[0-6]` | all seven running |
| `GTSC_GBTJ_Ks030_LONG` | `20261008T165136Z-nogit-bc3997` | `57937722_0` | running |

Each task requests one CPU, 16 GiB memory, 120 GiB node-local temporary space,
and the standard partition's 14-day limit. A planned interrupt at 13 days 21
hours preserves an exact checkpoint for continuation without replay. There is
no self-imposed two-worker ceiling; the seven-task array may run seven
conditions concurrently. HPC3 bills each 16 GiB task as three CPU slots. The
request retains more than a two-times margin over the measured 7.8 GiB peak.

The initial submissions `57937688` and `57937710` requested 32 GiB on the small
`SDILLON1` allocation and were held for `AssocGrpBillingMinutes`; both remained
scientifically unstarted and were cancelled. Their replacements use the
user-accessible `SDILLON1_LAB` allocation, which reported 90,565 available SU.
Four unrelated jobs (`57937556_[0-3]`) on that account were not modified. The
unused prepared plan `20261008T165121Z-nogit-bcedab` retained the old resources
after a local path error and must not be submitted.

All eight replacements entered `RUNNING` concurrently. At about one minute,
their batch steps had accumulated approximately 52 CPU seconds each and used
4.0-4.6 GiB average RSS, with 7.32-7.80 GiB maximum RSS. This confirms that the
scientific processes entered their workloads and that the 16 GiB requests have
adequate initial headroom. No full production trajectory is yet complete;
preflight completion must not be reported as production completion.
