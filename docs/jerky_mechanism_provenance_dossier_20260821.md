# Jerky grain-growth provenance and mechanism dossier

## Executive status

The corrected 900 K evidence supports C as the minimum existing mechanism for
reproducing GTSC kinetics and release-coupled intermittency.  It does not yet
support a distinct `C_GBTJ` defect-sink model: all free-volume demand and climb
accommodation are local to GB domains, even in TC/GTC/GTSC.  TJ compatibility
is a motion gate and a work sampler, not a point-defect sink.

Accordingly, direct final-checkpoint occupancy and reconstructed time occupancy
are reported separately below.  The 36-run temperature/seed campaign remains
on hold.

## Source SHAs and evidence generations

| Evidence | Recorded source SHA | Status and use |
|---|---|---|
| Legacy baseline/climb/long/video | `b2278ec695ba91a3ed6afa042b615585110975e1-dirty` | Pre-gate-only, pre-arclength legacy evidence; continuation context only |
| Migration-closure comparison | `10f4ebdab65967315956cd679729856196765d3a-dirty` | Gate-only versus hybrid closure audit; before final climb/TJ integrity fixes |
| Preproduction convergence | `238ea9d7c518958e95e86b8e05290a1c63ee9fbf-dirty` | Gate-only resolution/time-step/size checks |
| Physical TJ-width movie | `e58d03838b246301582d8ae8595216eafe55dd69-dirty` | Physical-width TJ gate validation |
| Corrected integrity smoke and factorial | `027b91cc417b0be17bf98f168a770613ef04dff8-dirty` | Corrected C gate separation, local arclength kinematics, TJ work |
| Corrected discriminating screen | `2265d3212a337e1857da10edba14257288cb4615-dirty` | G/T barriers, S packets, and C stage-barrier sweeps |
| Expanded causal analysis implementation | `08e47a0d74a98198740fb2291146d9dea9b06139` | Event-null, survival, work, factorial, topology, occupancy reconstruction |
| Expanded scientific report | `913982da30b38bb2babeef528e4f29db494aa729` | Prior integrated conclusions |

`-dirty` is part of each simulation manifest's recorded provenance.  It means
the base commit is known but the run cannot be represented as a clean checkout
of that SHA alone.  The complete run manifests and embedded resolved configs
must accompany every simulation claim.

## Full test result

Command:

```text
PYTHONPATH=src /opt/anaconda3/bin/python -m pytest -vv \
  --junitxml=results/validation/jerky_mechanism_full_tests_20260821.xml
```

Result: **124 passed, 0 failed, 0 errors, 0 skipped** in 193.39 s.

- Geometry: 8
- Integration: 43
- PF kernels/kinetics: 13
- Stochastic clocks: 8
- Unit/analysis/mechanics: 52
- JUnit record:
  `results/validation/jerky_mechanism_full_tests_20260821.xml`
- JUnit SHA-256:
  `c574ea5f40981ff9556a9d06021f0c43e29c7b07083fb14cdde699eb89b8ff44`
- Runtime: Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3, pandas 2.2.3,
  PyArrow 19.0.0, pytest 8.3.4.

The record includes the local free-volume conservation/scaling test, exact
serial-stage mean test, gate-only single climb-release test, GB/TJ separation,
physical TJ-width gate test, activation-work tests, restarts, compiled PF
equivalence, and the new causal-null/factorial/survival tests.

## Climb quota conservation derivation

For GB arclength domain (d), the corrected local swept measure over one
physics update is

\[
\Delta A_d^{swept}=|v_{n,d}|\Delta t\max(L_d,\Delta x).
\]

With excess volume per GB measure \(\delta_V^{GB}\) and point-defect formation
volume \(\Omega_{pd}\), the required defect quota is

\[
\Delta N_{req,d}=\frac{\delta_V^{GB}}{\Omega_{pd}}
\Delta A_d^{swept}.
\]

The state variables are cumulative required quota (R_d), accommodated quota
(A_d), and deficit (q_d=R_d-A_d).  Migration updates

\[
R_d^{k+1}=R_d^k+\Delta N_{req,d}.
\]

At serial climb completion, the accepted quota is

\[
a_d=\min(Q_{release},\max(q_d,0)),\qquad
A_d^{k+1}=A_d^k+a_d.
\]

Therefore the local invariant is exact by construction:

\[
R_d-A_d-q_d=0.
\]

The stored free-volume energy and chemical potential are

\[
E_q=\tfrac12K_q q_d^2,\qquad \Delta\mu_d=K_q q_d.
\]

The serial residence time is

\[
\langle t_C\rangle=r_{nuc}^{-1}+r_{ex}^{-1}+r_{tr}^{-1},
\]

not the inverse of the summed rates.  The gate-only correction accommodates
the quota once in `_advance_climb`; the release-summary event passes zero
additional quota into the legacy base writer.  This is exercised by
`test_gate_only_climb_completion_does_not_accommodate_quota_twice`.

This is **local quota conservation**, not global material conservation.  The
model uses `abs` of swept area, has no vacancy/interstitial sign or global
reservoir, and retires local domain state when the geometric entity disappears.
It therefore cannot yet demonstrate global point-defect conservation across
topology changes or exchange between GB and TJ sinks.

## Density and excess-volume parameters

The corrected smoke, factorial, and discriminating runs use:

| Parameter | Value | Meaning |
|---|---:|---|
| `event_domain_length` | 12.0 | Target physical connected-GB arclength per stochastic domain |
| `encounter_density` | 0.20 | Geometric encounter hazard per unit local swept measure; not a calibrated 3-D number density |
| `excess_volume_per_area` | 0.01 | \(\delta_V^{GB}\) |
| `point_defect_formation_volume` | 0.02 | \(\Omega_{pd}\) |
| ratio \(\delta_V^{GB}/\Omega_{pd}\) | 0.5 | Required quota per unit local swept measure |
| `climb_trigger_quota` | 0.25 | Blocking threshold; equivalent to 0.5 swept-measure unit from zero deficit |
| `climb_release_quota` | 1.0 | Maximum quota accommodated per completed serial cycle |
| `free_volume_stiffness` | 0.05 | \(K_q\), giving \(E_q=0.025q^2\) and \(\Delta\mu=0.05q\) |
| nucleation barrier/prefactor | 0.45 eV / \(10^5\) | Reference serial stage |
| exchange barrier/prefactor | 0.55 eV / \(10^5\) | Reference serial stage |
| transport barrier/prefactor | 0.65 eV / \(10^5\) | Reference serial stage |
| `tj_correlation_length` | 2.0 | Physical TJ gate half-width |

No `sink_density`, `GB_sink_fraction`, `TJ_sink_fraction`, or equivalent
partition parameter exists in the current source or resolved configs.

## Sink partitioning

Current partitioning is structural rather than configurable:

- Free-volume demand: 100% assigned to the migrating GB arclength domain.
- Serial nucleation/exchange/transport cycle: 100% owned by that GB domain.
- Accommodated quota: 100% returned to that GB domain's `FreeVolumeState`.
- TJ point-defect accommodation: 0% because no TJ free-volume reservoir or
  sink flux exists.
- TJ activation work can average chemical potential from adjacent GB domains,
  but this changes a TJ compatibility hazard; it does not consume quota.

Consequently, TC/GTC/GTSC test climb plus TJ compatibility, not a partitioned
GB+TJ defect-sink model.  A genuine `C_GBTJ` result requires a conserved global
or shared reservoir and explicit GB/TJ allocation counters.

## Legacy continuation roots

The following are the legacy generation retained for comparison:

- Baseline reconciliation:
  `results/campaigns/20260818T041229Z-49171551f9`
- Climb backpressure sweep:
  `results/campaigns/20260818T121303Z-cf04b62e5c`
- Long-horizon selected mechanisms:
  `results/campaigns/20260818T220844Z-bed34ea2da`
- Representative video archive:
  `results/video_runs/20260818T190827Z-409bd74b2d`

Their manifests report no `resumed_runs`; they are standalone campaigns built
from shared equilibrated initial-condition files, not checkpoint continuations
in the strict storage sense.  They are called continuation roots only in the
project-history sense and must not be pooled with corrected gate-only results.

Intermediate migration-closure roots are:

- Matched closure campaign:
  `results/migration_closure/20260819T153213Z-5f14d1e295`
- Matched closure movies:
  `results/migration_closure_video/20260819T213538Z-76352e8475`

Preproduction nonreplicate checks are:

- Reference: `results/preproduction_convergence/20260820T032540Z-8d9f46eb85`
- Half dt: `results/preproduction_convergence/20260820T062834Z-092c4a8f29`
- Fine grid: `results/preproduction_convergence/20260820T074800Z-7180f14944`
- Larger size: `results/preproduction_convergence/20260820T083720Z-5d5bacaef7`

## Corrected C_GB roots

`C_GB` is a descriptive alias here for the existing GB-local C implementation,
not a literal regime name.

- Corrected frame-preserving smoke campaign:
  `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0`
- Smoke C run:
  `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0/C-T900-s5101`
- Corrected 17-case factorial:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c`
- Factorial C run:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c/C-T900-s5101`
- Corrected stage-barrier screen:
  `results/jerky_discriminating_900K/20260821T080200Z-124f03cd0d`
- Stage variants in that root: `C_FAST-T900-s5101`, `C_REF-T900-s5101`,
  and `C_SLOW-T900-s5101`.

## Corrected C_GBTJ roots

There is no corrected root for a strict partitioned `C_GBTJ` sink model.
Operational climb-plus-TJ-gate proxies in the corrected factorial are:

- Pure TC proxy:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c/TC-T900-s5101`
- GTC proxy:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c/GTC-T900-s5101`
- Full GTSC proxy:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c/GTSC-T900-s5101`
- Frame-preserving GTC and GTSC proxies:
  `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0/GTC-T900-s5101`
  and
  `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0/GTSC-T900-s5101`.

Calling these `C_GBTJ` without the qualifier “TJ compatibility proxy” would
incorrectly imply TJ point-defect sink partitioning.

## Movies and indexes

The complete path-level inventory is in
`docs/jerky_mechanism_movie_index_20260821.md`.  Its authoritative indexes are:

- `results/video_runs/20260818T190827Z-409bd74b2d/video_manifest.json`
- `results/migration_closure_video/20260819T213538Z-76352e8475/video_manifest.json`
- `results/tj_gate_radius_minimal/20260820T142900Z-aa3c2c9e8f/video_manifest.json`
- `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0/video_manifest.json`

The corrected factorial and discriminating roots contain no frame archives and
therefore no movies.

## Corrected stage-resolved climb kinetics at 900 K

Fractions below are event-conditioned sums of expected residence `1/r`; they
are diagnostic resistance shares, not causal time occupancy.

| Case | Stage | Events | Median rate | Mean expected wait | Resistance share |
|---|---|---:|---:|---:|---:|
| C factorial | nucleation | 17,247 | 302.072 | 0.00331 | 0.124% |
| C factorial | exchange | 17,239 | 13.651 | 0.07164 | 2.685% |
| C factorial | transport | 16,913 | 0.4677 | 2.6432 | 97.191% |
| GTSC | nucleation | 15,820 | 302.072 | 0.00331 | 0.124% |
| GTSC | exchange | 15,811 | 13.637 | 0.07179 | 2.685% |
| GTSC | transport | 15,521 | 0.4677 | 2.6469 | 97.191% |
| C_FAST | nucleation/exchange/transport | 17,859 / 17,855 / 17,624 | 1096.696 / 49.960 / 1.698 | 0.00091 / 0.01961 / 0.7312 | 0.123% / 2.642% / 97.235% |
| C_REF | nucleation/exchange/transport | 8,082 / 8,077 / 7,740 | 302.072 / 13.666 / 0.4677 | 0.00331 / 0.07149 / 2.5587 | 0.131% / 2.830% / 97.039% |
| C_SLOW | nucleation/exchange/transport | 3,254 / 3,235 / 2,831 | 83.202 / 3.758 / 0.1288 | 0.01202 / 0.2603 / 8.6876 | 0.154% / 3.305% / 96.541% |

The near-identical C and GTSC shares demonstrate that TJ compatibility does not
partition the climb sink or alter the reference serial bottleneck.  Transport
controls about 97% of event-conditioned serial resistance at all three barrier
levels.

## Event-triggered responses

Saved boundary histories were aligned from the first observed pinned frame to
the exact release row.  Median effective-barrier changes are:

| Case/event | Histories | Median duration | Median Delta effective barrier | p10 to p90 |
|---|---:|---:|---:|---:|
| G compatibility | 6,297 | 1.189 | -0.00069 eV | -0.00968 to +0.00508 eV |
| GT compatibility | 6,045 | 1.165 | -0.00064 eV | -0.00926 to +0.00458 eV |
| C quota completion | 16,562 | 1.554 | -0.00049 eV | -0.00874 to +0.00566 eV |
| GTSC compatibility | 1,751 | 1.296 | -0.00036 eV | -0.00716 to +0.00361 eV |
| GTSC quota completion | 15,162 | 1.532 | -0.00097 eV | -0.00915 to +0.00450 eV |

At 900 K, `kBT=0.0776 eV`; median changes are under 1.3% of kBT.  The current
output cadence does not show systematic stress loading before most releases.
The full binned response is
`results/production_summaries/jerky_mechanism_expanded/pre_release_barrier_trajectories.csv`.

## Causal-null results

The large-burst definition is the strictly positive top 5% of grain-radius-rate
increments.  The null preserves each grain's number of event-associated
intervals and randomizes their locations over 200 deterministic shuffles.

| Case | Actual release risk ratio | Shuffled ratio | Excess | Top-5% preceded-fraction excess | Release-leading cross-correlation excess |
|---|---:|---:|---:|---:|---:|
| G | 1.126 | 0.991 | 0.135 | 0.0155 | 0.294 |
| T | 1.143 | 1.087 | 0.055 | 0.0019 | 0.207 |
| C | 1.561 | 0.997 | 0.564 | 0.1125 | 0.169 |
| GT | 1.206 | 0.989 | 0.218 | 0.0280 | 0.244 |
| GTSC | 1.546 | 1.004 | 0.542 | 0.1121 | 0.128 |
| B0/S/QIU | no releases | -- | 0 | 0 | 0 |

C and GTSC are essentially identical by the top-motion causal metric, while
high CV/stationarity in S and QIU has no event cause.  The detailed table is
`results/production_summaries/jerky_mechanism_expanded/event_burst_causality.csv`.

## Kinetic crossover

The legacy fully physical temperature analysis placed a change in dominant
microscopic resistance near 900 K, but that inference belongs to the
`b2278ec...-dirty` generation and predates gate-only/arclength/TJ-work fixes.  It
must not be treated as a corrected crossover temperature.

The corrected data contain only 900 K and therefore cannot locate a temperature
crossover.  They do establish a parameter crossover:

- G_LOW/REF/HIGH: linear K `0.04316 / 0.03395 / 0.01938`; median GB pin time
  `0.48 / 1.28 / 2.88`.
- T_LOW/REF/HIGH: linear K `0.04744 / 0.04007 / 0.02268`; stationary fraction
  `0.746 / 0.772 / 0.852`.
- GT_LOW/REF/HIGH: linear K `0.04162 / 0.02928 / 0.01302`; the high-barrier
  finite window begins to prefer the parabolic comparator.
- C_FAST/REF/SLOW: linear K `0.01530 / 0.00845 / 0.00818`; median pin time
  `0.48 / 1.44 / 4.64`.

In the full factorial, GT has series-fit `KT=0.03355` with no finite fitted KG,
whereas C and GTSC have finite series-fit `KG=0.10026` and `0.10071` with no
finite KT.  This finite-window switch marks climb-controlled kinetics at 900 K;
it is not yet a temperature-dependent crossover.

## Direct and reconstructed occupancy

### Direct final-checkpoint snapshot

The final checkpoint persists exact `compatibility_pending`, free-volume totals,
climb stage, GB `blocked`, and TJ `blocked` state.

| Case | Final GB domains | Free | G only | C only | G+C | Blocked TJs / all TJs |
|---|---:|---:|---:|---:|---:|---:|
| G | 388 | 340 | 48 | 0 | 0 | 0 / 0 |
| T | 356 | 356 | 0 | 0 | 0 | 8 / 140 |
| C | 616 | 332 | 0 | 284 | 0 | 0 / 0 |
| GT | 393 | 341 | 52 | 0 | 0 | 13 / 164 |
| GC | 604 | 322 | 11 | 271 | 0 | 0 / 0 |
| TC | 615 | 315 | 0 | 300 | 0 | 4 / 320 |
| GTC | 606 | 313 | 15 | 276 | 2 | 5 / 319 |
| GTSC | 621 | 354 | 11 | 254 | 2 | 8 / 328 |

Among the 267 directly blocked GTSC GB domains at the final checkpoint, 95.13%
are C-only, 4.12% G-only, and 0.75% simultaneously G+C.  TJ blockage is a
separate 2.44% of final TJ domains and cannot be added to the GB-domain
denominator.

### Reconstructed blocked-domain time

The event-classified time reconstruction gives GTSC 85.87% C-limited, 5.21%
G-limited, 1.97% simultaneous, and 6.95% unresolved.  This is not direct state
occupancy; it uses generic blocked tracks plus observed completion events.  The
agreement in dominance with the direct final snapshot is reassuring, but exact
time-integrated G/T/C occupancy requires reason-specific fields at every output.

## Revised factorial comparisons

The design contains 16 G/T/S/C subsets including B0 plus QIU, for 17 completed
cases total.  Descriptive high-minus-low contrasts are:

- Jerkiness CV: C `+1.252`, S `+0.642`, SxC `-0.597`.
- Stationary fraction: C `+0.0802`, S `+0.0676`, SxC `-0.0651`.
- Linear K: C `-0.02137`, S `-0.01868`, SxC `+0.01846`.
- Event Fano: C `+10.84`, G `+6.87`, GxC `-7.35`.
- Top-5% causal release excess: C `+0.0934`, TxC `-0.0112`, G `+0.0094`.

At common topology N=190 to 160, C versus GTSC gives:

- `dR/dt`: `0.005486` versus `0.005241`.
- Stationary fraction: `0.9634` versus `0.9630`.
- Pinned fraction: `0.4946` versus `0.4756`.
- Release Fano: `1.148` versus `1.118`.

Thus C approximately equals GTSC is not an equal-duration/topology artifact.
The large negative SxC interactions show that the serial C gate masks most
incremental shear response.

## Revised minimum model

The standardized distance uses linear K, stationary fraction, top-5% motion
concentration, event Fano, median pin duration, and top-5% causal release
coupling.

- C: discrepancy `0.133`, Pareto-optimal one-mechanism model.
- TC proxy: `0.130`; adding a TJ gate changes little and does not add a TJ sink.
- SC: `0.095`, Pareto-optimal two-mechanism model.
- TSC: `0.067`, Pareto-optimal three-mechanism model.
- GTSC: `0` by definition.

C remains the revised minimum existing model.  The result should be described
as “GB-local serial climb reproduces the selected GTSC observables at 900 K,”
not as proof that a physical C_GBTJ sink network is unnecessary.

## Recommendation for the next temperature/seed campaign

Do not launch the original 36 runs yet.

First add:

1. explicit signed point-defect inventory or a documented unsigned surrogate;
2. configurable GB/TJ sink allocation with fractions summing to one;
3. persisted cumulative demand, GB accommodation, TJ accommodation, retired
   inventory, and conservation residual;
4. reason-specific G/C/T pending flags and TJ tracks at ordinary output cadence;
5. stable event arclength/spatial coordinates for branching analysis.

Then run two 900 K/seed-5101 diagnostic replays, `C_GB` and true `C_GBTJ`, and
require global conservation plus direct time occupancy before a temperature
ensemble.

If they pass, the smallest informative temperature screen is 18 runs:

- regimes: `C_GB`, `C_GBTJ`, and GTSC;
- temperatures: the lowest and highest planned temperatures;
- seeds: three paired seeds.

Use the existing 900 K B0/C/GTSC evidence as the center anchor.  Add B0 at the
two extremes only if intrinsic mobility is made temperature dependent; with the
current constant intrinsic mobility it does not need six repeated temperature
controls.  Fill the two interior temperatures only after the extreme screen
shows whether C_GB or C_GBTJ continues to track GTSC.
