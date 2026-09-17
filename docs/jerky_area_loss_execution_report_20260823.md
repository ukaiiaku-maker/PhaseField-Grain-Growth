# Corrected area-loss climb execution report

## Scope and outcome

This campaign replaces the former swept-area climb source with a conserved
grain-boundary-area-loss source and keeps the former swept-area calculations
only as explicitly named `LEGACY_*` controls.  Corrected diagnostics, long
factorials, compatible non-C continuations, strict legacy continuations, and a
selected barrier/packet discrimination campaign were executed at 900 K and
seed 5101 from archived checkpoints or the matched clean initial state as
appropriate.  No multi-temperature or multi-seed campaign was launched.

## Source lineage

The work was performed on `fix/jerky-mechanism-integrity-20260820` without
merging `main`.  The principal pushed commits are:

| Commit | Purpose |
|---|---|
| `92816e2` | conserved GB-area-loss source, signed inventory, and GB/TJ climb |
| `124f639` | localized competing TJ defect sinks |
| `396be9e` | exact stage residence and long-run cadence |
| `837a91a` | ten-case barrier/packet discriminating configuration |
| `36d41c0` | stage, causal, occupancy, and spatial analysis extensions |
| `e3586c7` | reproducible execution/movie indexer |
| `c45ad06` | global-inventory release-state selection correction |
| `be2aeea` | factorial, crossover, and minimum-model comparisons |
| `361fccc` | conservation derivation and parameter record |

Run manifests append `-dirty` because preserved historical result directories
were already untracked; the tracked launch commit is recorded separately in
each manifest.  Raw simulation roots remain untracked under repository policy.

## Corrected climb physics and conservation

For unit out-of-plane thickness, the implemented excess volume is

\[
V_{ex}^{GB}(t)=\alpha\sum_\beta v_{ex,\beta}^{GB}L_\beta(t),
\]

and the released defect quota is

\[
\Delta N_{req}=\frac{\max[V_{ex}^{GB}(t)-V_{ex}^{GB}(t+\Delta t),0]}
{\Omega_{def}}.
\]

No swept area, absolute length change, or rigid-translation displacement is
used as a source.  The production configuration supplies
`excess_volume_per_area=0.01` directly and
`point_defect_formation_volume=0.02`, so one unit of total GB-length loss
releases 0.5 defect quota.  No density/thickness pair is supplied for these
runs.  The general optional material expression remains
`delta_GB * (1 - rho_GB/rho_lattice)`.

The material-wide signed ledger enforces

\[
N_{required}=N_{GB}+N_{TJ}+q_{stored}+N_{external}
\]

through domain split, merge, retirement, grain disappearance, accommodation,
and checkpoint/restart.  Production runs use no external sink.  `C_GB`
permits only GB-disconnection sinks.  `C_GBTJ` lets GB and TJ sinks compete;
TJ completion additionally requires correctly signed available flux,
three-boundary geometric compatibility, Burgers/disconnection compatibility
or the finite residual-energy model, and atomic exactly-once quota removal.

Frozen production parameters are: area-loss coefficient 1; free-volume
stiffness 0.05; climb trigger 0.25; GB and TJ release quota 1; stage barriers
0.45/0.55/0.65 eV; stage prefactors 100000; TJ compatibility tolerance 0.25;
TJ step length 1; Burgers magnitude 0.25; and residual stiffness 1 eV.
The complete derivation is in
[area_loss_climb_conservation_20260822.md](area_loss_climb_conservation_20260822.md).

## Regression and full-suite validation

The new tests cover translation with zero demand, analytic shortening demand,
zero release for length increase, split/merge and retirement conservation,
grain disappearance, exactly-once GB/TJ consumption, wrong-sign and
incompatible TJ rejection, finite Burgers residual behavior, exact restart,
and approximate physical grid invariance.

The final valid invocation used `/opt/anaconda3/bin/python` with
`PYTHONPATH=src` and passed **143/143 tests in 33.80 s**.  JUnit artifact:
`results/validation/area_loss_climb_full_tests_20260823.xml`; SHA-256:
`7856b28418fd056c30cdf211838890f4e4c423b826f5b53ce096dbf6f17b7559`.
An immediately preceding invocation without `PYTHONPATH=src` failed during
collection only; its JUnit output was overwritten and is not evidence.

## Canonical execution roots

| Campaign | Root | Source |
|---|---|---|
| Valid corrected diagnostic | `results/jerky_area_loss_diagnostic_20260822/20260823T002904Z-457cf9e522` | `124f6398cc0ce7bb11b02184f056a73a57bc6152-dirty` |
| Corrected long 16-case factorial | `results/jerky_area_loss_factorial_long_20260822/20260823T012141Z-b907865386` | `396be9e48fdf85279833de1004ecf2b9c74065a0-dirty` |
| Compatible non-C continuations | `results/jerky_long_non_c_20260822/20260823T012450Z-364b2906f9` | `396be9e48fdf85279833de1004ecf2b9c74065a0-dirty` |
| Swept-area legacy continuations | `results/jerky_long_legacy_20260822/20260823T012535Z-793b9644d0` | `396be9e48fdf85279833de1004ecf2b9c74065a0-dirty` |
| Barrier/packet discrimination | `results/jerky_area_loss_discriminating_20260822/20260823T102804Z-2ab688a139` | `361fcccfddcbe9829288775425066074b3607bce-dirty` |

The invalid interrupted diagnostic development trace
`results/jerky_area_loss_diagnostic_20260822/20260823T002055Z-385881403d`
is excluded from all conclusions.

Non-C continuations originate from
`results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0`.
Legacy continuations originate from
`results/jerky_factorial_900K/20260820T233721Z-f5d010440c` at step 3500 and
retain the old swept-area physics.  No corrected C run resumes a legacy C
checkpoint.

## Diagnostic integrity evidence

The six valid diagnostics are B0, G, C_GB, C_GBTJ, GTSC_GB, and GTSC_GBTJ.
Their executed analysis is in
`results/production_summaries/area_loss_diagnostic_20260822_v3`.

| Condition | Required | GB accommodated | TJ accommodated | Stored | max abs residual |
|---|---:|---:|---:|---:|---:|
| C_GB | 413.5 | 377.0 | 0 | 36.5 | 7.96e-13 |
| C_GBTJ | 401.5 | 157.0 | 205.226323 | 39.273677 | 1.36e-12 |
| GTSC_GB | 292.5 | 263.0 | 0 | 29.5 | 8.53e-13 |
| GTSC_GBTJ | 339.5 | 131.0 | 178.345927 | 30.154073 | 1.25e-12 |

The diagnostic C_GBTJ sink split is 43.34% GB / 56.66% TJ; GTSC_GBTJ is
42.35% GB / 57.65% TJ.  Pure C_GB has zero persisted G occupancy and zero TJ
accommodation.  Duplicate GB consumption residual is exactly zero and TJ
roundoff is approximately `-2e-13`.

Persisted C occupancy is 3.43% on C_GB GBs; C_GBTJ has 0.86% on GBs and 1.68%
on TJs.  Mean roughness quality remains approximately 1.78--1.79 across the
diagnostics, with no runaway waviness.  The v3 release-state calculation uses
the material-wide inventory, so `DeltaMu_Nv` is nonzero where physically
appropriate.

## Barrier/packet discrimination

The executed result tables are in
`results/production_summaries/area_loss_discriminating_20260823`; roughness is
in `area_loss_discriminating_20260823_roughness.csv`.

| Pair/case | End step | End grains | Outcome |
|---|---:|---:|---|
| G_LOW | 981 | 130 | target |
| G_HIGH | 1400 | 150 | censored |
| T_LOW | 890 | 130 | target |
| T_HIGH | 1400 | 139 | censored |
| C_GB_FAST | 874 | 130 | target |
| C_GB_SLOW | 1400 | 134 | censored |
| C_GBTJ_FAST | 877 | 130 | target |
| C_GBTJ_SLOW | 1006 | 130 | target |
| GS_PACKET_05 | 1400 | 152 | censored |
| GS_PACKET_2 | 1400 | 151 | censored |

C_GB_FAST required 579.5, accommodated 552 at GB sinks, and stored 27.5.
C_GB_SLOW required 610, accommodated 537, and stored 73.  C_GBTJ_FAST split
accommodation 60.80% GB / 39.20% TJ; C_GBTJ_SLOW split 33.73% GB / 66.27% TJ.
All maximum conservation residuals are no larger than `2.62e-12`.

The observed mean GB-sink cycle time changes from 1.430 (C_GB_FAST) to 5.098
(C_GB_SLOW).  For C_GBTJ, the observed means are 2.475/1.049 for fast GB/TJ
paths and 6.548/2.473 for slow GB/TJ paths.  Mean boundary roughness spans
1.782--1.820 across all ten conditions, so the discriminating extremes do not
produce runaway waviness.

## Final run and movie index

All 49 required conditions are terminal and have movies: 20 reached their
grain-count target and 29 are explicitly censored at the step ceiling.  There
are no duplicate run paths or missing movies.  The complete per-run record of
source SHA, starting and ending step/time/grain count, status, frame count, and
movie path is
[jerky_area_loss_movie_index_20260823.md](jerky_area_loss_movie_index_20260823.md).
The machine-readable equivalent is
`results/production_summaries/area_loss_execution_index_20260823.csv`.

The five campaign contact sheets are:

- `results/jerky_area_loss_diagnostic_20260822/20260823T002904Z-457cf9e522/contact_sheet.png`
- `results/jerky_area_loss_factorial_long_20260822/20260823T012141Z-b907865386/contact_sheet.png`
- `results/jerky_long_non_c_20260822/20260823T012450Z-364b2906f9/contact_sheet.png`
- `results/jerky_long_legacy_20260822/20260823T012535Z-793b9644d0/contact_sheet.png`
- `results/jerky_area_loss_discriminating_20260822/20260823T102804Z-2ab688a139/contact_sheet.png`

The corrected long factorial reached 70 grains for C_GB, GC_GB, TC_GB,
GTC_GB, C_GBTJ, GC_GBTJ, TC_GBTJ, and GTC_GBTJ.  The eight S-containing
conditions that did not reach 70 were retained at step 10000 and ended with
83--166 grains.  The compatible non-C continuations genuinely increased the
checkpoint time/step; B0, G, T, and GT reached 69--70 grains, while S, GS, TS,
GTS, and QIU were censored at step 10000.  All eight strict legacy controls
continued from step 3500/time 140 to step 10000/time 400.

## Corrected long-run sink partition and occupancy

The executed long-run sink table is part of
`results/production_summaries/area_loss_factorial_long_20260823`.  GB-only
conditions have `f_GB=1` and `f_TJ=0` exactly.  The competing-sink results are:

| Regime | Required | GB | TJ | Stored | f_GB | f_TJ | max abs residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| C_GBTJ | 1254.5 | 520 | 674.867 | 59.633 | 0.435 | 0.565 | 3.18e-12 |
| GC_GBTJ | 1321 | 688 | 585.623 | 47.377 | 0.540 | 0.460 | 4.09e-12 |
| TC_GBTJ | 1274.5 | 542 | 679.565 | 52.935 | 0.444 | 0.556 | 6.25e-12 |
| SC_GBTJ | 3129 | 1293 | 1831.716 | 4.284 | 0.414 | 0.586 | 3.55e-11 |
| GTC_GBTJ | 1325.5 | 688 | 598.490 | 39.010 | 0.535 | 0.465 | 4.55e-12 |
| GSC_GBTJ | 2100 | 1018 | 1077.698 | 4.302 | 0.486 | 0.514 | 5.25e-11 |
| TSC_GBTJ | 3601 | 1407.75 | 2187.552 | 5.698 | 0.392 | 0.608 | 2.01e-10 |
| GTSC_GBTJ | 2502 | 1186 | 1288.085 | 27.915 | 0.479 | 0.521 | 4.32e-11 |

These fractions use accommodated quota in the denominator; stored quota is
reported separately.  No external sink is active.

Direct occupancy was integrated from persisted mechanism-state flags rather
than reconstructed from release rows.  Representative all-entity fractions
are:

| Regime | G | T | C |
|---|---:|---:|---:|
| C_GB | 0 | 0 | 0.0362 |
| C_GBTJ | 0 | 0.0010 | 0.0121 |
| GTC_GB | 0.1103 | 0.0568 | 0.0176 |
| GSC_GBTJ | 0.0359 | 0.00004 | 0.0088 |
| GTSC_GB | 0.0313 | 0.0476 | 0.0110 |
| GTSC_GBTJ | 0.0323 | 0.00013 | 0.0100 |

The exact GB/TJ entity-resolved values are in `direct_occupancy.csv`.

## Common-window comparisons, crossover, and factorial effects

All comparisons use identical grain-count windows; the executed tables are in
`results/production_summaries/area_loss_comparisons_final_20260823`.
Over N=190 to 160, C_GB grows at 0.032910 versus 0.040376 for B0
(ratio 0.815) and 0.004708 for LEGACY_C (ratio 6.990).  GTSC_GBTJ grows at
0.018473 versus 0.014743 for GTS (ratio 1.253) and 0.004299 for LEGACY_GTSC
(ratio 4.297).  The corrected/legacy GTSC_GBTJ ordering reverses in the late
N=120 to 100 window, where the ratio is 0.544; this is the clearest measured
kinetic crossover and is recorded as such rather than hidden by a trajectory-
wide average.

The corrected within-sink-family factorial contrasts show:

| Sink model | Metric | G | T | S |
|---|---|---:|---:|---:|
| GB | N190--160 growth rate | +0.00545 | -0.00577 | -0.01491 |
| GB | N190--160 jerkiness CV | -1.30465 | -0.25287 | +1.42381 |
| GBTJ | N190--160 growth rate | +0.00265 | -0.00022 | -0.02218 |
| GBTJ | N190--160 jerkiness CV | -1.32159 | -0.03518 | +1.40403 |
| GBTJ | TJ sink fraction | -0.08897 | +0.00638 | +0.04585 |

Thus S is the dominant negative early-growth and positive-jerkiness main
effect, while G reduces jerkiness and partly restores growth.  This supersedes
the old swept-area factorial comparison.

Against the full four-mechanism reference, the revised complexity/fidelity
compromises are GTC_GB (standardized RMS 0.449 with three mechanisms) and
GSC_GBTJ (RMS 0.326 with three mechanisms).  C_GB and C_GBTJ remain
Pareto-optimal one-mechanism models, but no longer justify the earlier claim
that C alone reproduces GTSC.  Full Pareto membership and metrics are in
`revised_minimum_models.csv`.

## Executed expanded analyses

The final 16-case output directory contains `run_summary.csv`,
`stage_residence.csv`, `cycle_time_audit.csv`,
`local_event_response.csv`, `causal_nulls.csv`,
`kinetic_resistance_windows.csv`, `release_state_selection.csv`,
`direct_occupancy.csv`, `sink_partition.csv`,
`spatial_release_correlations.csv`, and `spatial_release_summary.csv`.
These are executed tables, not empty templates: they contain 96 stage rows,
32 path-cycle audits, 140951 local responses, 64 causal-null rows, 80
resistance fits, 39214 matched release/control states, 28 occupancy rows, and
542050 spatial event pairs.

### Stage-resolved kinetics

Representative observed mean residence times are:

| Regime/path | Nucleation | Exchange | Transport | Observed cycle | rate-mixture prediction | KS p |
|---|---:|---:|---:|---:|---:|---:|
| C_GB / GB | 0.00330 | 3.49e-5 | 4.23485 | 4.23824 | 2.00896 | 1.68e-3 |
| C_GBTJ / GB | 0.00331 | 6.33e-5 | 9.31078 | 9.31393 | 0.90560 | 6.58e-49 |
| C_GBTJ / TJ | 0.24578 | 1.58e-5 | 2.69618 | 2.93702 | 0.23496 | 2.12e-110 |
| GTSC_GBTJ / GB | 13.8996 | 6.18e-5 | 12.4610 | 19.5012 | 1.21830 | 3.16e-85 |
| GTSC_GBTJ / TJ | 1.03812 | 2.83e-5 | 3.37069 | 4.40322 | 0.52049 | 2.72e-231 |

The sink application itself is atomic immediately after transport, so its
measured residence `t_sink` is zero rather than an unconfigured fourth
stochastic wait.  All 24 exercised path distributions reject the simple
conditional three-stage hypoexponential mixture at p<=0.0017.  The observed
residence includes competition, repinning, geometry gating, and topology
changes; the bare rate sum is therefore not an adequate cycle model for these
coupled trajectories.

### Local response and causal nulls

Release-aligned persistent-domain motion is nonzero and sink-specific.  For
C_GB GB completions, mean signed/absolute velocity changes at 1 frame are
-0.0134/+0.0789 and at 8 frames are +0.0033/+0.0458.  For C_GBTJ TJ
completions they are -0.0108/+0.0421 and -0.0070/+0.0700.  GTSC_GBTJ TJ
completions have small signed changes (-0.00157 at 1 frame) but positive
absolute response (+0.00737); all event-level signed displacement and
1/2/4/8-frame rows are retained in `local_event_response.csv`.

The stronger circular-shift and block-shuffle nulls do **not** show an
aggregate growth-rate excess at the saved-frame cadence.  Release indicators
occupy 57--201 output frames per regime and, in the long ceiling runs, nearly
every frame.  Across all cases/windows, maximum absolute circular and block
excesses are only 2.81e-6 and 6.39e-6.  Thus these data establish local
release-aligned responses but do not establish a frame-resolved causal burst
excess; a future campaign needs finer event-adjacent output for that test.

### Resistance and release-state selection

Successive topology-window fits of `1/Rdot = a R + b` are weak: median R-squared
is 0.0717 and only 2.5% of windows reach 0.5.  Several selected cases reverse
the sign of `a` from early to late topology (for example C_GB, +361.9 to
-664.7, and C_GBTJ, +34.5 to -1519.1).  A single trajectory-wide linear
resistance law is therefore not supported.

Matched release-minus-pinned mean work differences are positive for the
minimum/reference pairs: C_GB has `pVn=+0.0480`,
`DeltaMu_Nv=W_total=+0.1465`; C_GBTJ has +0.0184/+0.2292; GTC_GB has
+0.0643/+0.2291; GSC_GBTJ has +0.0371/+0.00531; GTSC_GB has
+0.0589/+0.0909; and GTSC_GBTJ has +0.0317/+0.0383.  The corresponding shear
differences are zero or below 0.001 in magnitude for these cases.  The full
release and matched-pinned rows are retained rather than reducing the result
to means.

### Spatial release correlations

The spatial table resolves time lag, periodic Euclidean distance, approximate
GB arclength, same entity/GB, shared-TJ connectivity, and sink-path pairing.
At lag <=0.16 and distance <=8, shared-TJ fractions are 0.54--0.81 for the six
selected cases, versus 0.0005--0.0074 beyond distance 32.  For example,
GTSC_GBTJ gives 0.718 in 443 near pairs versus 0.00066 in 7621 far pairs.
Same-GB arclength pairs at this short lag are almost entirely repeated events
on the same persistent domain; there are no populated >24 arclength pairs, so
the present archive cannot support a separate long-arclength branching claim.
These correlations are descriptive connectivity evidence, not a randomized
spatial-causality test.

The roughness analyses cover every saved frame.  Mean boundary quality spans
1.748--1.833 in the corrected long factorial, 1.750--1.823 in compatible
non-C continuations, 1.835--1.846 in legacy continuations, and 1.782--1.820 in
the discriminating campaign.  No run shows runaway boundary waviness.

## Recommendation for the next temperature/seed campaign

Do not immediately launch a broad four-temperature matrix.  The next staged
screen should compare the revised minimum model with its full reference inside
each sink architecture: GTC_GB versus GTSC_GB, and GSC_GBTJ versus GTSC_GBTJ.
Run 800 and 1100 K first with three matched seeds (24 conditions).  If the two
minimum/reference pairs retain their ordering and sink partition, add 900 K
with the same seeds (12 conditions; the exact compatible seed-5101 results may
be reused with provenance), making a 36-condition inference set.  Add B0 at
the extremes only if intrinsic mobility is made temperature-dependent.

This design estimates seed uncertainty, directly tests whether the revised
minimum models remain adequate away from 900 K, and preserves separate GB and
GBTJ interpretations.  It should use N=190--160 as the mandatory common
window, retain later windows for crossover detection, and require the same
conservation, occupancy, stage, event-null, roughness, and movie checks used
here.
