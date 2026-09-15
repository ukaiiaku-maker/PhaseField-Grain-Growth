# 900 K shear-scale validation execution report

## Outcome

Phase A is complete.  All ten matched 900 K, seed-5101 shear-screen runs are
terminal, all ten diagnostic movies exist, event-triggered solver-step traces
were analyzed, conservation and morphology checks pass, and the full repository
suite passes 158/158 tests.

The decision is **NO-GO for the 24-run temperature campaign**.  No 800 K,
1100 K, or additional 900 K seed calculation was launched.  The strongest
shear effect occurs only at `shear_stiffness=0.40`, where every tested family
ceiling-censors at 155--157 grains after 10,000 steps.  In addition, maximum
stored shear energy is non-monotone between 0.20 and 0.40 in both full GTSC
families.  This fails the directive's requirement that the large S effect not
be a one-value pathological-stagnation result, even though the early response
curves are graded, the reduced/full ordering is stable, and 0.20 is not an
abrupt threshold.

## Provenance and validation

The work was performed on `exp/shear-temperature-screen-20260823` without
merging `main`.

| SHA | Purpose |
|---|---|
| `7253fc0a92db9114e1b5443ff0c22daf5e5ab343` | frozen completed 900 K reference |
| `f23aa55a84e51723890e669b187f7712f452aae6` | density/resistance corrections and invariant event tracing |
| `dbca0609e16e196f37bc3fd03346252fe9faac16` | restart-safe Zstd Parquet event traces; exact simulation SHA |
| `1faa4f7b5a67c773b814dfec038299ef9cecbf8d` | shear-screen analysis, recovery, and gate implementation |

The final full suite passed `158 passed in 47.81s`.  A separate deterministic
event-recorder replay passed trajectory identity, exact conservation, zero TJ
accommodation for `C_GB`, competing GB/TJ sinks for `C_GBTJ`, no hidden G route
in pure C, and no duplicate quota consumption.  Its report is
`/Users/sdillon/PF-graingrowth-output/event_trace_replay_parquet_20260823/replay_report.json`.

The valid campaign root is:

```text
results/jerky_shear_scale_900K_20260823/20260823T153124Z-39f3ffde01
```

An earlier storage-diagnostic partial root is preserved but excluded from all
analysis:

```text
results/jerky_shear_scale_900K_20260823/20260823T152452Z-e26005f341
```

## Physics corrections retained

The optional density-derived GB excess volume is now

\[
v_{ex}^{GB}=\delta_{GB}\left(1-\frac{\rho_{GB}}{\rho_{lattice}}\right).
\]

All completed production calculations, including this screen, bypass the
density branch by supplying `excess_volume_per_area=0.01` directly.  They use
`point_defect_formation_volume=0.02`, so one unit of total GB-length loss
creates 0.5 defect quota.  Therefore neither the frozen 49-run evidence set nor
this campaign requires a physics rerun for the density-expression correction.

The diagnostic local resistance law is correctly labeled

\[
\frac{1}{\dot R}=aR+b.
\]

It remains a diagnostic only; the five-window fits in
`analysis/area_loss/kinetic_resistance_windows.csv` are generally weak and are
not used to force a kinetic law.

The climb ledger enforces

\[
N_{required}=N_{GB}+N_{TJ}+q_{stored}+N_{external}.
\]

No external sink is active.  The full derivation and frozen material parameters
remain in [area_loss_climb_conservation_20260822.md](area_loss_climb_conservation_20260822.md).

## Terminal run and movie index

All runs start at solver step 0, simulation time 0, and 200 grains.  A terminal
status at 10,000 steps is a valid ceiling-censored completion, not target
completion.

| Regime | End step | End time | End grains | Status | Movie |
|---|---:|---:|---:|---|---|
| `GTC_GB` | 1728 | 69.12 | 119 | completed | `GTC_GB-T900-s5101/microstructure.gif` |
| `GTSC_GB_Ks010` | 2265 | 90.60 | 120 | completed | `GTSC_GB_Ks010-T900-s5101/microstructure.gif` |
| `GTSC_GB_Ks020` | 3561 | 142.44 | 120 | completed | `GTSC_GB_Ks020-T900-s5101/microstructure.gif` |
| `GTSC_GB_Ks040` | 10000 | 400.00 | 155 | completed/censored | `GTSC_GB_Ks040-T900-s5101/microstructure.gif` |
| `GSC_GBTJ_Ks010` | 2049 | 81.96 | 120 | completed | `GSC_GBTJ_Ks010-T900-s5101/microstructure.gif` |
| `GSC_GBTJ_Ks020` | 3656 | 146.24 | 120 | completed | `GSC_GBTJ_Ks020-T900-s5101/microstructure.gif` |
| `GSC_GBTJ_Ks040` | 10000 | 400.00 | 156 | completed/censored | `GSC_GBTJ_Ks040-T900-s5101/microstructure.gif` |
| `GTSC_GBTJ_Ks010` | 1979 | 79.16 | 120 | completed | `GTSC_GBTJ_Ks010-T900-s5101/microstructure.gif` |
| `GTSC_GBTJ_Ks020` | 3763 | 150.52 | 120 | completed | `GTSC_GBTJ_Ks020-T900-s5101/microstructure.gif` |
| `GTSC_GBTJ_Ks040` | 10000 | 400.00 | 157 | completed/censored | `GTSC_GBTJ_Ks040-T900-s5101/microstructure.gif` |

Movie paths above are relative to the valid campaign root.  This environment
has Pillow but no ffmpeg writer, so the final animated six-panel diagnostics
are GIFs.  The machine-readable index is
[shear_scale_movie_index_20260823.csv](shear_scale_movie_index_20260823.csv).
The final-state contact sheet is
`analysis/shear_contact_sheet.png`.

## Matched topology-window kinetics

The mandatory 190-to-160-grain results are:

| Regime | Rdot | CV jerk | stationary fraction | fG | fT | fC |
|---|---:|---:|---:|---:|---:|---:|
| `GTC_GB` | 0.017091 | 1.446 | 0.071 | 0.750 | 0.125 | 0.125 |
| `GTSC_GB_Ks010` | 0.017044 | 1.206 | 0.133 | 0.719 | 0.140 | 0.141 |
| `GTSC_GB_Ks020` | 0.015465 | 1.334 | 0.000 | 0.661 | 0.187 | 0.152 |
| `GTSC_GB_Ks040` | 0.004677 | 2.051 | 0.000 | 0.452 | 0.322 | 0.226 |
| `GSC_GBTJ_Ks010` | 0.022478 | 1.084 | 0.000 | 0.739 | 0.000 | 0.261 |
| `GSC_GBTJ_Ks020` | 0.018095 | 1.377 | 0.000 | 0.672 | 0.000 | 0.328 |
| `GSC_GBTJ_Ks040` | 0.004604 | 2.382 | 0.136 | 0.369 | 0.000 | 0.631 |
| `GTSC_GBTJ_Ks010` | 0.024056 | 0.976 | 0.083 | 0.706 | 0.046 | 0.247 |
| `GTSC_GBTJ_Ks020` | 0.018473 | 0.864 | 0.000 | 0.644 | 0.047 | 0.309 |
| `GTSC_GBTJ_Ks040` | 0.009565 | 1.510 | 0.028 | 0.415 | 0.098 | 0.487 |

Reduced/full ordering is stable at all three stiffnesses.  For GB sinks,
`Rdot(GTSC_GB)-Rdot(GTC_GB)` is -0.000047, -0.001626, and -0.012413.
For competing sinks, `Rdot(GTSC_GBTJ)-Rdot(GSC_GBTJ)` is +0.001579,
+0.000378, and +0.004961.  Thus model identity does not oscillate.

The late 160-to-140 window exposes the pathological branch: Rdot falls to
0.000332 (`GTSC_GB_Ks040`), 0.000322 (`GSC_GBTJ_Ks040`), and 0.000036
(`GTSC_GBTJ_Ks040`).  The 0.10/0.20 cases reach 120 grains, whereas none of the
0.40 cases reaches 155.  Full window tables are in
`analysis/shear_screen/topology_window_metrics.csv`.

## Shear mechanism audit

Across 0.10, 0.20, and 0.40, event-trace RMS shear state decreases smoothly
from about 0.209 to 0.196 to 0.143, while p95 absolute internal stress rises
from about 0.045 to 0.082--0.090 to 0.114--0.118.  The p95 normalized
instantaneous shear work also rises monotonically:

| Family | p95 abs(tau Vtau)/kBT at 0.10 / 0.20 / 0.40 | events changed >2x at 0.10 / 0.20 / 0.40 |
|---|---|---|
| `GTSC_GB` | 0.244 / 0.423 / 0.565 | 0.17% / 1.77% / 3.44% |
| `GSC_GBTJ` | 0.266 / 0.490 / 0.817 | 0.06% / 2.31% / 6.94% |
| `GTSC_GBTJ` | 0.258 / 0.458 / 0.713 | 0.09% / 2.12% / 5.32% |

Maximum total stored shear energy is 2.266/3.775/3.506 for `GTSC_GB` and
2.202/3.709/3.532 for `GTSC_GBTJ`; the 0.40 value is lower than 0.20 despite
higher stress/work.  This is consistent with fewer active migrating boundaries
during the plateau, but it means the constitutive scale is not yet calibrated
well enough to treat total stored energy as a monotone material response.

The evidence therefore favors S acting mainly through prolonged incompatible
state residence and accessible-mode suppression, with secondary changes in
sink competition, rather than through ubiquitous large instantaneous
`tau*Vtau` spikes.  Even at 0.40, only 3.4--6.9% of events change by more than
2x from instantaneous shear work and essentially none change by 5x or 10x.

## Event-triggered response and causal nulls

Every run contains merged, unique solver-step traces: 0.866--1.244 million
entity-step rows and 4,250--5,821 indexed events per run, with zero duplicate
entity-step rows.  The output evaluates 1, 2, 4, 8, 16, and 32-step windows
against circular-shift, block-shuffle, same-entity, and same-topology-window
nulls.  All 60 run/window rows and all four null families are complete.

Here `absolute_velocity_response` means the signed change
`abs(v_post)-abs(v_pre)`, not `abs(v_post-v_pre)`; negative values are valid
post-event deceleration.  Most regimes show strong positive one-step and
eight-step enrichment against all nulls.  Responses decay by 32 steps.  The
notable exception is `GSC_GBTJ_Ks010` at one step, whose mean change is -0.555;
its circular/block/same-entity controls are still more negative, while the
same-topology control is near zero.  Consequently the solver-step traces prove
event-linked, heterogeneous local response but do not justify claiming that a
large instantaneous shear-work spike directly causes every release.

Detailed event rows and null statistics are in
`analysis/shear_screen/solver_step_event_responses.csv` and
`analysis/shear_screen/solver_step_causal_nulls.csv`.  The older ordinary-frame
nulls remain saturated or nearly null and are retained separately in
`analysis/area_loss/causal_nulls.csv`.

## Sink partition, occupancy, stage kinetics, and conservation

GB-only runs have exactly zero TJ accommodation.  Competing-sink cases use TJ
sinks substantially:

| Regime | GB sink fraction | TJ sink fraction | final stored quota |
|---|---:|---:|---:|
| `GSC_GBTJ_Ks010` | 0.438 | 0.562 | 45.956 |
| `GSC_GBTJ_Ks020` | 0.435 | 0.565 | 35.749 |
| `GSC_GBTJ_Ks040` | 0.369 | 0.631 | 0.000 |
| `GTSC_GBTJ_Ks010` | 0.395 | 0.605 | 35.778 |
| `GTSC_GBTJ_Ks020` | 0.443 | 0.557 | 31.417 |
| `GTSC_GBTJ_Ks040` | 0.398 | 0.602 | 1.375 |

Direct occupancy falls sharply as stiffness increases.  For example,
`GTSC_GB` G occupancy is 0.0530, 0.0297, and 0.00534; its C occupancy is
0.0141, 0.0125, and 0.00552.  In `GTSC_GBTJ`, G occupancy is 0.0608, 0.0319,
and 0.00486 and C occupancy is 0.00729, 0.00985, and 0.00382.  This supports
mode-population/state-duration control rather than a pure instantaneous-work
mechanism.

Observed stage residence strongly rejects bare independent-stage kinetics.
GB nucleation residence is roughly 5.2--11.7 time units and GB transport
2.4--11.3, versus milliscale or smaller exchange residence except in the 0.40
GBTJ cases.  TJ nucleation is 0.19--0.72 and TJ transport 2.11--2.98.  The
hypoexponential KS p-values are extremely small, confirming that geometry,
competition, repinning, and topology govern cycle timing.

The largest absolute material conservation residual across the screen is
`3.14e-11`, comfortably below the gate tolerance `1e-9`.  There is no hidden
sink consumption, duplicate quota use, or forbidden TJ accommodation.

## Morphology

Visual review of all ten movies and the contact sheet finds bounded polygonal
grains without progressive numerical waviness.  Mean pixel isoperimetric
quality lies in 1.767--1.815.  Cross-stiffness family ratios are at most 1.0242
and maximum progressive ratios at most 1.0360.  High-wave-number boundary power
is similarly narrow, about 0.2895--0.2960.  Morphology passes the gate.

## Gate and next campaign recommendation

The automated gate passes 11 of 12 checks and fails
`stress_and_energy_scale_continuously` because the full-family maximum stored
energy is non-monotone.  Scientific review retains **NO-GO** because all three
0.40 cases also share severe late stagnation; the strongest S effect is
therefore inseparable from a one-value censored branch.

The next campaign should remain at 900 K and seed 5101 and localize the
constitutive crossover with `Ks={0.25,0.30,0.35}` for `GTSC_GB`, `GSC_GBTJ`,
and `GTSC_GBTJ`, reusing `GTC_GB` as the no-S control.  Add stored shear energy
normalized per active GB length/entity and compare it in matched topology
windows, not only by global trajectory maximum.  If a smooth admissible band is
found, confirm the selected value with at least seeds 5102 and 5103 before the
24-run 800/1100 K matrix.  Do not vary beta simultaneously; the optional
two-point beta check remains secondary.

Accordingly, no temperature roots or 900 K multi-seed completion roots exist
for this phase.  The prior corrected and legacy datasets remain untouched.

## Artifact index

Within the valid campaign root:

- `analysis/shear_gate/shear_gate.json` and `.md`: automated gate.
- `analysis/shear_screen/`: shear, topology-window, event-response, null, and plot outputs.
- `analysis/area_loss/`: sink, occupancy, stage, cycle, resistance, spatial, and release-state tables.
- `analysis/frame_roughness.csv`: all-frame morphology audit.
- `analysis/boundary_spectral_roughness.csv`, `_all_frames.csv`, and `.md`: spectral audit.
- `analysis/shear_contact_sheet.png`: ten-run final-state sheet.
- `video_manifest.json`: ten terminal run paths/statuses.
- `recovery_manifest.json`: completed recovery inventory.

