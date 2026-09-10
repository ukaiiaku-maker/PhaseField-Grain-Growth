# Phase-1 long-time 900 K: initial integrated audit and scientific report

## Executive decision

The initial 100,000-step Phase-1 matrix is operationally complete and locally
integrated: all nine required regimes are present in the canonical campaign
root, all nine passed checkpoint/manifest/config/frame/morphology checks, all
four climb-bearing runs close the material inventory, and all nine movies were
reconstructed.  No simulation is currently running.

The scientific gate is mixed:

- `B0`, `G`, `T`, `GT`, and `GTC_GB` reached 100 grains and show finite
  cumulative growth exponents near 1.5--1.8.  None supports a cubic law; local
  terminal-window exponents are noisy, but their fixed-law coefficients do not
  show progressive collapse.
- `GSC_GBTJ_Ks025`, `GTSC_GBTJ_Ks025`, and
  `GTSC_GBTJ_Ks030_LONG` show progressive slowing.  Their fitted exponent rises
  while `K2`, `K3`, and normalized growth rate fall.  This is not merely a
  smaller constant coefficient.
- `QIU` undergoes a late, non-self-similar avalanche from 494 grains at step
  9,800 to 99 grains at step 10,246, an interval of 446 solver steps.  Its
  terminal mean compactness is 4.37,
  versus about 1.34--1.36 for the other eight cases.  It is a transient/
  morphology diagnostic, not a converged growth-law result.

Consequently, the matrix is complete as an initial ceiling evaluation, but the
scientific campaign should not yet be called final.  Continue the two active
0.25 shear/climb cases to 150,000 steps and diagnose QIU before launching a
temperature or seed campaign.

## Canonical provenance and integration

- Campaign root:
  `results/long_time_kinetics_900K_20260824/20260825T012949Z-6d83ee5c82`
- Frozen production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`
- Corrected analysis source: `2785ff62c13b1c9b36338fb99ca50f29047399de`
- Initial-state generation source: `f75e44c40f59b57174aa988c51f13a4e7b870fcf`
- Initial NPZ: `results/initial_conditions/seed-5101-5ba2f47025f7d22a.npz`
- Initial NPZ SHA-256:
  `106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`
- Initial metadata SHA-256:
  `ff66cb5f938b5b70fe9e611674956f1e45f30f0979714507a853750f61e1a9bb`
- Initial state: 384 by 384, periodic, 800 grains, seed 5101, 900 K.
- Campaign manifest: `campaign_manifest.json` (`status=completed`).
- Checksummed integration manifest: `integration_manifest.json`.
- External ownership ledger:
  `/Users/sdillon/HPC3/phase1-state/ownership.json`, revision 153 at finalization;
  every record is terminal `integrated` with audit status `passed`.
- Audit/integration implementation commit in the external orchestrator:
  `36a0ec6` (`phase1-orchestrator`, 24/24 tests passed).

The integration manifest records an essential-file SHA-256 inventory for each
run plus an aggregate name/content digest for every source-frame sequence.
Original HPC3 retrieval archives remain preserved outside the scientific
repository.

## GTC bookkeeping reconciliation

HPC3 runner run `20260827T150438Z-fd59673-ed038d`, Slurm job `55622081`, ended
with `application_exit=1`, `finalization_exit=0`, and `complete=false`.  The only
terminal wrapper error was `production configuration checksum mismatch`.  The
cause was operational: HPC3 rewrote `initial_state_file` to a node-local
absolute path, while the old wrapper compared that runtime hash with the
campaign-relative production hash.

The scientific result is valid:

- retrieved archive SHA-256 expected and actual:
  `b83726aee8e8fee107c6e5dc74cc81bbdb453dc7fc59041703604ec9037c2b6e`;
- manifest `status=completed`, no scientific failure;
- source SHA exactly `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`;
- normalized production config SHA exactly
  `45419cc85f3d893122401a4f0f757d19818ff17782a1b82cbd734adc39a93b39`;
- runtime config SHA exactly matches the untouched run manifest;
- checkpoint and manifest agree at step 12,868, time 514.72, 100 grains;
- all 55 frames span step 0 through 12,868 with maximum step gap 600;
- maximum conservation residual is `2.7831e-10`.

It is therefore classified as
`scientific_success_wrapper_bookkeeping_failure`, not rerun or silently
relabeled.  Evidence is in
`/Users/sdillon/HPC3/phase1-state/evidence/formal-gtc-wrapper-reconciliation-20260910T175844Z.json`.

## Formal run status

| Regime | Owner | End step/time | N | G/G0 | Frames | Scientific classification |
|---|---|---:|---:|---:|---:|---|
| B0 | local | 8,200 / 328.00 | 100 | 2.828 | 37 | apparently converged finite regime |
| GTSC_GBTJ_Ks025 | local | 100,000 / 4000.00 | 333 | 1.550 | 122 | progressive slowing; continue |
| GTSC_GBTJ_Ks030_LONG | local | 100,000 / 4000.00 | 406 | 1.404 | 114 | progressive slowing; long plateau |
| GTC_GB | HPC3 | 12,868 / 514.72 | 100 | 2.828 | 55 | apparently converged finite regime |
| GSC_GBTJ_Ks025 | HPC3 | 100,000 / 4000.00 | 312 | 1.601 | 119 | progressive slowing; continue |
| G | HPC3 | 11,250 / 450.00 | 100 | 2.828 | 52 | apparently converged finite regime |
| T | HPC3 | 8,865 / 354.60 | 100 | 2.828 | 44 | apparently converged finite regime |
| GT | HPC3 | 11,751 / 470.04 | 100 | 2.828 | 53 | apparently converged finite regime |
| QIU | HPC3 | 10,246 / 409.84 | 99 | 2.843 | 21 | non-self-similar avalanche transient |

The complete machine-readable classification, including last grain-count
change, latest/cumulative fits, compactness, distribution distance, and next
action, is `analysis/phase1_regime_classification.csv`.

## Full test results

The authoritative invocation was:

```text
PYTHONPATH=src:scripts /opt/anaconda3/bin/python -m pytest -vv \
  --junitxml=results/validation/phase1_full_tests_20260910.xml
```

Result: **178 passed, 0 failed, 0 errors, 0 skipped in 53.00 s**.  JUnit
SHA-256:
`5ffdee46f27e159073ea53775b55fe61cc6550ab57e62def9e1a1adef223ce73`.

An immediately preceding `PYTHONPATH=src` invocation stopped during collection
because `tests/unit/test_shear_transition_analysis.py` imports a sibling module
from `scripts`; it executed no tests and its JUnit file was overwritten by the
authoritative passing run.  The external integration orchestrator separately
passes 24/24 tests in 0.48 s.

## Corrected climb conservation derivation and parameters

For unit out-of-plane thickness, the material excess volume is

\[
V_{ex}^{GB}(t)=\alpha\sum_\beta v_{ex,\beta}^{GB}L_\beta(t),\qquad
v_{ex}^{GB}=\delta_{GB}\left(1-\rho_{GB}/\rho_{lattice}\right).
\]

The direct calibration takes precedence in production.  The released signed
vacancy-equivalent quota is

\[
\Delta N_{req}^{n+1}=
\frac{\max[V_{ex}^{GB,n}-V_{ex}^{GB,n+1},0]}{\Omega_{def}}.
\]

A translating straight boundary creates no demand; only loss of total physical
GB length does.  The global, species-resolved invariant is

\[
N_{required,s}=N_{GB,s}+N_{TJ,s}+N_{external,s}
+q_{active,s}+q_{retired,s}.
\]

Split, merge, retirement, sink transfer, and checkpoint/restart preserve this
identity.  Same-sign availability, continuous completion-time ordering, and
atomic consumption prevent double use of quota.  The detailed proof is
[area_loss_climb_conservation_20260822.md](area_loss_climb_conservation_20260822.md).

Production uses `excess_volume_per_area=0.01`,
`point_defect_formation_volume=0.02`, and area-loss coefficient 1, hence 0.5
quota per unit GB-length loss.  No density pair is substituted.  Other frozen
values are: climb trigger 0.25; GB/TJ release quota 1; free-volume stiffness
0.05; GB and TJ serial barriers 0.45/0.55/0.65 eV; every stage prefactor
`1e5`; TJ tolerance 0.25; TJ step length 1; Burgers magnitude 0.25; residual
stiffness 1 eV; event-domain length 12; encounter density 0.20.

## Conservation and sink partitioning

| Regime | Required | GB accommodated | TJ accommodated | Stored | GB/TJ fraction | Max residual |
|---|---:|---:|---:|---:|---:|---:|
| GTC_GB | 7,599.0 | 7,345.0 | 0 | 254.0 | 1.000 / 0.000 | 2.78e-10 |
| GSC_GBTJ_Ks025 | 15,238.5 | 6,548.0 | 8,445.90 | 244.60 | 0.437 / 0.563 | 2.27e-9 |
| GTSC_GBTJ_Ks025 | 11,482.5 | 3,946.19 | 7,501.57 | 34.74 | 0.345 / 0.655 | 1.45e-9 |
| GTSC_GBTJ_Ks030_LONG | 9,218.5 | 3,582.97 | 5,635.03 | 0.50 | 0.389 / 0.611 | 2.12e-9 |

All are far inside the formal `1e-6` integration tolerance.  The fractions are
cumulative accommodated-quota shares; the zero terminal frame sink flags mean
no sink fired in that particular frame and must not be confused with cumulative
partition.  Full values are in `analysis/phase1_sink_partition.csv`.

## Long-time kinetics and crossover

At common progress `x=G/G0=1.15--1.40`, the population-size secant rates are:

| Regime | dG/dt | K2 | Elapsed time across window |
|---|---:|---:|---:|
| B0 | 0.11790 | 4.606 | 32.49 |
| G | 0.08413 | 3.287 | 45.52 |
| T | 0.10389 | 4.058 | 36.86 |
| GT | 0.07508 | 2.933 | 51.01 |
| GTC_GB | 0.06991 | 2.731 | 54.78 |
| GSC_GBTJ_Ks025 | 0.01114 | 0.435 | 343.73 |
| GTSC_GBTJ_Ks025 | 0.01298 | 0.507 | 295.11 |
| GTSC_GBTJ_Ks030_LONG | 0.001724 | 0.0674 | 2,220.95 |
| QIU | 0.01055 | 0.412 | 362.95 |

For the finite-growth group, cumulative profile fits at the terminal state give
`n=1.8` (B0), `1.5` (G), `1.5` (T), `1.5` (GT), and `1.6` (GTC_GB).  The moving
late-window optimum is less stable because only about 100--200 grains remain,
so the report gives fixed `K2/K3` and profile width equal weight rather than
overinterpreting a single local `n`.

The shear/climb group behaves differently:

- `GTSC_GBTJ_Ks025`: moving `n=7.2 -> 9.9 -> 30.0`; `K2=1.55 -> 0.565 -> 0.059`.
- `GSC_GBTJ_Ks025`: moving `n=7.3 -> 11.0 -> 23.2 -> 45.1`;
  `K2=1.54 -> 0.542 -> 0.072 -> 0.0266`.
- `GTSC_GBTJ_Ks030_LONG`: moving `n=10.6 -> 27.3 -> >=50`;
  `K2=1.16 -> 0.0888 -> 0.0252`; grain count is unchanged for the final
  42,200 steps.

Thus the earlier high exponents in the shear branch were not merely early
transients that recovered to `n=2--3`.  Over the available range they grow and
the coefficient collapses.  Increasing `Ks` from 0.25 to 0.30 reduces the
common-window rate by 86.7% and lengthens the crossing time by a factor 7.53.
This is a **shear-stiffness/kinetic-regime crossover**, not a temperature
crossover.  Only one temperature has been run, so no corrected crossover
temperature or apparent activation energy can be inferred.

## Direct occupancy and mechanism residence

The direct persisted-state time-sample fractions are in
`analysis/phase1_direct_occupancy.csv`.  Representative GB/TJ values are:

- `G`: GB G-pending 0.1187.
- `T`: TJ T-pending 0.07375.
- `GT`: GB G-pending 0.1131; TJ T-pending 0.05560.
- `GTC_GB`: GB G/C pending 0.1054/0.01934; TJ T-pending 0.05161.
- `GSC_GBTJ_Ks025`: GB G/C 0.00499/0.00222; TJ C 0.00581.
- `GTSC_GBTJ_Ks025`: GB G/C 0.00445/0.00112; TJ T/C
  0.00023/0.00327.
- `GTSC_GBTJ_Ks030_LONG`: GB G/C 0.00274/0.00095; TJ T/C
  0.00011/0.00269.

These low pending-state fractions do not contradict the strong slowdown:
continuous shear backstress acts outside discrete pending intervals.  In the
latest available windows, the slow cases have 88--94% stationary trajectory
fractions and 56--87% top-5% motion concentration.  Approximately 5% of active
GB length has `chi_s>0.8` in their terminal frames.  Joint G/T/C overlap is
reported explicitly in `analysis/mechanism_joint_occupancy_vs_grain_size.csv`;
independent residence fractions are not forced to sum to one.

## Stage-resolved kinetics

Exact stage records show broad, heavy-tailed residence distributions.  Median
GB nucleation/transport times are:

| Regime | GB nucleation | GB transport | TJ nucleation | TJ transport |
|---|---:|---:|---:|---:|
| GTC_GB | 5.522 | 0.755 | -- | -- |
| GSC_GBTJ_Ks025 | 0.00359 | 2.517 | 0.01988 | 0.1463 |
| GTSC_GBTJ_Ks025 | 0.00473 | 3.582 | 0.01431 | 0.1472 |
| GTSC_GBTJ_Ks030_LONG | 0.00433 | 3.052 | 0.01576 | 0.1505 |

Exchange is normally completed within the recording resolution.  Means and
upper tails are much larger than medians (for example full-Ks0.25 GB transport
mean 47.0 and p99 798.3), demonstrating competition, repinning, geometry
gating, and topology censoring rather than a simple sum of three independent
exponential means.  Full stage counts, means, medians, p90, and p99 are in
`analysis/phase1_stage_resolved_kinetics.csv`; stage occupancy versus grain size
is in `analysis/mechanism_residence_vs_grain_size.csv`.

## Event-triggered response and causal null

The corrected event-trace join accounts for repeated domain-local event IDs by
using each event's saved trace interval and nearest event step.  All aligned
rows now lie in the specified `-25..+50` solver-step window.

Weighted post-minus-pre absolute local-velocity responses are positive in all
available regimes: GB/TJ responses are 0.0306/0.1056 for full Ks0.25,
0.0224/0.0793 for full Ks0.30, 0.0354/0.1143 for GSC Ks0.25,
0.0694/0.2323 for GTC, and 0.0764/0.2536 for GT.

The new causal null retains each sampled event neighborhood and local velocity
history, but moves the alignment center outside +/-5 steps within the same
trace.  With 200 deterministic shuffles, 142/163 regime/window/event-family/
entity rows have positive excess at one-sided `p<0.05`.  Restricting to rows
with at least 20 usable events, weighted excess is 0.0666 (full Ks0.25), 0.0627
(full Ks0.30), 0.0778 (GSC Ks0.25), 0.1393 (GTC), 0.0789 (G), 0.2688 (T), and
0.1549 (GT).

This resolves, rather than contradicts, the prior coarse-frame null.  The old
ordinary-frame circular/block shuffle found no aggregate excess because
release flags occupied nearly every saved frame.  The present fine solver-step
trace supports temporal release alignment, but remains conditioned on sampled
event neighborhoods; it does not compare event-bearing neighborhoods with the
bulk.  Tables are `analysis/event_triggered_response_summary.csv` and
`analysis/causal_nulls_vs_grain_size.csv`.

## Distribution collapse, morphology, and jerkiness

For the five finite-growth cases, the last successive normalized equivalent-
diameter KS distance is about 0.049--0.071.  The three slow shear/climb cases
are also only 0.065--0.080, so their kinetic slowdown is not accompanied by an
obvious runaway size-distribution change over the limited `x=1.4--1.6` range.
That supports a constitutive/backstress origin, while the limited range still
precludes an asymptotic statement.

QIU is the exception.  Its final two outputs jump through a large topology
change, the normalized distributions do not evolve smoothly, and terminal
mean/max compactness are 4.37/7.65.  The QIU movie and contact sheet make the
avalanche directly inspectable; its terminal `N=99` must not be treated as a
successful conventional-growth endpoint.

Jerkiness remains separate from mean kinetics.  In the latest available slow
windows, event Fano factors are 2.08--7.13, top-5% motion concentration is
0.561--0.869, and stationary fraction is 0.885--0.940.  The finite-growth
cases can still show release-centered motion without progressive coefficient
collapse.

## Revised factorial interpretation and minimum model

At common `x=1.15--1.40`, G lowers `dG/dt` relative to B0 by 0.03376 and T by
0.01400; their interaction restores 0.00495.  Adding GB climb to GT lowers the
rate by only 0.00517.  In contrast, the combined switch from `GTC_GB` to full
`GTSC_GBTJ_Ks025` lowers it by 0.05693, while raising `Ks` from 0.25 to 0.30
lowers it by another 0.01125.  Full T added to `GSC_GBTJ_Ks025` changes the
rate by only +0.00184.

Therefore the revised long-time minimum model for the progressive-slowing full
response is **`GSC_GBTJ_Ks025`**.  It is much closer to full
`GTSC_GBTJ_Ks025` than `GTC_GB` in common-window rate, progressive exponent,
sink partition, stationary fraction, and event response.  T remains useful for
local response detail but is not required to reproduce the primary long-time
slowdown.  `GTC_GB` remains the correct no-shear, GB-only-sink control, not the
minimum surrogate for the full progressive-slowing branch.

The complete secant contrasts are in
`analysis/phase1_revised_factorial_comparisons.csv`.  These Phase-1 comparisons
supersede short fixed-step endpoint rankings but do not replace replicated
uncertainty estimates.

## Corrected and legacy roots

Strict legacy/project-history roots retained without pooling are:

- `results/campaigns/20260818T041229Z-49171551f9`
- `results/campaigns/20260818T121303Z-cf04b62e5c`
- `results/campaigns/20260818T220844Z-bed34ea2da`
- `results/video_runs/20260818T190827Z-409bd74b2d`
- strict swept-area continuations:
  `results/jerky_long_legacy_20260822/20260823T012535Z-793b9644d0`

Corrected `C_GB` roots are:

- diagnostic:
  `results/jerky_area_loss_diagnostic_20260822/20260823T002904Z-457cf9e522/C_GB-T900-s5101`
- long factorial:
  `results/jerky_area_loss_factorial_long_20260822/20260823T012141Z-b907865386/C_GB-T900-s5101`
- discriminating fast/slow:
  `results/jerky_area_loss_discriminating_20260822/20260823T102804Z-2ab688a139/C_GB_{FAST,SLOW}-T900-s5101`
- Phase-1 GTC GB-only control:
  `results/long_time_kinetics_900K_20260824/20260825T012949Z-6d83ee5c82/GTC_GB-T900-s5101`

Corrected true competing-sink `C_GBTJ` roots are:

- diagnostic:
  `results/jerky_area_loss_diagnostic_20260822/20260823T002904Z-457cf9e522/C_GBTJ-T900-s5101`
- long factorial:
  `results/jerky_area_loss_factorial_long_20260822/20260823T012141Z-b907865386/C_GBTJ-T900-s5101`
- discriminating fast/slow:
  `results/jerky_area_loss_discriminating_20260822/20260823T102804Z-2ab688a139/C_GBTJ_{FAST,SLOW}-T900-s5101`
- Phase-1 competing-sink cases: `GSC_GBTJ_Ks025`,
  `GTSC_GBTJ_Ks025`, and `GTSC_GBTJ_Ks030_LONG` under the canonical root.

These replace the earlier pre-area-loss dossier statement that no true
`C_GBTJ` model existed.  Historical reports remain unchanged for provenance.

## Movie paths and indexes

The Phase-1 machine-readable index is `analysis/phase1_movie_index.csv`; it
contains exact paths, source-frame counts, sizes, and SHA-256 hashes for 9/9
MP4 files:

- `B0-T900-s5101/microstructure.mp4`
- `GTSC_GBTJ_Ks025-T900-s5101/microstructure.mp4`
- `GTSC_GBTJ_Ks030_LONG-T900-s5101/microstructure.mp4`
- `GTC_GB-T900-s5101/microstructure.mp4`
- `GSC_GBTJ_Ks025-T900-s5101/microstructure.mp4`
- `G-T900-s5101/microstructure.mp4`
- `T-T900-s5101/microstructure.mp4`
- `GT-T900-s5101/microstructure.mp4`
- `QIU-T900-s5101/microstructure.mp4`

Each path is relative to the canonical campaign root.  Contact sheets are
`analysis/terminal_contact_sheet.png` and
`analysis/progress_matched_contact_sheet.png`.  Earlier authoritative indexes
are [jerky_area_loss_movie_index_20260823.md](jerky_area_loss_movie_index_20260823.md),
[jerky_mechanism_movie_index_20260821.md](jerky_mechanism_movie_index_20260821.md),
and [shear_scale_movie_index_20260823.csv](shear_scale_movie_index_20260823.csv).

## Analysis products

The common figure set contains `G`, `G2-G02`, `G3-G03`, `Gdot`, moving `n`,
`K2/K3`, best `Kn`, residual profiles, normalized distributions, distribution
collapse, G/T/C occupancy, sink partition, shear/backstress, jerkiness,
mechanism residence, activation work, and event-triggered response figures.
The principal tables are all under the campaign `analysis` directory:

- `grain_size_trajectory.csv`
- `growth_law_windows.csv` and `growth_law_profiles.csv`
- `normalized_distributions.csv` and `distribution_collapse.csv`
- `mechanism_residence_vs_grain_size.csv`
- `mechanism_joint_occupancy_vs_grain_size.csv`
- `frame_mechanisms_vs_grain_size.csv`
- `jerkiness_vs_grain_size.csv`
- `activation_work_vs_grain_size.csv`
- `event_triggered_velocity_vs_grain_size.csv`
- `event_triggered_response_summary.csv`
- `causal_nulls_vs_grain_size.csv`
- `grain_velocity_summary_vs_grain_size.csv`
- `grain_velocity_distribution_vs_grain_size.csv`
- `burst_wait_distribution_vs_grain_size.csv`
- the seven `phase1_*` audit/summary tables and
  `phase1_analysis_manifest.json`.

## Recommended next temperature/seed campaign

Do not start temperature variation yet.  First:

1. checkpoint-continue `GSC_GBTJ_Ks025` and `GTSC_GBTJ_Ks025` to 150,000
   steps, then reevaluate `n`, `K2/K3`, distribution collapse, residence, and
   recent grain loss;
2. retain `GTSC_GBTJ_Ks030_LONG` as the strong progressive-slowing diagnostic,
   extending it only if the purpose is to test persistence of the plateau;
3. diagnose/replay the QIU step-10,000 avalanche with finer frames and a
   morphology guard before using QIU in any inference set.

After those gates, the smallest informative temperature screen is 18 extreme-
temperature runs:

- regimes: `GTC_GB` (no-shear GB-sink control), `GSC_GBTJ_Ks025` (revised
  minimum), and `GTSC_GBTJ_Ks025` (full reference);
- temperatures: 800 and 1100 K;
- seeds: three paired seeds.

If ordering, sink partition, and progressive-versus-finite classification are
reproducible, add the two missing 900 K seeds for these three regimes (six new
runs) and reuse the present seed-5101 triplet, yielding a 27-condition
three-temperature inference set from 24 new runs.  Add B0 temperature controls
only if intrinsic mobility is made temperature-dependent.  Do not fit
`Q_app` until the same conservation, direct occupancy, stage, fine-trace null,
morphology, and movie audits pass at every temperature.
