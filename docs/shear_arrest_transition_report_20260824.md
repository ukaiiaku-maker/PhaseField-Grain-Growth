# Shear-backstress arrest-transition campaign

## Scope and provenance

This campaign resolves the 900 K shear-stiffness transition before any
temperature variation.  No 800 K or 1100 K production runs were launched.

| Role | Source SHA |
|---|---|
| Phase-A closeout and retained automated NO-GO provenance | `a7cd9f7b4f2f26da6251295cbed7e077fe9b01fc` |
| first force-balance instrumentation | `84dec7d561c55d25f30c4186f7c66b84485c1bd3` |
| corrected production separation of net backstress and event impulse | `c61f5f3e661cb45469aa0903a807bb94b821c4b2` |
| frozen two-seed crossover campaign specification | `10f0c58f96e5c6b13c93face7aa6685f57126609` |
| unbiased ordinary-frame/event-conditioned analysis separation | `718cebd8147d126501c129dbeae6c73bcc94d933` |
| experimental-observable, movie-index, and required-plot analysis | `71ef9e99c61aeecf55c677875dcb889d2c604b98` |

The retained legacy Phase-A source root is
`results/jerky_shear_scale_900K_20260823/20260823T153124Z-39f3ffde01`
(source `dbca0609e16e196f37bc3fd03346252fe9faac16`).  The corrected nine-run
localization root is
`results/jerky_shear_transition_900K_20260823/20260824T012827Z-75bfb9de8b`.
An earlier interrupted diagnostic root,
`results/jerky_shear_transition_900K_20260823/20260824T010800Z-9d5a8348b1`,
is superseded because its preliminary `p_net` definition included the
transient event impulse; it is excluded from every scientific comparison.

## Executed force pathway

The executed-code derivation is documented in
[shear_force_pathway_20260823.md](shear_force_pathway_20260823.md).  For the
production `local_memory` backend,

\[
\tau_{int}=-K_s s,\qquad p_{shear}=\beta\tau_{int},
\]

and the sharp-interface diagnostic balance is

\[
p_{net}=p_{cap}+p_{chem}+p_{shear},\qquad
p_{cap}=\gamma\kappa,\qquad p_{chem}=0.
\]

Capillarity is already supplied by the diffuse PF functional; only
`p_shear+p_event` is written as the additional antisymmetric PF driving field.
The executed total is therefore recorded separately as

\[
p_{applied,total}=p_{net}+p_{event}.
\]

The transient stochastic release impulse `p_event` is deliberately excluded
from the arrest ratio

\[
\chi_s=-\frac{p_{shear}}{p_{cap}+p_{chem}}.
\]

This continuous force is distinct from transition-state shear work,

\[
W_s^\ddagger=\tau_{resolved}V_\tau^\ddagger,
\]

which changes activation hazards and mode selection.  Compatibility filtering,
G/T/C residence, and GB/TJ sink competition remain separate coupled routes.
No source branch or hard gate is keyed directly on `K_s`.

## Frozen production parameters

All localization and seed-test runs use a 192 by 192 periodic grid,
`dx=1`, interface width 4, `dt=0.04`, GB energy 1, intrinsic mobility 4,
adaptive stepping, 267 generated grains equilibrated/compacted to 200, and a
120-grain target with a 10,000-step ceiling.  The event domain length is 12,
encounter density 0.20, core barrier 0.25 eV, `b` and `h` coefficients 0.02
eV, attempt frequency 10, packet size 2, packet window 1, shear coupling
`beta=0.35`, and shear trigger 0.25.

The climb/excess-volume values are fixed at excess volume per lost area 0.01,
point-defect formation volume 0.02, free-volume stiffness 0.05, climb trigger
quota 0.25, and release quota 1.0.  GB and TJ sink prefactors are 100000 with
nucleation/exchange/transport barriers 0.45/0.55/0.65 eV; TJ step length is 1,
TJ Burgers magnitude 0.25, compatibility tolerance 0.25, and TJ residual
stiffness 1 eV.  Only model family, seed, and `K_s` change in the robustness
screen.

## Verification before execution

The corrected production SHA passed 163 tests.  After separating ordinary
frame and event-conditioned analysis, the full suite passed 164 tests in
47.86 seconds.  With the experimental-observable and target-step analyses, the
final full suite passes 165 tests in 32.18 seconds.  The force decomposition,
no-double-capillarity pathway, event-impulse separation, `chi_s` sign
convention, diagnostic persistence, and analysis change-point/arrest helpers
all have direct regression coverage.

## Seed-5101 localization screen

| Family | `K_s` | end step | final grains | outcome | maximum absolute conservation residual |
|---|---:|---:|---:|---|---:|
| GTSC_GB | 0.25 | 7173 | 120 | target | 1.14e-11 |
| GTSC_GB | 0.30 | 10000 | 131 | censored | 1.32e-11 |
| GTSC_GB | 0.35 | 10000 | 150 | censored | 1.80e-11 |
| GSC_GBTJ | 0.25 | 6382 | 120 | target | 2.96e-11 |
| GSC_GBTJ | 0.30 | 10000 | 121 | censored | 5.80e-11 |
| GSC_GBTJ | 0.35 | 10000 | 149 | censored | 2.66e-11 |
| GTSC_GBTJ | 0.25 | 5417 | 120 | target | 1.59e-11 |
| GTSC_GBTJ | 0.30 | 10000 | 132 | censored | 5.34e-11 |
| GTSC_GBTJ | 0.35 | 10000 | 152 | censored | 3.66e-11 |

Thus all three families have a practical terminal first-passage cliff between
0.25 and 0.30 for seed 5101.  This is not a failure-status artifact: every
ceiling case has a terminal completed manifest, finite continued releases,
bounded morphology, and conservation better than `6e-11`.

The topology-window growth rates are:

| Family | `K_s` | 190 to 160 | 160 to 140 | 140 to 120 |
|---|---:|---:|---:|---:|
| GTSC_GB | 0.10 | 0.017044 | 0.015087 | 0.028765 |
| GTSC_GB | 0.20 | 0.015465 | 0.010196 | 0.018232 |
| GTSC_GB | 0.25 | 0.014218 | 0.005520 | 0.004651 |
| GTSC_GB | 0.30 | 0.010122 | 0.002624 | unavailable |
| GTSC_GB | 0.35 | 0.007499 | unavailable | unavailable |
| GTSC_GB | 0.40 | 0.004677 | unavailable | unavailable |
| GSC_GBTJ | 0.10 | 0.022478 | 0.024033 | 0.023025 |
| GSC_GBTJ | 0.20 | 0.018095 | 0.010603 | 0.010887 |
| GSC_GBTJ | 0.25 | 0.016420 | 0.003245 | 0.005151 |
| GSC_GBTJ | 0.30 | 0.011676 | 0.002774 | unavailable |
| GSC_GBTJ | 0.35 | 0.009529 | unavailable | unavailable |
| GSC_GBTJ | 0.40 | 0.004604 | unavailable | unavailable |
| GTSC_GBTJ | 0.10 | 0.024056 | 0.022763 | 0.028228 |
| GTSC_GBTJ | 0.20 | 0.018473 | 0.007778 | 0.008201 |
| GTSC_GBTJ | 0.25 | 0.016941 | 0.005382 | 0.005055 |
| GTSC_GBTJ | 0.30 | 0.010839 | 0.002258 | unavailable |
| GTSC_GBTJ | 0.35 | 0.009574 | unavailable | unavailable |
| GTSC_GBTJ | 0.40 | 0.009565 | unavailable | unavailable |

Segmented kinetic fits put the seed-5101 middle-window crossover at 0.20 for
GTSC_GB and 0.25 for both GBTJ families.  Those point estimates are subordinate
to the directly observed 0.25/0.30 terminal bracket and require the cross-seed
test below.

## Force-balance and morphology interpretation

The corrected ordinary-frame results show a smooth rise in backstress.  In the
190-to-160 window, length-weighted `p95(chi_s)` rises from approximately
0.50-0.51 at `K_s=0.25` to 0.68-0.76 at 0.35.  Only about 1.1-1.7% of active
GB length lies in `0.8 < chi_s < 1.2`; the fraction above one is approximately
2.4-3.8%.  Mechanical cancellation is therefore real but localized, not a
spatially uniform frozen interface.  The large kinetic effect must be
interpreted as the coupled result of continuous velocity reduction,
transition-state work/mode accessibility, compatibility residence, and finite
stochastic release.

All nine terminal structures remain compact and polygonal.  Mean pixel
compactness lies between 1.766 and 1.792, and matched-state high-frequency
boundary-power fractions lie near 0.288-0.297.  The roughness outputs and
terminal sheet are under the corrected root's `analysis_fast_ordinary`
directory.  All nine six-panel GIFs are present as
`<run>/microstructure.gif` under the corrected localization root.

## Cross-seed decision gate

The seed-5101 bracket is `K_s,below=0.25` and `K_s,above=0.30`.  GSC_GBTJ is
included because its 0.30 result (121 grains) differs enough from the two full
models (131-132 grains) to affect the minimum-model conclusion.  The frozen
campaign therefore contains 12 runs: three families by two stiffnesses by
seeds 5102 and 5103.  Its root is
`results/jerky_shear_seed_test_900K_20260824/20260824T101946Z-075b81ced4`.

The completed first-passage results are:

| Family | `K_s` | seed 5101 | seed 5102 | seed 5103 | target fraction |
|---|---:|---|---|---|---:|
| GTSC_GB | 0.25 | 7173 / 120 | 7079 / 120 | 5503 / 120 | 3/3 |
| GTSC_GB | 0.30 | 10000 / 131 | 10000 / 123 | 10000 / 127 | 0/3 |
| GSC_GBTJ | 0.25 | 6382 / 120 | 5396 / 120 | 3791 / 120 | 3/3 |
| GSC_GBTJ | 0.30 | 10000 / 121 | 10000 / 128 | 10000 / 128 | 0/3 |
| GTSC_GBTJ | 0.25 | 5417 / 120 | 5170 / 120 | 4803 / 120 | 3/3 |
| GTSC_GBTJ | 0.30 | 10000 / 132 | 10000 / 121 | 10000 / 131 | 0/3 |

Each entry is `end step / final grain count`.  Every 0.25 realization reaches
the target and every 0.30 realization is terminally censored.  The residual
terminal grain count is seed sensitive, but the classification is not.

## Corrected crossover and stage-resolved kinetics

The robust practical bracket is

\[
0.25 < K_s^{crit} \leq 0.30.
\]

A convenient descriptive midpoint is `0.275 +/- 0.025`; the half-width is the
tested bracket resolution, not a statistical confidence interval.  The
segmented fits select 0.25 as the critical tested stiffness in all three
families.  A force-only fitted critical point is deliberately not reported:
only three corrected ordinary-force stiffness levels are available, so that
fit is underdetermined.  The binary first-passage bracket is the defensible
estimate.

Cross-seed mean rates in matched topology windows are:

| Family | `K_s` | 190 to 160 | 160 to 140 | 140 to 120 |
|---|---:|---:|---:|---:|
| GSC_GBTJ | 0.25 | 0.017657 | 0.009173 | 0.006677 |
| GSC_GBTJ | 0.30 | 0.013907 | 0.005174 | 0.001633 partial |
| GTSC_GB | 0.25 | 0.016445 | 0.008779 | 0.004249 |
| GTSC_GB | 0.30 | 0.014280 | 0.004228 | 0.002149 partial |
| GTSC_GBTJ | 0.25 | 0.018269 | 0.008678 | 0.005599 |
| GTSC_GBTJ | 0.30 | 0.013119 | 0.003652 | 0.001435 partial |

All three seeds contribute to the first two windows.  The 0.30 late-window
values are partial summaries because no 0.30 run reaches 120 grains.  In the
middle window, the 0.30/0.25 rate ratio is approximately 0.56 for GSC_GBTJ,
0.48 for GTSC_GB, and 0.42 for GTSC_GBTJ.  Thus the terminal cliff is preceded
by a smooth, family-dependent kinetic slowdown.

## Ordinary-frame force balance and normalized energy

Cross-seed means in the 160-to-140 window are:

| Family | `K_s` | p95 `chi_s` | length `chi_s>0.8` | length `0.8<chi_s<1.2` | length `chi_s>1` | mean `abs(tau_int)` | energy / active GB length |
|---|---:|---:|---:|---:|---:|---:|---:|
| GSC_GBTJ | 0.25 | 0.6178 | 0.0371 | 0.0159 | 0.0266 | 0.03591 | 0.000655 |
| GSC_GBTJ | 0.30 | 0.6578 | 0.0393 | 0.0131 | 0.0313 | 0.04067 | 0.000660 |
| GTSC_GB | 0.25 | 0.6105 | 0.0358 | 0.0138 | 0.0277 | 0.03577 | 0.000676 |
| GTSC_GB | 0.30 | 0.7655 | 0.0472 | 0.0150 | 0.0380 | 0.04074 | 0.000658 |
| GTSC_GBTJ | 0.25 | 0.6142 | 0.0372 | 0.0120 | 0.0302 | 0.03605 | 0.000671 |
| GTSC_GBTJ | 0.30 | 0.6870 | 0.0418 | 0.0126 | 0.0352 | 0.04128 | 0.000677 |

Stress and the high percentile of `chi_s` rise smoothly, while near-unity
cancellation occupies only 1.2-1.6% of active length and local drive reversal
occupies 2.7-3.8%.  The stored shear energy per active length is nearly flat
across the bracket, exactly as expected when the strain required for balance
falls with stiffness.  Per-domain, length-weighted, active-domain, and local
p50/p90/p95/p99 forms are retained in the analysis tables and plots.  The
ordinary-frame metrics exclude the stochastic impulse; its trajectory remains
available separately in the event-trace outputs.

## Activation work, compatibility residence, and arrest

The p95 transition-work magnitude `abs(tau V_tau^dagger)/kBT` rises from
0.5003 to 0.5429 for GTSC_GB, 0.6073 to 0.7107 for GSC_GBTJ, and 0.5916 to
0.6948 for GTSC_GBTJ between 0.25 and 0.30.  The fraction of events receiving
more than a factor-of-two shear hazard change rises only from 0.0259 to 0.0288,
0.0367 to 0.0528, and 0.0346 to 0.0516, respectively.  Transition-state work
therefore contributes, but its small affected fraction cannot alone explain
the macroscopic slowdown.

The median detected arrest episode is two solver steps in every group.  Across
the three seeds, the longest-episode time ranges are:

| Family | `K_s=0.25` | `K_s=0.30` |
|---|---:|---:|
| GSC_GBTJ | 3.20 to 24.64 | 6.00 to 8.76 |
| GTSC_GB | 5.48 to 14.20 | 8.16 to 14.36 |
| GTSC_GBTJ | 5.92 to 22.24 | 10.20 to 104.56 |

Only the full GBTJ family develops the extreme 0.30 residence tail (the two
new seeds reach 104.56 and 63.20).  This is evidence for an interaction between
continuous backstress and TJ compatibility residence, not for a universal
solver freeze.

## Sink partitioning and direct occupancy

Cross-seed middle-window means are:

| Family | `K_s` | GB sink fraction | TJ sink fraction | G occupancy | T occupancy | C occupancy |
|---|---:|---:|---:|---:|---:|---:|---:|
| GTSC_GB | 0.25 | 1.000 | 0.000 | 0.0187 | 0.0144 | 0.00884 |
| GTSC_GB | 0.30 | 1.000 | 0.000 | 0.0109 | 0.0124 | 0.00745 |
| GSC_GBTJ | 0.25 | 0.388 | 0.612 | 0.0236 | 0.00003 | 0.00931 |
| GSC_GBTJ | 0.30 | 0.362 | 0.638 | 0.0114 | 0.00002 | 0.00728 |
| GTSC_GBTJ | 0.25 | 0.382 | 0.618 | 0.0230 | 0.00003 | 0.00933 |
| GTSC_GBTJ | 0.30 | 0.365 | 0.635 | 0.0110 | 0.00005 | 0.00703 |

There is no wholesale sink switch across the transition.  G and C occupancy
fall, TJ occupancy remains tiny, and the GBTJ families retain a mild TJ-sink
majority.  The sharp first-passage change therefore is not explained by an
abrupt sink-partition discontinuity.

## Event-triggered response and causal-null tests

At an eight-step response window, actual event-linked changes in `abs(v_n)`
range from approximately 0.104 to 0.136 across family/stiffness groups.
Enrichment over circular-time, block-time, entity, and topology-matched nulls
ranges from approximately 0.100 to 0.134 and is positive for every seed.  In
the event trajectories, mean post-event steps 1 through 5 exceed mean pre-event
steps -5 through -1 by 0.209 to 0.269 in every family/stiffness/seed group.

Finite local release from arrested states therefore survives all causal-null
comparisons.  Its amplitude is similar at 0.25 and 0.30; the crossover is driven
by reduced continuous velocity and altered residence/accessibility/frequency,
not by loss of the event response itself.  The exact event-step value can be
lower than the subsequent release, which is why event traces are not mixed
into ordinary-frame force statistics.

## Experimental-observable consequences

Cross-seed middle-window means, with the late-window stationary/wait measures
shown in the final two columns, are:

| Family | `K_s` | growth rate | exponential-fit rate | top-10% burst share | p95 grain `abs(v_r)` | stationary fraction | maximum wait |
|---|---:|---:|---:|---:|---:|---:|---:|
| GTSC_GB | 0.25 | 0.00878 | 0.0737 | 0.272 | 0.0590 | 0.0000 | 0.00 |
| GTSC_GB | 0.30 | 0.00423 | 0.0499 | 0.239 | 0.0348 | 0.0031 | 0.67 |
| GSC_GBTJ | 0.25 | 0.00917 | 0.0624 | 0.254 | 0.0637 | 0.0000 | 0.00 |
| GSC_GBTJ | 0.30 | 0.00517 | 0.0524 | 0.248 | 0.0431 | 0.0175 | 2.00 |
| GTSC_GBTJ | 0.25 | 0.00868 | 0.0631 | 0.290 | 0.0649 | 0.0000 | 0.00 |
| GTSC_GBTJ | 0.30 | 0.00365 | 0.0401 | 0.251 | 0.0412 | 0.0302 | 2.00 |

The retained table also contains burst/wait distributions, event-linked local
motion, grain-level velocity summaries, and effective topology-window
kinetics.  Ordinary output cadence under-resolves short waits, so the
solver-step arrest analysis is the preferred duration measure.  No external
quantitative experimental calibration metric is present in the repository;
consequently a stiffness range must be preserved rather than declaring a
unique calibrated value.

## Morphology, movies, and minimum model

All 12 seed-test runs have completed terminal manifests, event traces,
conservation audits, roughness records, and six-panel GIFs at each run's
`microstructure.gif`.  The combined movie index contains 30 entries covering
the retained legacy, corrected localization, and seed-test campaigns; every
indexed movie exists.  The seed-test terminal contact sheet is
`analysis/seed_contact_sheet.png`.  Cross-seed terminal compactness is about
1.773-1.811 and high-frequency boundary power about 0.287-0.290.  There is no
morphological discontinuity or grid-scale roughening at the crossover.

GSC_GBTJ reproduces the terminal crossover, topology-window slowdown, sink and
occupancy trends, and causal event response.  It is therefore the revised
minimum model for screening the existence and location of the transition.  It
does not reproduce the extreme arrest-duration tail of full GTSC_GBTJ; TJ
compatibility is required to study that residence/release-tail mechanism.  The
temperature comparison pairs remain GTC_GB versus GTSC_GB and GSC_GBTJ versus
GTSC_GBTJ.

## Scientific decision

The strong-shear branch is an emergent backstress-limited mechanical-arrest
regime, not a numerical artifact.  The evidence is mutually consistent:
smooth force growth, localized near-balance and drive reversal, exact
conservation, bounded morphology and roughness, finite causal release events,
no `K_s`-conditioned source gate, and a perfectly replicated 0.25/0.30
first-passage classification across three seeds and all three families.  The
mechanism is coupled: continuous backstress supplies the primary velocity
reduction, while activation work, mode accessibility, and compatibility
residence control access to and escape from locally arrested states.

## Recommendation for the temperature/seed campaign

Do not launch the previously contemplated 24-run temperature campaign as a
single-stiffness production sweep.  Preserve both sides of the robust bracket:

- use `K_s=0.25` as the mobile, near-crossover condition;
- use `K_s=0.30` as the mechanically arrested, intentionally censored condition;
- begin with one seed at 800 K and 1100 K for the two retained reduced/full
  comparisons, reusing the no-S control where diagnostics are compatible;
- require interpretable family ordering, conservation, morphology, and event
  response at that pilot gate before adding seeds 5102 and 5103.

If only one stiffness can be afforded for a conventional normal-growth study,
choose 0.25.  A study of arrest physics should carry both 0.25 and 0.30 and
treat censoring at 0.30 as a scientific outcome rather than a failed run.  No
temperature simulations were launched during this campaign.

## Artifact index

The definitive combined analysis is under
`results/jerky_shear_seed_test_900K_20260824/20260824T101946Z-075b81ced4/analysis`.
It contains the required target-step/time, kinetics, arrest, force-balance,
normalized-energy, transition-work, sink, occupancy, event-response,
causal-null, morphology, contact-sheet, and experimental-observable outputs.
The canonical combined movie index is `analysis/movie_index.csv` under that
root.  The corrected localization roughness audit remains under
`results/jerky_shear_transition_900K_20260823/20260824T012827Z-75bfb9de8b/analysis_fast_ordinary`.
