# Weekend development report: stochastic disconnection-mediated grain growth

> Live campaign report. Sections marked **pending production inference** are
> intentionally not interpreted until their immutable campaign manifest is
> complete. The final archival commit and `v1.0-weekend-campaign` tag will
> replace this notice.

## 1. Reproducibility snapshot

- Development branch: `codex/weekend-disconnection-pf-2026-08-14`.
- Validated curvature-scaling milestone: tag `v0.8-scaling-validated`.
- Original mechanism simulation launch SHA:
  `9e8796e6e9e13bc4bb418662cfc33601adeaa502`.
- Exact mechanism composite: 165 matched completed conditions in
  `results/campaigns/20260817T042553Z-composite-28a9a5a185`, assembled with
  whole-condition repairs for the rejected/interrupted SC4, P1, J1, J2, and
  J3 entries and immutable source-run provenance.
- Independent 256-square convergence campaign:
  `results/campaigns/20260815T083416Z-49423c71d0` (completed).
- Current fast suite: 103 tests passed in 89.98 seconds under the live
  production load.
- Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3.

Every run manifest contains its launch SHA, canonical configuration hash,
software versions, seed, initial-condition provenance, stopping criteria, and
restart-artifact checksums. Dense checkpoints remain outside Git. Compact
summaries, validation results, plots, failure records, and manifests are
versioned. Manifests, checkpoint archives, and their companion JSON metadata
are replaced atomically. The archive embeds the authoritative metadata generation, so an
interruption between replacements cannot pair new arrays with stale state on
resume; legacy split checkpoints remain readable.

## 2. Architecture implemented

The implementation has four deliberately separated layers:

1. A constrained, pairwise Qiu/Steinbach-type multiphase-field kernel with
   permanent phase extinction, adaptive stepping, periodic/Neumann boundaries,
   and compact active support.
2. Persistent grain, finite GB-domain, and TJ entities, with geometry-driven
   encounter hazards independent of activation residence time.
3. Discrete `(h,b,Nv)` modes, cumulative first-passage clocks, multihit memory,
   local and nonlocal shear mechanics, free-volume storage, serial climb, and
   entity-attached obstacles.
4. Immutable campaign execution, exact restart, scaling/Arrhenius inference,
   intermittency diagnostics, publication plots, and final-summary aggregation.

Campaign workers receive one run per multiprocessing chunk. This prevents a
long coupled-mechanism case from trapping later runs behind it while other
workers sit idle near the end of a heterogeneous matrix.

The complete module map is in `docs/model_hierarchy.md`; the equations actually
executed are in `docs/physics_equations.md`.

## 3. Mathematical model

The baseline evolves constrained order parameters with the pairwise Qiu
functional and recovers the sharp-interface target

\[
v_n=M\Gamma\kappa.
\]

Mode-specific hazards use

\[
r_m=N_{s,m}\nu_m\exp[-\max(0,\Delta G_{0,m}-W_m)/(k_BT)],
\]

with an attempt-limited zero-barrier state. The work contains independently
signed normal, resolved-shear, and vacancy-chemical-potential contributions.
Parallel admissible modes compete by summed hazard; necessary combinations and
climb stages add residence times.

All stochastic events are sampled by cumulative hazard. Persistent `K`-hit
completion is Erlang; packet-reset completion requires `K` hits within the
finite, checkpointed `packet_window_time`. A separate geometric clock creates
difficult configurations from GB measure change or TJ travel.

Local shear stores \(E_s=K_s s^2/2\) and feeds back with the negative energy
gradient. The nonlocal backend applies periodic FFT incompatibility projection
to event eigenstrain. Free volume stores \(E_q=K_q q^2/2\), drives
Butler--Volmer exchange, and uses \(\tau=C\ell^2/D\) for transport. Full climb
is the serial chain nucleation, exchange, transport, quota completion.

## 4. Source and provenance map

`docs/source_manifest.md` records hashes and uses for all supplied archives,
documents, MATLAB/PF prototypes, and the Qiu code. `docs/qiu_code_audit.md`
documents the Qiu PF, reference construction, internal-stress kernel, elastic
driving, numerics, and limitations. `docs/external_sources.md` records licenses
and primary literature. Qiu source remains pristine under ignored external
storage; production code is an independently tested implementation.

## 5. Validation results

### Numerical PF

`results/validation/numerical_validation.json` passes. Across six shrinking
circle cases, `R^2(t)` has `R² >= 0.99973`; doubling mobility changes the slope
by factors 1.99944--2.00069. Sharp-interface slope error decreases from about
11.25% at `dx=1` to 2.48% at `dx=0.5`. Surface-energy grid anisotropy is 0.847%,
and equal-energy TJ angles differ from 120 degrees by at most 4.29 degrees.

The independent 256-square case evolved 300 equilibrated grains to 90 at step
3270. Its fixed-`n=2` coefficient is 1.17534 versus 1.13852 for the 192-square
20-run ensemble, a 3.23% difference. Its one-realization free fit gives
`n=2.442`, but gains only 0.00093 in `R²` and has residual autocorrelation 0.997;
that exponent is retained as sensitivity rather than interpreted as a
one-sample confidence result. Restart hashes and final plots are audited.

### Stochastic engine

`results/validation/stochastic_validation.json` passes. For rate 2.5, the
single-hit mean is 0.40366 versus 0.4 and the exponential KS p-value is 0.824.
Paired event times remain invariant across an eightfold timestep change to
`3.8e-14`. Persistent `K=5` gives mean/CV 1.9980/0.44658 versus
2.0/0.44721. Packet completion is 0.37649 versus the exact Poisson tail
0.37729. Parallel and serial means match 1/7 and 0.7.

A four-temperature fixed-geometry isolation test recovers an imposed 0.65 eV
barrier as 0.6516 +/- 0.0004 eV for single hit and
0.6490 +/- 0.0013 eV for persistent `K=5`.

### Geometry, modes, mechanics, and climb

The fast suite covers persistent state retirement, TJ motion, geometric measure
invariance, exact restart, isotropic zero-driving sampling, shear and curvature
selection, combined competing hazards, finite minimum Burgers magnitude,
angular quadrature convergence, feasible TJ combinations, Burgers closure,
free-volume balance, event-strain sums, attempt limiting, Butler--Volmer/Onsager
limits, diffusivity and length-squared transport, and distinct serial-stage
rates. Production integration tests cover finite PF event release, quenched
barriers, mixed shear/climb release, renewal-window statistics and restart, and
output-cadence invariance.

### Qiu-type benchmarks

`results/validation/qiu_regression_benchmarks_213728b2.json` passes on a clean
post-sign-correction revision. In the matched
24-grain polycrystal, full-field shear reduces velocity-curvature correlation
from 0.4358 to 0.3897 and raises active reverse-curvature motion from 12.70% to
15.35%. The four-grain benchmark raises active reverse motion from 3.09% to
4.18%. Eigenstrain, nonlocal stress, elastic feedback, and phase-field
divergence are all finite.

### Baseline ensemble scaling

The exact-checkpoint 20-realization baseline extension reaches 60--61 grains.
Over the objective population window, radius-space regression gives
`n=1.9150`, bootstrap 95% CI `[1.7079,2.1503]`, `K=0.8931`, and
`R²=0.99953`. Fixing `n=2` changes `R²` by only `-3.21e-5`; independent mean-
radius and perimeter-radius fits give `n=2.056` and `n=2.009`. The earlier
`n≈1` short-window result is retained as a failed/transient inference, not
deleted.

## 6. Production mechanism matrix

The finalized matrix contains B0--B1, G1--G3, T1--T3, S1--S3, C1--C5,
SC1--SC4, P1--P5, Q0--Q1, E0--E2, and J1--J3, each with five matched
realizations at 900 K. The original campaign contributes 146 valid completed
conditions; corrected whole-condition campaigns supply the remaining 19.
The exact 165-condition composite is complete, audited, analyzed, and plotted.
It contains 33 regimes with five matched realizations at 900 K and zero
tracebacks. Fifteen superseded Q1/E1/SC3 duplicates are replaced by the later
corrected runs; corrected SC4/J1/J2 and P1/J3 conditions supply the previously
missing inference set. The authoritative table is
`results/production_summaries/mechanism_composite_165_summary.csv` and 256
PNG/PDF plot files are under `results/plots/mechanism_composite_165/`.

Four superseded partial matrices remain preserved. The first diagnosed tracker
throughput; the second stopped before climb runs after detecting unwired
Butler--Volmer/transport production paths and an incomplete packet-renewal
interpretation. The third exposed same-step first passages sampled after a
single-hit release had already changed the physical blocked state. Its five-seed
G1 ledger contained 98,540 primitive hits for 79,289 releases. The corrected
clock stops at the physical completion time; unit, packet-window, and actual
blocked-GB ledger regressions enforce that behavior while continuous easy-mode
flux retains multiple passages. The fourth exposed the corresponding encounter-
side boundary: multiple geometric passages were consumed after the first had
already created a blocked state. The corrected geometric clock stops at the
exact first encounter coordinate for gated GB/TJ callers, while its ungated
multi-passage behavior remains available. Audit records are under
`results/validation/`; none of the partial matrices is included in final
inference. A later audit found two additional production boundaries: explicit
GB compatibility was missing from the encounter predicate and then remained
reachable through an ungated explicit-mode fallback, while per-event CSV
flushing made valid high-rate J ledgers grow to multiple gigabytes.
Geometry-gate regressions now cover both paths. Event ledgers support
checkpoint-member gzip and fixed-schema Zstandard Parquet; the latter commits
atomic parts and truncates by the authoritative checkpoint part count. The
first Parquet J1/J2 rerun then exposed a separate physical conservation error:
both endpoints received the same Burgers increment, TJ relaxation omitted the
configured packet factor, strict unresolved TJs continued accepting GB events,
and the finite-residual model lacked its declared positive-definite back stress.
Those zero-traceback but residual-runaway campaigns are preserved and excluded.
The repaired model applies opposite endpoint signs, packet-consistent relaxing
modes, nearest-domain TJ adjacency, strict-state gating, and event-by-event
residual-energy rate feedback. The same audit found positive eigenstrain
self-work in the FFT backend, opposite to both elastic-energy stability and the
explicit minus sign in Qiu's PF elastic driving. The corrected self-stress is
the negative energy derivative. Q1, E1, SC3, SC4, and J2 are therefore rerun
whole alongside J1. A corrected 64x64 production-path smoke generated 266 J1
and 1,002 J2 ledger rows in 20 steps, with J2 residual norm bounded by 0.707,
rather than the prior attempt-limited explosion. The now-superseded SC4 repair completed
five seeds with zero tracebacks and gives `n=1.3248` (95% CI
`[1.0000, 2.2126]`), `K=0.1763`, jerkiness CV `1.162`, and Fano factor
`280.95`; it is retained as diagnostic evidence but excluded from inference
because it used the autocatalytic full-field sign. The corrected full-field/TJ
replacements are included in the exact composite; superseded runs remain excluded.

The corrected-sign SC4 replacement is now complete for five seeds. Its interim
fit gives `n=1.000` at the scan bound, `K=0.15522`, fit `R²=0.96374`, CV
`1.17168`, Fano `584.80`, reverse-motion fraction `0.65717`, and velocity--
curvature `R²=6.3e-5` across 1,254,987 primitive events. Within its serial climb
stages, expected `1/r` resistance is 99.15% transport, 0.74% exchange, and
0.11% nucleation. This corrected result is one of the five-seed conditions
admitted to the exact composite.

The corrected explicit-compatibility J1 replacement is complete for five
seeds at 900 K. It is physically stagnant under the objective observability
gate: characteristic radius changes by only 1.023% (minimum required 2%) while
the mean grain population falls from 200.0 to 186.2, so the analysis reports
`n=NaN` and `K=0` rather than a bound-dependent exponent. The retained event
statistics show strong intermittency (CV `3.23145`, Fano `291.36`, burstiness
`0.44233`), 71.70% stationary intervals, 60.45% reverse motion, and velocity--
curvature `R²=6.73e-6` across 3,606,707 primitive events. There are 426,780
TJ endpoint-failure rows, an incidence of 13.43% per completed GB mode event;
60.55% of those failures involve nominally easy modes. TJ activation accounts
for 97.07% of the event-conditioned expected residence. This demonstrates a
strong intrinsic compatibility arrest, but not the requested coexistence of
strong jerkiness with sustained smooth scaling in this parameterization.

The corrected strict persistent-multihit J2 replacement is also complete for
five seeds. It is fully stagnant: characteristic radius does not change and
the mean grain population moves only from 200.0 to 199.4 over physical time
14--140, so `n=NaN` and `K=0`. Nevertheless, it records 68,815,481 primitive
events, CV `4.24167`, Fano `3208.40`, burstiness `0.11194`, 58.18% reverse
motion, and velocity--curvature `R²=2.21e-5`. Its 66,545,846 TJ endpoint-
failure rows have incidence 1.094 per completed GB event (two endpoints can
fail); 70.68% involve nominally easy modes. Relative to J1, persistent
multihit kinetics raises CV and Fano but lowers waiting-time burstiness, so the
intermittency effect is metric dependent. All 30 corrected repair runs are now
complete with zero tracebacks and audited restart checksums.

The exact composite resolves the matched mechanism contrasts. With geometric
GB encounters held fixed, G1→G2→G3 changes the event history from single hit
to persistent and packet-reset multihit: all three remain at `n=1`, while `K`
falls 0.17450→0.10734→0.07307 and CV rises 1.330→1.516→1.686. Across the
distinct geometric scaling controls, however, the GB-area analogue P1 has
`n=1`, TJ pinning P4 has `n=1.3935` (95% CI `[1.0000,1.7691]`), and swept-area
spatial controls P2/P3 reach the upper scan bound `n=6`. The latter are lower-
bound statements on steepness, not uniquely identified sixth-order laws.

Climb/free-volume accommodation alone is strongly decorrelating but only
moderately intermittent: C1--C5 have CV 0.994--1.178, reverse-motion fraction
0.704--0.802, and velocity--curvature `R²=1.6e-4`--`8.2e-4`. In serial C5,
the event-conditioned residence fraction over all observed primitive types is
88.88% transport, 0.84% exchange, 0.10% nucleation, and 10.19% generic
activation; within the three climb stages it is 98.96%, 0.93%, and 0.11%.
The coupled SC results are architecture-dependent. Independently necessary
SC1 has `K=0.07496`, close to slow serial-climb C5 (`0.07226`), and 99.28% of
its climb-stage residence is transport. The mixed co-relaxing event SC2 is
much faster (`K=0.62315`), parallel competition SC3 is intermediate
(`K=0.34784`), and strict SC4 returns to `K=0.15522` with 98.75% of all
event-conditioned residence in transport and 0.39% in TJ activation. Thus
shear and climb are neither universally additive nor governed by one fixed
fraction: serial constraints are climb-dominated, whereas mixed/parallel
coupling removes or bypasses that residence bottleneck.

The discrete isotropic spectrum generates signed event shear even without
feedback (E0 accumulated shear 158.50). Full-field selection E1 raises that
to 4256.82 and lowers `K` from 0.59314 to 0.37073; local-memory E2 gives
171.27 and `K=0.53170`. This establishes stress/curvature mode selection and
strong backend sensitivity, but not quantitative equality to a continuously
coupled Qiu field. Explicit TJ conservation is more decisive: J1/J2 are
stagnant, whereas geometric surrogate J3 has `n=1`, `K=0.13380`, CV 1.449,
and Fano 107.25. The surrogate therefore does not reproduce explicit-mode
arrest or its extreme intermittency. In SC1--SC4, absolute signed accumulated
strain is 99.920--99.997% shear component and 0.003--0.080% volumetric
component; these are component fractions of signed strain, not energy or
causal resistance fractions.

Final analysis reads checkpointed event ledgers in bounded-memory batches. A
production regression over 31,255,548 primitive rows produced a byte-identical
summary and identical diagnostics apart from sub-`8e-12` summation-order roundoff,
while reducing observed RSS by at least 68.2% (5.66 GB to 1.84 GB).
Event-overlay, waiting-time, type-count, and TJ-failure plots use the same
streaming path. A visually inspected 3.23-million-event E1 plot pass used only
358 MB maximum RSS.

TJ residual Burgers vectors and the other persistent TJ entity fields are now
part of the authoritative restart archive. The live J1/J2 repair runs were
verified to have started fresh after the earlier campaign restart, rather than
resuming from a checkpoint lacking those fields; they remain valid while
uninterrupted. The scope audit is in
`results/validation/tj_checkpoint_state_20260816.json`.

## 7. Growth-law results

The exact composite reports free and fixed-`n=2`
fits, objective topology windows, realization bootstrap intervals, exponent
profiles, local exponents, residual autocorrelation, radius-measure sensitivity,
population-band sensitivity, and source-manuscript Class-B/Class-C comparators.
An ensemble whose fitted-window characteristic radius changes by less than 2%
is classified as stagnant/censored before nonlinear fitting: its summary reports
`n=NaN`, `K=0`, and the explicit observability reason. Event, boundary,
intermittency, and TJ-failure diagnostics are still retained. This prevents an
optimizer from assigning a bound-dependent exponent and tiny positive growth
coefficient to an exactly flat trajectory.

At 900 K, B0/B1 give `n=1.6308/1.6345`; P4 gives `n=1.3935`; 23 activated
regimes reach the lower scan bound `n=1`; P2/P3 reach the upper bound `n=6`;
and J1/J2 are censored as physically stagnant. The mechanism coefficients span
more than five orders of magnitude. Among activated families, event-rich
E1/SC3 have Fano factors 16,588/15,700 and velocity--curvature `R²` below
`2.3e-5`, while SC4 has CV 1.172 and Fano 584.8. This breadth is a result, not
a claim that every bound-limited exponent is uniquely identified.

## 8. Temperature and activation-energy campaigns

Two distinct four-temperature experiments are configured:

- mechanism isolation: 11 selected barrier regimes at 800, 900, 1000, 1100 K
  with intrinsic PF mobility held constant, plus a B1 intrinsic Arrhenius control;
- fully physical: 11 `FP-*` regimes at the same temperatures, including a
  0.45 eV intrinsic-mobility barrier normalized to `M(900 K)≈4`.

Each condition uses the same five initial structures. A common exponent is fit
across temperature before extracting `K_n(T)`. Event activation is estimated
separately from primitive event counts divided by integrated GB-domain-time
exposure, including zero-event exposure. Adjacent-temperature slopes and
Arrhenius curvature are retained. A grain-growth activation fit is suppressed
when any series changes radius by less than 2%, rather than assigning a
spurious activation energy to stagnation. The completed 20-run B1 control uses
five matched seeds at each temperature and recovers its imposed 0.45 eV
intrinsic barrier as `0.450000 eV` with bootstrap 95% CI
`[0.437646, 0.462170]`, Arrhenius `R²=1.0`, and all three adjacent-temperature
slopes equal to 0.45 eV within `4e-11`. Its common exponent is `1.63449`
(95% CI `[1.43851,1.82699]`) over a 70.2% radius increase. The remaining
selected and fully physical families remain **pending production inference**.

The selected campaign has completed 204/240 runs: B1, G2, T2, S2, C2--C5,
SC3, and P2 each have all 20 matched temperature/seed conditions. The remaining
P4 conditions continue with zero tracebacks. The fully
physical campaign is running at `results/campaigns/20260817T082103Z-01b68e843b` with
220 enumerated conditions and 10 workers; all initially verified manifests
record clean SHA `878081a206a4b75f243ddfd25821b252c954a87c`. An earlier 10-run
launch that recorded a dirty SHA is preserved but excluded from inference.
The completed G2 persistent-multihit family provides the first production
event/growth separation. Its event rate per integrated GB-domain time recovers
the imposed 0.25 eV barrier as `0.250122 eV` (95% CI
`[0.240313,0.260878]`, `R²=0.99920`), while its coarse-grained growth
coefficient has only `Q_app=0.143697 eV` (95% CI
`[0.127936,0.159105]`, `R²=0.99981`). The common exponent is at the lower
scan bound `n=1`; trajectory CV falls from 1.705 to 1.361 over 800--1100 K
while event-count Fano rises from 43.3 to 230.8.

The completed T2 persistent-TJ family also gives common `n=1.000` at the
lower scan bound, with coefficients increasing from `0.041912` at 800 K to
`0.126131` at 1100 K. Grain-growth activation is `0.275837 eV` (95% CI
`[0.251763,0.304763]`, Arrhenius `R²=0.99605`); event-rate activation is
`0.266257 eV` (95% CI `[0.257348,0.274163]`, `R²=0.99926`). The local growth
slopes range from 0.238 to 0.317 eV, whereas local event slopes decline from
0.278 to 0.244 eV. This curvature accompanies a changing event-conditioned
resistance mix: the TJ share falls from 85.25% at 800 K to 66.67% at 1100 K.
Grain-scale CV falls from 2.266 to 1.450 and stationary fraction from 66.51%
to 45.68%, while velocity--curvature `R²` remains below `7.1e-4` at every
temperature. Persistent TJ compatibility therefore preserves measurable
ensemble kinetics while producing strong, temperature-dependent intermittency.

The completed S2 stochastic shear-relaxation family has common `n=1.000`,
with growth `Q_app=0.101674 eV` (95% CI `[0.083684,0.118337]`, Arrhenius
`R²=0.98044`) and event-rate activation `0.196481 eV` (95% CI
`[0.186468,0.206294]`, `R²=0.99506`). Its local growth slope falls from
0.131 to 0.069 eV and local event slope from 0.221 to 0.157 eV. Across
800--1100 K, reverse-motion fraction rises from 56.94% to 72.00% while
velocity--curvature `R²` remains between 0.0014 and 0.0028. Velocity has a
negative correlation with stored internal shear (`-0.074` to `-0.085`),
consistent with the implemented negative energy derivative, and accumulated
signed shear-event strain rises from 18.15 to 29.03. This is direct production
evidence that reduced local memory can generate Qiu-like decorrelation and
reverse-curvature migration without a full-field solve.

The completed C2 climb-quota family retains common `n=1.000`, with growth
activation `0.331391 eV` (95% CI `[0.308080,0.353101]`, `R²=0.99273`) and
event-exposure activation `0.379942 eV` (95% CI `[0.376539,0.383542]`,
`R²=0.99492`). Local growth slopes decline 0.375→0.338→0.242 eV and
local event slopes decline 0.436→0.360→0.317 eV. CV falls 1.139→0.981,
reverse-curvature motion remains 67.63--79.94%, and velocity--curvature
`R²<=0.00221`. This mechanism-isolation curvature is opposite to FP-C2's
low-temperature shoulder and approximately 0.75 eV high-temperature slopes.

The completed C3 exchange-limited family also retains common `n=1.000`.
Growth activation is `0.166992 eV` (95% CI `[0.139985,0.189155]`,
`R²=0.95759`) and event-exposure activation is `0.220058 eV` (95% CI
`[0.217210,0.222619]`, `R²=0.99186`). Local growth slopes decline
0.243→0.130→0.097 eV and local event slopes decline
0.252→0.223→0.158 eV. Pairing this with FP-C3 shows that intrinsic
mobility reverses the event-slope curvature and raises the apparent scale,
while the common growth exponent remains at one.

The completed C4 transport-limited family retains common `n=1.000`, and its
growth and event activation energies agree: `0.264693 eV` (95% CI
`[0.241650,0.288602]`) and `0.264346 eV` (95% CI
`[0.262326,0.266883]`). The event Arrhenius fit is nearly linear
(`R²=0.99984`) with local slopes 0.270, 0.263, and 0.255 eV, whereas local
growth slopes rise 0.220→0.286→0.306 eV. CV falls 1.324→1.016 while
reverse-curvature motion remains 77.49--79.23%. This is the isolated
transport reference for the fully physical FP-C4 comparison.

The completed C5 serial-climb family retains common `n=1.000`, with growth
activation `0.258759 eV` (95% CI `[0.239275,0.283113]`) and event-exposure
activation `0.309896 eV` (95% CI `[0.307792,0.312735]`); both Arrhenius
`R²` values exceed 0.9995. Local growth/event slopes are
0.254→0.252→0.279 and 0.300→0.312→0.325 eV. Within completed serial
stages, expected residence remains 99.17--98.39% transport, while exchange
rises 0.745→1.437% and nucleation 0.087→0.168%. Across all observed event
types, transport falls 92.86→79.33% as generic activation rises
6.36→19.37%. These are event-conditioned expected-residence diagnostics,
not causal probabilities.

The completed P2 swept-area spatial control is identical at all four
temperatures and has no activated primitive events. Its common fit reaches the
upper search bound `n=6`, with `K=7388.828819` and a nominal growth activation
of only `9.9e-12 eV` (95% CI `[-0.0782,0.0857]`). This is intentionally
reported as bound-censored: the evidence identifies a scaling class at or
above the analyzed range, not a uniquely determined sixth-power law. Its
temperature invariance independently verifies that the measured behavior is
geometric rather than Arrhenius.

The completed SC3 full-field parallel-mode family has common `n=1.000` and
essentially no grain-growth activation: `Q_app=-0.007028 eV` (95% CI
`[-0.023510,0.007732]`, `R²=0.3326`). Its event-exposure activation is only
`0.005290 eV` (95% CI `[0.004939,0.005657]`, `R²=0.99675`), with rates
206.35→210.68 events per integrated GB-domain time. A dedicated streamed
occupation analysis over 66,720,712 completed modes nevertheless resolves a
small temperature-dependent coupling shift. Mean `|beta|` decreases
2.36756→2.35827 (0.392%; non-overlapping five-seed bootstrap intervals), easy-
shell occupation rises 31.716%→32.079%, and high-shell occupation falls
34.236%→33.953%. Thus Arrhenius occupation creates a measurable effective-
coupling temperature dependence even though the spectrum is fixed and the
coarse growth coefficient is temperature invariant.

The first completed fully physical family, FP-G2, retains common `n=1.000`
but has growth `Q_app=0.500831 eV` (95% CI `[0.479760,0.521590]`) and event-
exposure activation `0.531786 eV` (95% CI `[0.525019,0.539315]`). These are
far above mechanism-isolation G2 growth/event slopes of 0.143697/0.250122 eV.
Both fully physical Arrhenius fits are curved: adjacent growth slopes rise
0.351, 0.573, 0.639 eV and event slopes rise 0.352, 0.629, 0.680 eV. Thus the
global apparent barrier reflects coupled intrinsic mobility, evolving domain
exposure, and compatibility kinetics rather than a simple sum of input
barriers. CV falls from 1.776 to 1.340, Fano rises from 46.0 to 721.2, and the
common exponent remains unchanged across 800--1100 K.

The completed FP-T2 family gives common `n=1.000`, growth activation
`0.658646 eV` (95% CI `[0.629913,0.682759]`, `R²=0.99943`), and event-
exposure activation `0.582746 eV` (95% CI `[0.576193,0.589511]`,
`R²=0.98892`). Local growth slopes are 0.659, 0.624, and 0.719 eV, while
event slopes show a low-temperature transition from 0.445 to 0.658 and 0.693
eV. The event-conditioned TJ resistance share falls from 84.07% at 800 K to
66.64% at 1100 K. Despite this changing balance, the common exponent remains
one; CV falls from 2.262 to 1.429 and velocity--curvature `R²` remains below
`6.6e-4`.

The completed FP-S2 shear-memory family also retains common `n=1.000`.
Its growth activation is `0.479876 eV` (95% CI `[0.455151,0.505639]`,
`R²=0.99717`) and its event-exposure activation is `0.495505 eV` (95% CI
`[0.485472,0.504766]`, `R²=0.97869`). Local growth slopes rise
0.443→0.470→0.564 eV, while event slopes rise more sharply
0.336→0.575→0.638 eV. Reverse-motion fraction remains 63.64--70.75%,
velocity--curvature `R²` is at most 0.00379, and velocity correlates
negatively with stored shear (`-0.083` to `-0.080`) at every temperature.
Intrinsic Arrhenius mobility therefore shifts the apparent activation scale
without removing the shear-memory signatures seen in mechanism isolation.

The completed FP-C2 climb-quota family likewise retains common `n=1.000`.
Its growth activation is `0.675428 eV` (95% CI `[0.643993,0.701673]`,
`R²=0.99563`) and its event-exposure activation is `0.658582 eV` (95% CI
`[0.656294,0.660891]`, `R²=0.99077`). Both fits resolve a low-temperature
shoulder: local growth slopes are 0.572, 0.743, and 0.739 eV, while local
event slopes are 0.514, 0.746, and 0.761 eV. CV falls 1.159→0.983 as
temperature rises, reverse-curvature motion remains 70.33--77.65%, and
velocity--curvature `R²` never exceeds 0.00205. Thus the fully physical
climb reference remains intermittent and decorrelated while its apparent
barrier approaches the high-temperature climb scale.

The completed FP-C3 exchange-limited family retains common `n=1.000`, with
growth activation `0.541221 eV` (95% CI `[0.523330,0.558270]`,
`R²=0.99989`) and event-exposure activation `0.522374 eV` (95% CI
`[0.519825,0.524536]`, `R²=0.98871`). Its local growth slopes remain nearly
constant at 0.531, 0.559, and 0.529 eV, but local event slopes rise from
0.398 eV below 900 K to 0.591 and 0.621 eV above it. CV remains
0.991--1.086, reverse-curvature motion 69.42--72.86%, and
velocity--curvature `R²<=0.00158`. Exchange-controlled event kinetics thus
show a transition near 900 K without changing the coarse growth exponent.

The plotting path applies the same distinction: it retains event-level
Arrhenius and local-slope panels even when every growth coefficient is zero.

## 9. Intermittency and constrained optimization

The diagnostics include grain-level `A(t)`, `R(t)`, area/radius rates,
stationary fraction, CV, event-window Fano factor, waiting-time burstiness,
motion concentration in the top 1/5/10% intervals, burst increments/durations,
burst CCDF, per-entity event waiting times, primitive event types, active-domain
fraction, and velocity correlations with curvature, shear, free-volume deficit,
neighbor number, and simultaneous topological neighbors.

The completed 60-run jerkiness search varies correlation length, `K`, encounter
density, barrier, and packet renewal duration. Pareto ranking maximizes
CV/Fano/burstiness subject to positive growth and an exponent near 2 (or a CI
containing 2). None of the 12 five-seed candidates passes that strict scaling
gate, so no candidate is relabeled as physically admissible. The closest is
`JK-L24-K3-SPARSE`, with `n=1.188` (95% CI `[1.000,1.798]`), trajectory
CV `1.811`, and Fano factor `4.09`; its CI still excludes 2. At fixed 24-pixel
correlation length, persistent `K=1,3,5` raises trajectory CV from
`1.435` to `1.704` to `1.905` and waiting-time burstiness from `-0.075` to
`-0.018` to `0.025`, while event-count Fano falls from `53.28` to `38.94`
to `32.87`. Thus increasing `K` is not a scalar increase in “jerkiness”:
trajectory concentration rises while count overdispersion falls. The complete
ranking and immutable hashes are in
`results/campaigns/20260816T025109Z-c0b1b0f774/pareto_summary.csv` and
`results/validation/jerkiness_search_200_completion.json`.

## 10. Shear, climb, and resistance attribution

For serial climb, every primitive stage transition records its instantaneous
rate and hazard threshold. Analysis converts rates to expected residence times
and reports the nucleation/exchange/transport fraction of total expected
resistance. It also reports event-conditioned `1/r` distributions and normalized
shares for every observed primitive type, separating GB/shear-compatible,
TJ-compatible, nucleation, exchange, and transport passages. Those shares are
explicitly diagnostics at observed events, not causal probabilities; matched
blocked-time comparisons remain necessary for causal attribution. Event shear
and volumetric strain are summed separately. The exact 165-run mechanism
composite rejects a universal resistance rule: C5 is transport-dominated,
T1--T3 and explicit J1/J2 are TJ-dominated, while SC1--SC4 change qualitatively
with constraint architecture. The selected C5 temperature series further shows
that all-event transport residence falls from 92.86% to 79.33% over
800--1100 K as generic activation rises from 6.36% to 19.37%.

## 11. Scientific questions Q1--Q24

The final evidence map is explicit so a missing answer cannot be hidden in a
grouped narrative. Entries marked pending are not interpreted before their
immutable source campaign completes.

| Question | Required inference | Current evidence gate |
|---|---|---|
| Q1 | Can high-barrier compatibility events yield strong grain-scale jerkiness with smooth ensemble scaling? | Yes, with a scaling caveat. Corrected SC4 gives CV `1.172`, Fano `584.80`, and smooth ensemble fit `R²=0.9637`, but its `n=1` is at the scan bound and is not conventional parabolic scaling. Explicit J1/J2 show that stronger compatibility can instead arrest growth entirely. |
| Q2 | How does geometric encounter measure set the coarse-grained exponent? | It changes the scaling class. P1's GB-area analogue has `n=1`, TJ pinning P4 has `n=1.3935` (CI `[1.0000,1.7691]`), and swept-area spatial P2/P3 exceed the scan range at `n=6`. The selected P2 temperature control independently repeats the same upper-bound class at all four temperatures. At fixed G-series geometry, changing single-hit to multihit lowers `K` 0.17450→0.10734→0.07307 without changing the bound-limited `n=1`. |
| Q3 | How do single-hit Poisson and multihit kinetics differ? | Corrected explicit J1/J2 are both stagnant, but persistent multihit J2 raises CV from 3.231 to 4.242 and Fano from 291 to 3208 while lowering burstiness from 0.442 to 0.112 and increasing primitive events 19-fold. The effect depends on the diagnostic. |
| Q4 | Does increasing `K` raise or lower intermittency at fixed spatial correlation? | In the completed L24 series, K=1→3→5 raises trajectory CV (1.435→1.704→1.905) and burstiness (-0.075→-0.018→0.025), but lowers event-count Fano (53.28→38.94→32.87); the answer is metric dependent, consistent with multihit self-averaging of event counts. |
| Q5 | Can reduced shear memory reproduce Qiu-like velocity-curvature decorrelation without full elasticity? | Yes. In S2, velocity--curvature `R²` is only 0.0014--0.0028 and velocity correlates negatively with stored shear (`-0.074` to `-0.085`). FP-S2 independently retains `R²<=0.00379` and negative correlation (`-0.083` to `-0.080`) after adding intrinsic Arrhenius mobility. |
| Q6 | Can shear memory generate reverse-curvature migration? | Yes. S2 reverse-motion fraction is 56.94--72.00%; FP-S2 remains 63.64--70.75%. Negative velocity--stored-shear correlation in both series verifies reversal from mechanical back force. |
| Q7 | Can climb alone generate intermittency and velocity-curvature decorrelation? | Yes. The exact 900 K C1--C5 composite gives CV 0.994--1.178, reverse motion 70.4--80.2%, and velocity--curvature `R²=1.6e-4`--`8.2e-4`; selected C4 spans CV 1.016--1.324 and `R²<=0.00162` across temperature. FP-C2/FP-C3 independently give CV 0.977--1.159, reverse motion 69.42--77.65%, and `R²<=0.00205`. The decorrelation is substantial; intermittency is moderate compared with strict TJ cases. |
| Q8 | What resistance fractions arise from nucleation, exchange, transport, shear, and TJ compatibility? | Fractions are condition-specific, not global. Across 800--1100 K, selected C5 all-event residence changes from 92.86% to 79.33% transport, 0.70% to 1.16% exchange, 0.08% to 0.14% nucleation, and 6.36% to 19.37% generic activation; the exact 900 K composite value is 88.88% transport. T1--T3 are 78.04--78.92% TJ; J1/J2 are 97.07%/75.18% TJ. SC4 is 98.75% transport and 0.39% TJ. These are event-conditioned mean-`1/r` fractions, not causal probabilities. |
| Q9 | Are simultaneous shear and climb additive in residence time, strongly coupled, or dominated by one process? | Architecture-dependent. Independent-and SC1 (`K=0.07496`) follows slow C5 (`0.07226`) and is transport-dominated; mixed-event SC2 is much faster (`0.62315`); parallel SC3 is intermediate (`0.34784`); strict SC4 is again transport-dominated (`K=0.15522`). A single additive rule is rejected. |
| Q10 | How does apparent grain-growth activation compare with imposed microscopic barriers? | B1 recovers 0.45 eV exactly within uncertainty. Mechanism-isolation event/growth Q is 0.250122/0.143697 eV for G2, 0.266257/0.275837 eV for T2, 0.196481/0.101674 eV for S2, 0.379942/0.331391 eV for C2, 0.220058/0.166992 eV for C3, 0.264346/0.264693 eV for C4, and 0.309896/0.258759 eV for C5. Fully physical event/growth Q becomes 0.531786/0.500831 eV (G2), 0.582746/0.658646 eV (T2), 0.495505/0.479876 eV (S2), 0.658582/0.675428 eV (C2), and 0.522374/0.541221 eV (C3). Coarse-grained Q is emergent rather than a direct microscopic-barrier readout. |
| Q11 | Can similar apparent activation energies coexist with very different mobilities? | Yes provisionally: G2 event Q=0.2501 eV and T2 event Q=0.2663 eV are close, while their growth coefficients differ substantially (for example 0.10734 versus 0.06991 at 900 K) because their compatibility statistics differ. Final paired inference awaits all selected regimes. |
| Q12 | Under what conditions does physical stagnation occur? | Explicit TJ compatibility can cause physical stagnation even under enormous event activity: J1 changes characteristic radius by only 1.02%, while strict persistent J2 changes it by 0% across 3,500 steps. High event count is not evidence of coarsening. |
| Q13 | Which parameters maximize jerkiness without destroying realistic mean scaling? | None of the 12 completed search candidates passes the strict scaling gate. `JK-L24-K3-SPARSE` is closest (`n=1.188`, CI `[1.000,1.798]`, CV 1.811) but remains rejected; the high-barrier candidate has the largest CV (3.111) and also fails scaling. |
| Q14 | Can anisotropic stress/curvature selection from an isotropic discrete mode spectrum generate effective shear coupling? | Yes in the model sense. E0 generates signed shear 158.50 without feedback; full-field E1 gives 4256.82 and local-memory E2 gives 171.27. The magnitude is strongly backend-dependent, so quantitative equality to continuous Qiu coupling is not claimed. |
| Q15 | Can effective shear coupling acquire temperature dependence through Arrhenius mode occupation rather than an imposed coupling factor? | Yes, weakly. Across 66.7 million completed SC3 modes, mean `|beta|` falls 2.36756→2.35827 from 800→1100 K (0.392%, with non-overlapping endpoint bootstrap intervals) as easy-shell occupation rises and high-shell occupation falls. The spectrum itself is unchanged. |
| Q16 | How often do low-barrier modes fail explicit TJ compatibility? | Corrected J1 has 426,780 endpoint failures across 3,178,918 GB events (13.43% incidence), 60.55% nominally easy. J2 has 66,545,846 failures across 60,804,329 GB events (1.094 incidence; two endpoints may fail), 70.68% nominally easy. |
| Q17 | What barrier distribution is sampled during TJ compatibility failures? | Corrected J1 bare barriers have median 0.29 eV, 75th percentile 0.35 eV, and 99th percentile 0.59 eV; residual work broadens effective barriers from zero to 0.769 eV. |
| Q18 | Does explicit Burgers conservation yield long waits and abrupt TJ motion? | Yes. TJ activation contributes 97.07% of J1 and 75.18% of J2 event-conditioned expected residence; both are stagnant yet strongly intermittent (CV 3.231 and 4.242). |
| Q19 | How sensitive are growth and jerkiness to minimum Burgers magnitude and mode discreteness? | Minimum-Burgers and angular-quadrature component tests pass, but the 165-run matrix holds the mode library fixed. A matched production sensitivity is still required for a macroscopic answer. |
| Q20 | Can the geometric TJ surrogate reproduce explicit-mode scaling and intermittency? | No for the tested parameters. Explicit J1/J2 are stagnant with CV 3.231/4.242 and Fano 291/3208; surrogate J3 coarsens with `n=1`, `K=0.13380`, CV 1.449, and Fano 107.25. |
| Q21 | What accumulated strain fractions are shear, climb, and mixed? | In SC1--SC4, absolute signed strain components are 99.920--99.997% shear and 0.003--0.080% volumetric. SC2 assigns both components to one mixed event. These are dimensionless signed-strain component fractions, not energy or resistance fractions. |
| Q22 | Does apparent grain rotation emerge from event-wise tangential displacement? | Not directly identifiable: orientations are fixed and no grain-rigid-body rotation state is evolved. Event tangential strain is retained, but this remains a stated model limitation. |
| Q23 | At what temperature does the dominant difficult event cross from shear/TJ compatibility to climb/point-defect control? | C2/C3 supply declining isolated nucleation/exchange trends; C4 transport event slopes remain 0.270→0.255 eV, and serial C5 slopes rise only 0.300→0.325 eV while transport remains dominant. In FP-C2, slopes instead stabilize near 0.74--0.76 eV above 900 K; the C3/FP-C3 pair reverses event curvature from 0.252→0.158 to 0.398→0.621 eV. A coupled-process crossover still requires FP-C4/FP-C5 and FP-SC3. |
| Q24 | Does that crossover curve the Arrhenius plot or change the growth exponent? | C2/C3 and FP-G2/FP-T2/FP-S2/FP-C2/FP-C3 all retain common `n=1` while local Arrhenius slopes change. FP-C2 growth/event slopes change 0.572→0.743→0.739 and 0.514→0.746→0.761 eV, opposite C2's declining slopes; FP-C3 reverses C3's event curvature while its growth remains nearly linear. Coupled physics curves Arrhenius behavior without changing the fitted common exponent. |

## 12. Unexpected physics and preserved failures

The major resolved failures were diffuse-tail grain resurrection, physical-time
contamination by initial interface relaxation, finite-statistics exponent bias,
an early-window `n≈1` baseline, Python geometry throughput, raw diffuse-interface
reverse-motion noise, and the pre-climb wiring audit. Each has an immutable JSON
record or rejected summary. No result was deleted to improve agreement.

## 13. Known limitations

- Two-dimensional PF geometry with optional 2-D/3-D analytical mechanism
  scaling; no production 3-D field solver.
- Synthetic isotropic modes rather than a material-specific DSC catalog.
- Fixed grain orientations; no direct grain-rotation evolution.
- Reduced domain-length transport closure rather than spatial GB diffusion.
- The FFT mechanics backend is an independently implemented Qiu-type surrogate,
  not a bitwise port of the real-space line kernel.
- GB domains are deterministic pairwise chunks; topology is persistent by
  pair/triplet identity but does not solve a global crystallographic integer
  compatibility problem at every PF update.

## 14. Reproduction commands and result locations

```bash
python -m pip install -e '.[analysis,test]'
PYTHONPATH=src pytest -q
PYTHONPATH=src python scripts/run_validations.py
PYTHONPATH=src python scripts/validate_stochastic_engine.py
PYTHONPATH=src python scripts/run_qiu_regressions.py

PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/mechanism_scaling_200.yaml --processes 10
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/temperature_selected_200.yaml --processes 10
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/temperature_fully_physical_200.yaml --processes 10
PYTHONPATH=src python scripts/run_campaign.py \
  configs/production/jerkiness_search_200.yaml --processes 10

PYTHONPATH=src python scripts/compose_campaigns.py \
  results/campaigns/20260815T115547Z-af440f638c \
  results/campaigns/20260816T020915Z-e91786f3cd \
  results/campaigns/20260816T022756Z-876f84d9c9 \
  results/campaigns/20260816T160354Z-188fe232c6 \
  --expected-runs 165 --prefer-later-duplicates
PYTHONPATH=src python scripts/analyze_campaign.py \
  results/campaigns/20260817T042553Z-composite-28a9a5a185 \
  --output results/production_summaries/mechanism_composite_165_summary.csv
PYTHONPATH=src python scripts/plot_campaign.py \
  results/campaigns/20260817T042553Z-composite-28a9a5a185 \
  --summary results/production_summaries/mechanism_composite_165_summary.csv \
  --output results/plots/mechanism_composite_165

PYTHONPATH=src python scripts/analyze_campaign.py results/campaigns/<campaign> \
  --output results/production_summaries/<name>.csv
PYTHONPATH=src python scripts/analyze_campaign.py results/campaigns/<campaign> \
  --allow-incomplete --regime B1 --output <completed-regime-summary.csv>
PYTHONPATH=src python scripts/plot_campaign.py results/campaigns/<campaign> \
  --summary results/production_summaries/<name>.csv \
  --output results/plots/<name>
PYTHONPATH=src python scripts/audit_restart_artifacts.py results/campaigns/<campaign>
PYTHONPATH=src python scripts/pareto_campaign.py <jerkiness-summary.csv> \
  --output <jerkiness-pareto.csv>
PYTHONPATH=src python scripts/aggregate_summaries.py <summaries...> \
  --output results/final_mechanism_summary.csv
```

Campaign manifests are under `results/campaigns/`; compact summaries are under
`results/production_summaries/`; validation and failure records are under
`results/validation/`; final figures are under `results/plots/`. The final
machine-readable deliverable is `results/final_mechanism_summary.csv`.
