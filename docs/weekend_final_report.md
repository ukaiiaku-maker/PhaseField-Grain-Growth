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
- Mechanism campaign: 146 valid completed conditions retained from
  `results/campaigns/20260815T115547Z-af440f638c`; whole-condition repair
  campaigns replace the rejected/interrupted SC4, P1, J1, J2, and J3 entries.
- Independent 256-square convergence campaign:
  `results/campaigns/20260815T083416Z-49423c71d0` (completed).
- Current fast suite: 102 tests passed in 89.14 seconds under the live
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
because it used the autocatalytic full-field sign. Final cross-regime
interpretation awaits the corrected full-field/TJ campaign and composite.

The corrected-sign SC4 replacement is now complete for five seeds. Its interim
fit gives `n=1.000` at the scan bound, `K=0.15522`, fit `R²=0.96374`, CV
`1.17168`, Fano `584.80`, reverse-motion fraction `0.65717`, and velocity--
curvature `R²=6.3e-5` across 1,254,987 primitive events. Within its serial climb
stages, expected `1/r` resistance is 99.15% transport, 0.74% exchange, and
0.11% nucleation. This is corrected interim evidence; matched cross-regime
inference still awaits the exact composite.

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
(95% CI `[1.43851,1.82699]`) over a 70.2% radius increase. Other selected and
fully physical results remain **pending production inference**.

The selected campaign has completed 180/240 runs: B1, G2, T2, S2, C2--C5,
and SC3 each have all 20 matched temperature/seed conditions. The remaining
P2/P3/final-family conditions continue with zero tracebacks. The fully physical
campaign is running at `results/campaigns/20260817T082103Z-01b68e843b` with
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
blocked-time comparisons remain necessary for Q8/Q9. Event shear and volumetric
strain are summed separately. Comparative mechanism results are **pending
production inference**.

## 11. Scientific questions Q1--Q24

The final evidence map is explicit so a missing answer cannot be hidden in a
grouped narrative. Entries marked pending are not interpreted before their
immutable source campaign completes.

| Question | Required inference | Current evidence gate |
|---|---|---|
| Q1 | Can high-barrier compatibility events yield strong grain-scale jerkiness with smooth ensemble scaling? | Corrected SC4 gives CV `1.172`, Fano `584.80`, and ensemble-fit `R²=0.9637`; final answer awaits the matched 165-condition composite. |
| Q2 | How does geometric encounter measure set the coarse-grained exponent? | Pending the G-series and explicit-versus-surrogate composite comparison. |
| Q3 | How do single-hit Poisson and multihit kinetics differ? | Corrected explicit J1/J2 are both stagnant, but persistent multihit J2 raises CV from 3.231 to 4.242 and Fano from 291 to 3208 while lowering burstiness from 0.442 to 0.112 and increasing primitive events 19-fold. The effect depends on the diagnostic. |
| Q4 | Does increasing `K` raise or lower intermittency at fixed spatial correlation? | In the completed L24 series, K=1→3→5 raises trajectory CV (1.435→1.704→1.905) and burstiness (-0.075→-0.018→0.025), but lowers event-count Fano (53.28→38.94→32.87); the answer is metric dependent, consistent with multihit self-averaging of event counts. |
| Q5 | Can reduced shear memory reproduce Qiu-like velocity-curvature decorrelation without full elasticity? | Yes. In the completed S2 series, velocity--curvature `R²` is only 0.0014--0.0028 across 800--1100 K, while velocity correlates negatively with stored shear (`-0.074` to `-0.085`). Matched full-field magnitude comparison still awaits the composite. |
| Q6 | Can shear memory generate reverse-curvature migration? | Yes. S2 reverse-motion fraction increases from 56.94% at 800 K to 72.00% at 1100 K; the negative velocity--stored-shear correlation verifies that reversal emerges from mechanical back force. |
| Q7 | Can climb alone generate intermittency and velocity-curvature decorrelation? | Pending C1--C5 composite results. |
| Q8 | What resistance fractions arise from nucleation, exchange, transport, shear, and TJ compatibility? | In corrected SC4, serial climb `1/r` resistance is 99.15% transport, 0.74% exchange, and 0.11% nucleation; final shear/TJ attribution awaits C/J/SC aggregation. |
| Q9 | Are simultaneous shear and climb additive in residence time, strongly coupled, or dominated by one process? | Corrected SC4 is transport-dominated in observed-event resistance, but causal comparison remains pending the SC1--SC4 matched composite. |
| Q10 | How does apparent grain-growth activation compare with imposed microscopic barriers? | B1 recovers 0.45 eV exactly within uncertainty. G2 gives event/growth Q of 0.250122/0.143697 eV; T2 gives 0.266257/0.275837 eV; S2 gives 0.196481/0.101674 eV. Coarse-grained Q is an emergent observable and can be far below the event-level slope. |
| Q11 | Can similar apparent activation energies coexist with very different mobilities? | Yes provisionally: G2 event Q=0.2501 eV and T2 event Q=0.2663 eV are close, while their growth coefficients differ substantially (for example 0.10734 versus 0.06991 at 900 K) because their compatibility statistics differ. Final paired inference awaits all selected regimes. |
| Q12 | Under what conditions does physical stagnation occur? | Explicit TJ compatibility can cause physical stagnation even under enormous event activity: J1 changes characteristic radius by only 1.02%, while strict persistent J2 changes it by 0% across 3,500 steps. High event count is not evidence of coarsening. |
| Q13 | Which parameters maximize jerkiness without destroying realistic mean scaling? | None of the 12 completed search candidates passes the strict scaling gate. `JK-L24-K3-SPARSE` is closest (`n=1.188`, CI `[1.000,1.798]`, CV 1.811) but remains rejected; the high-barrier candidate has the largest CV (3.111) and also fails scaling. |
| Q14 | Can anisotropic stress/curvature selection from an isotropic discrete mode spectrum generate effective shear coupling? | Mode-selection regressions pass; production magnitude awaits E/SC results. |
| Q15 | Can effective shear coupling acquire temperature dependence through Arrhenius mode occupation rather than an imposed coupling factor? | Pending event-resolved temperature series. |
| Q16 | How often do low-barrier modes fail explicit TJ compatibility? | Corrected J1 has 426,780 endpoint failures across 3,178,918 GB events (13.43% incidence), 60.55% nominally easy. J2 has 66,545,846 failures across 60,804,329 GB events (1.094 incidence; two endpoints may fail), 70.68% nominally easy. |
| Q17 | What barrier distribution is sampled during TJ compatibility failures? | Corrected J1 bare barriers have median 0.29 eV, 75th percentile 0.35 eV, and 99th percentile 0.59 eV; residual work broadens effective barriers from zero to 0.769 eV. |
| Q18 | Does explicit Burgers conservation yield long waits and abrupt TJ motion? | Yes. TJ activation contributes 97.07% of J1 and 75.18% of J2 event-conditioned expected residence; both are stagnant yet strongly intermittent (CV 3.231 and 4.242). |
| Q19 | How sensitive are growth and jerkiness to minimum Burgers magnitude and mode discreteness? | Discrete minimum-Burgers tests pass; production sensitivity awaits matched regime fits. |
| Q20 | Can the geometric TJ surrogate reproduce explicit-mode scaling and intermittency? | J3 replacement is complete; final answer awaits J1/J2 and the composite. |
| Q21 | What accumulated strain fractions are shear, climb, and mixed? | Signed shear and volumetric event sums are available; final fractions await SC aggregation. |
| Q22 | Does apparent grain rotation emerge from event-wise tangential displacement? | Not directly identifiable: orientations are fixed and no grain-rigid-body rotation state is evolved. Event tangential strain is retained, but this remains a stated model limitation. |
| Q23 | At what temperature does the dominant difficult event cross from shear/TJ compatibility to climb/point-defect control? | Pending the fully physical temperature campaign. |
| Q24 | Does that crossover curve the Arrhenius plot or change the growth exponent? | T2 retains common `n=1` but shows local growth Q from 0.238 to 0.317 eV as its event-conditioned TJ resistance falls from 85.25% to 66.67%; full coupled-regime inference remains pending. |

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
