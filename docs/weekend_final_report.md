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
- Current fast suite: 85 tests passed in 47.99 seconds under the live
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

`results/validation/qiu_regression_benchmarks.json` passes. In the matched
24-grain polycrystal, full-field shear reduces velocity-curvature correlation
from 0.4358 to 0.3877 and raises active reverse-curvature motion from 12.70% to
30.85%. The four-grain benchmark raises active reverse motion from 3.09% to
11.38%. Eigenstrain, nonlocal stress, elastic feedback, and phase-field
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
Results remain **pending production inference** until the five-seed J1 and J2
repairs complete and the exact 165-condition composite is assembled.

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
Geometry-gate regressions now cover both paths. Event ledgers are buffered to
authoritative checkpoints, optionally gzip-compressed, and exactly truncated
on resume from checkpointed byte extents. The definitive SC4 repair completed
five seeds with zero tracebacks and gives `n=1.3248` (95% CI
`[1.0000, 2.2126]`), `K=0.1763`, jerkiness CV `1.162`, and Fano factor
`280.95`; final cross-regime interpretation still awaits the composite.

## 7. Growth-law results

**Pending production inference.** Final tables will report free and fixed-`n=2`
fits, objective topology windows, realization bootstrap intervals, exponent
profiles, local exponents, residual autocorrelation, radius-measure sensitivity,
population-band sensitivity, and source-manuscript Class-B/Class-C comparators.

## 8. Temperature and activation-energy campaigns

Two distinct four-temperature experiments are configured:

- mechanism isolation: 12 selected regimes at 800, 900, 1000, 1100 K, with
  intrinsic PF mobility held constant;
- fully physical: 11 `FP-*` regimes at the same temperatures, including a
  0.45 eV intrinsic-mobility barrier normalized to `M(900 K)≈4`.

Each condition uses the same five initial structures. A common exponent is fit
across temperature before extracting `K_n(T)`. Event activation is estimated
separately from primitive event counts divided by integrated GB-domain-time
exposure, including zero-event exposure. Adjacent-temperature slopes and
Arrhenius curvature are retained. A grain-growth activation fit is suppressed
when any series changes radius by less than 2%, rather than assigning a
spurious activation energy to stagnation. Results are **pending production
inference**.

## 9. Intermittency and constrained optimization

The diagnostics include grain-level `A(t)`, `R(t)`, area/radius rates,
stationary fraction, CV, event-window Fano factor, waiting-time burstiness,
motion concentration in the top 1/5/10% intervals, burst increments/durations,
burst CCDF, per-entity event waiting times, primitive event types, active-domain
fraction, and velocity correlations with curvature, shear, free-volume deficit,
neighbor number, and simultaneous topological neighbors.

The 60-run jerkiness search varies correlation length, `K`, encounter density,
barrier, and packet renewal duration. Pareto ranking maximizes CV/Fano/burstiness
subject to positive growth and an exponent near 2 (or a CI containing 2).
Results are **pending production inference**.

## 10. Shear, climb, and resistance attribution

For serial climb, every primitive stage transition records its instantaneous
rate and hazard threshold. Analysis converts rates to expected residence times
and reports the nucleation/exchange/transport fraction of total expected
resistance. Event shear and volumetric strain are summed separately. Comparative
mechanism results are **pending production inference**.

## 11. Scientific questions Q1--Q24

Final answers are **pending production inference**. They will explicitly cover:

- grain-scale jerkiness versus ensemble smoothness and encounter-size scaling;
- single-hit, persistent multihit, and finite packet-reset differences;
- local versus nonlocal shear memory, decorrelation, and reverse motion;
- climb-only intermittency and nucleation/exchange/transport resistance;
- independent versus mixed shear/climb residence times;
- event-level versus growth-level activation and mechanism crossover;
- matched apparent activation energies with different mobilities;
- physical stagnation and constrained jerkiness optima;
- stress/curvature selection from the isotropic discrete mode spectrum;
- low-mode TJ compatibility failures, Burgers residuals, and barrier sampling;
- sensitivity to mode discreteness and comparison with geometric surrogates;
- accumulated shear, climb, and mixed strain fractions.

Apparent grain rotation (Q22) is not directly evolved in the current fixed-
orientation PF state. Event tangential strain is retained, but a spatial
grain-rigid-body rotation estimator requires centroid-resolved event history and
will be reported as a limitation unless added without invalidating the campaign.

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

PYTHONPATH=src python scripts/analyze_campaign.py results/campaigns/<campaign> \
  --output results/production_summaries/<name>.csv
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
