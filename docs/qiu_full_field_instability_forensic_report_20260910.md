# QIU full-field instability forensic report

Status: the legacy-forensic mechanism is closed; the independent corrected
full-scale production matrix is still running. This file must not be cited as
a final `FFT_EIGENSTRAIN_V2` production qualification until
`qualification_decision.json` is issued.

## Executive finding

The historical Phase-1 `QIU-T900-s5101` trajectory is a rejected,
non-self-similar legacy-surrogate transient. Its 494-to-99 grain transition
occurred from solver step 9,800 to 10,246: 446 steps, not 246. Nothing in this
work modifies or supersedes that canonical artifact.

The old backend is not the archived Qiu reference model. The archive rebuilds
orientation-dependent line-disconnection sources from the current GB geometry;
the historical repository backend accumulates point eigenstrain at one tracked
boundary representative and applies a different Fourier projection. The
model-identity classification is therefore `QIU_REFERENCE_MODEL_MISIDENTIFIED`.
The corrected model selected for qualification is explicitly
`FFT_EIGENSTRAIN_V2`, a periodic accumulated-eigenstrain surrogate, not a
faithful Qiu port.

The terminal replay and clean source/work continuation close the immediate
mechanism: boundary `gb:343-357:1` encounters an exactly zero radius-one
velocity-stencil denominator; the `1e-14` hard floor produces displacement
`2.5367e13`, source norm `6.2781e12`, and an elastic-energy jump of
`1.9707e25`. Stress explosion follows, the first disconnection occurs at step
9673, and the visible morphology/population avalanche occurs later. No further
production replay of this legacy failure is warranted. The separately defined
`QIU_SI_REFERENCE` reproduction is documented in
`docs/qiu_si_reference_reproduction_status.md`.

Three concrete defects are independently demonstrated in the old surrogate:

1. Its Fourier projection does not satisfy mechanical equilibrium (normalized
   residual about 0.69 in the fixed oracle case).
2. Whole-grain area motion is reused once per tracked GB domain and deposited
   at one point. Historical fragmentation raised the mean/max tracked domains
   per physical pair from 1.376/3 at step 9,800 to 4.341/34 at step 10,246,
   multiplying the source precisely as morphology fragmented.
3. The accumulated legacy source is timestep-sensitive: halving dt changes
   matched-time maximum stress by about 50% in the reduced causal matrix.

These establish defects, and the source/work continuation now identifies the
numerical source singularity that initiates the historical late avalanche.
Full corrected/refined trajectories answer the separate question of whether
the independent `FFT_EIGENSTRAIN_V2` model has any qualified transition.

## Immutable provenance

- Audited base: `9f66c8d7a5a266687284d8da35aefbc6062808f7`.
- Historical production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`.
- Qualification branch: `codex/qiu-full-field-qualification-v1`.
- Canonical initial-state SHA-256:
  `106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`.
- Canonical historical campaign remains under
  `results/long_time_kinetics_900K_20260824/20260825T012949Z-6d83ee5c82`.
- New local evidence root:
  `results/qiu_full_field_qualification_20260910`.

The baseline suite passed 178 tests. After read-only instrumentation it passed
181. The complete suite at corrected checkpoint step 2000 passed all 234
collected tests in 60.67 s, with no failure, error, or skip; JUnit SHA-256 is
`f959e6452185a76ac34e44dd6703c370f48bb09822a9619c6b3d18b953c9a93b`.
The subsequent raw-capture/transition-classification analyzer change passes 12
focused analyzer and renderer tests; the final complete suite remains reserved
for the terminal evidence revision.

The latest split-workstream revision passes the complete 251-test suite in
51.26 seconds with zero failure, error, or skip. Its JUnit SHA-256 is
`3e970b0a4220c5a9da20f7fe0c1da72ccfad77f142a1b1718a52d0ea393eadff`.
This software result does not promote either running corrected trajectory to a
terminal scientific result.

## Historical transition

The saved historical series shows interfacial energy increasing from 14,854 at
step 9,900 to 64,319 at step 10,000. Maximum tracked boundary shear rises from
52.7 at step 9,000 to about 3.24e13 at step 10,000, while the disconnected-grain
fraction changes from zero to 0.827. The available frame archives are not
restart states: they lack complete PF, active-phase, eigenstrain, and
bookkeeping state. The only defensible reconstruction is deterministic replay
from the hashed initial condition.

That immutable replay is HPC3 run `20260910T232437Z-nogit-f4803a`, Slurm job
`55930486`, commit `147141b`. Its 17 saved frames through step 8,000 are
byte-identical to the canonical historical frames. That source revision still
has the original absolute clipping marker, which is oversensitive to ordinary
double-obstacle projection, but its wrapper explicitly continues after a
marker. The first marker itself is excluded while the per-step scalars and
dense fields continue through 10,246. The replay has now closed checkpoint
step 9,500/t=380/N=495. Its 18 common frames through step 9,000 are
byte-identical to history, and the predicted clipping-only marker at step 9,001
has zero extinctions and compactness mean/max 1.337/1.448. All 500 exact rows
from 9,001 through 9,500 have N=495, zero extinctions, zero disconnected
grains, and no compactness threshold crossing. Nevertheless, accepted total
energy first rises at step 9,002 and rises on 311 of 499 differences. From the
first to last 20-step medians, source-increment L2 rises 1.149 to 37.553,
stress linf 163.797 to 172.523, eigenstrain linf 150.401 to 162.560, and
elastic energy by 6,480 while interfacial energy falls by 12.7. Thus the
legacy elastic/source positive feedback demonstrably precedes population and
morphology failure rather than being caused by the extinction cascade.

The first 16 scalar rows also revealed a read-only instrumentation omission in
the production closure override: source-energy/work values are undefined and
no per-boundary rows are emitted, although trajectory fields remain valid.
Commit `bebd53b` adds those recorder hooks without changing any numerical array;
a paired full-size two-step smoke is bitwise trajectory-identical and produces
finite source/work plus 4,078 boundary rows. Immutable prepared continuation
`20260911T111913Z-nogit-06fc67` resumes the clean exact step-9,000 checkpoint
and will provide the required terminal source/work record after a slot frees.

## Corrected model

For each nonzero Fourier mode the corrected solver minimizes elastic energy over
compatible displacement strain and enforces

\[
k_j C_{ijkl}(\epsilon^c_{kl}-\epsilon^*_{kl})=0.
\]

The zero mode uses traction-free mean strain. Plane strain is explicit and is
selected only to preserve the historical surrogate's constitutive choice; it
is not attributed to the archived Qiu line kernel. Tests cover equilibrium,
symmetry, linearity, translations, uniform and manufactured compatible fields,
nonnegative energy, three energy derivatives, a dense global displacement
oracle, and a single Fourier mode. Corrected residuals are below 1e-10.

The unique local donor/receiver transfer is

\[
q_{i\to j}=\frac{(-\Delta\eta_i)_+(\Delta\eta_j)_+}
{\sum_k(\Delta\eta_k)_+},\qquad
B_{ij}=\beta_{ij}\operatorname{sym}(t\otimes n).
\]

The same `B` deposits eigenstrain and returns stress work to the phase field,
so `dE_el/dq=-sigma:B` in the discrete implementation. Manufactured planar,
circular, multineighbor, disconnected-pair, cancellation, rigid-translation,
direction-reversal, topology-persistence, and grid-refinement cases pass. The
source-integral error is zero in the exact discrete refinements and the
finite-difference work error is below 1e-4 relative.

Accepted coupled steps use the minimum intrinsic/external timestep and are
rejected and halved if complete interfacial-plus-elastic energy increases. No
stress cap, eigenstrain deletion, artificial smoothing, modulus reduction, or
extinction-threshold change is used.

## Causal matrix

The accepted 48x48, 18-grain, 400-step matrix is
`20260910T171500Z-reduced-matrix-r2`. At matched t=8 all controls have N=16.

| Variant | Change | max stress | equilibrium residual | source-work error |
|---|---|---:|---:|---:|
| L0 | legacy | 0.129657 | 0.7062 | -0.00405 |
| L1 | legacy, dt/2 | 0.193962 | 0.7068 | -0.00241 |
| K | corrected operator | 0.096608 | ~8e-15 | -0.00607 |
| S | corrected source/force | 0.037754 | 0.7054 | ~2e-18 |
| KS/KST | corrected operator and source | 0.023760 | ~1e-14 | ~-9.5e-9 |

`T` is identical to L0 because its external bound is inactive in this mild
interval. This is correct limiting behavior, not evidence that the limiter is
unnecessary. Forced-drive tests independently activate it and converge across
four targets.

The longer reduced target=0.02/0.01 pair reaches N=5 at step 2,405/t=96.2
without a guard, rejection, nonfinite value, energy increase, or morphology
failure. The trajectories are bitwise identical because both external limits
remain far above dt=0.04. It passes observable comparisons but is explicitly
inconclusive as target-refinement evidence; production refinement is required.

## Spatial evidence

The accepted fixed-physical-domain dx=1, 0.5, and 0.25 comparison preserves
domain size, interface width, initial grain density, and normalized seed
coordinates. At matched t=4 all grids have N=17 and identical physical G. For
dx=0.5 versus 0.25, differences are 0.22% for interfacial/total energy density,
2.77% for stress p95, 0.14% for mean compactness, and 1.48% for compactness p95.
Elastic-energy density improves from a 3.73% coarse-pair difference to 2.59%
on the finest pair, still narrowly failing the requested 2% sensitivity
threshold despite its small absolute value (~4e-6). The trend is convergent,
but the narrow limitation remains open rather than being rounded into a pass.

## Full-scale campaign state

The corrected seed-5101 production run is immutable HPC3 plan
`20260911T003946Z-nogit-c53868`, Slurm `55932457`, source `8bb7837`, with a
10-day allocation and 16 GB memory. Its source bundle SHA-256 is
`e9926fc38ba41502d4deefadc57b344826924bfd2e4c0d75f4299ee1a6b8d5b0`.
The latest locally retrieved nonterminal snapshot is the atomic step
2000/t=80/N=317 state. It contains all 2,000 contiguous scalar rows, exact
tracker prefixes, 11 compact frames, and full fields at steps 500, 1000, 1500,
and 2000. Complete energy decreases every step; maximum equilibrium residual
is `3.42e-12`; maximum absolute source-work error is `1.42e-7`; there are no
rejected steps, no disconnected grains, and no guard capture. The external
limit remains inactive (minimum finite limit 4.03 versus dt 0.04). Population
loss slows from 41 to 27 to 20 grains across the last three 250-step stages.
Its 13-frame preview and initial/midpoint/latest contact sheet show connected,
substantially polygonal coarsening without runaway elongation. Matching
local/remote recovery-archive SHA-256 is
`2a8e3d1ab77404c5efb0aef778496fb2f9ea86bb498480f53b7d00845ad57812`.
This snapshot is recovery/status evidence only and is not pooled as a completed
run.

The two active Slurm jobs are the only full-scale workers. The corrected target
refinement is running. The two additional seeds remain prepared but will not
start until both seed-5101 jobs are terminal, locally retrieved, and pass the
combined mechanics/work/energy/timestep gate. Terminal movies, final comparison
plots, and the decision JSON remain pending terminal evidence.

## Current decision boundary

The reference-identity decision is closed: the historical backend was
misidentified, and future inference from the corrected surrogate must use the
`FFT_EIGENSTRAIN_V2` name and production config. The legacy numerical trigger
and downstream causal ordering are also closed. Whether the independent
corrected model develops any qualified analogous transition remains open and
will not be inferred from small systems or partial runs.
