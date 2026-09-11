# QIU full-field qualification status

Last update: 2026-09-10, Phase 6 reduced causal matrix complete; legacy replay running.

## Current source state

- Branch: `codex/qiu-full-field-qualification-v1`
- Audited base: `9f66c8d7a5a266687284d8da35aefbc6062808f7`
- Phase-0 commit: `66dce379ad6d9a0391ad76f4b653626891867d8a`
- Phase-1 instrumentation commit: `147141b54227e7c80c457e9d80c1a2ca3475bf61`
- Elastic-operator commit: `471412c`.
- Reference-identity commit: `13b43a6`.
- Worktree: `/private/tmp/qiu-full-field-qualification-v1`
- Historical production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`
- Historical QIU run: canonical, read-only, all integration-manifest hashes verified
- Current scientific decision: historical QIU remains an unresolved non-self-similar transient; source identity is provisionally `legacy FFT eigenstrain surrogate`, pending formal mathematical gates

## Commands completed

- Recorded `git status --short`, branch, HEAD, ten-commit log, and worktree list.
- Verified initial NPZ SHA-256 and every canonical QIU essential-file hash plus the exact 21-frame aggregate.
- Ran `PYTHONPATH=src:scripts /opt/anaconda3/bin/python -m pytest -vv --junitxml=results/validation/qiu_fix_baseline_tests.xml`.
- Ran `scripts/qiu_phase0_baseline.py` against the immutable historical QIU run.
- Inspected the archived `functions_2ref.py` pathways for `beta`, line-source construction, stress kernels, field extension, and PF feedback.
- Visually inspected the historical transition maps and terminal contact sheet.
- Added opt-in, per-step PF/mechanics/coupling/morphology diagnostics, partitioned
  Parquet streams, dense-field capture, guard capture, and restart state.
- Ran focused diagnostic invariance and restart tests, then the complete suite:
  `PYTHONPATH=src:scripts /opt/anaconda3/bin/python -m pytest -vv
  --junitxml=results/validation/qiu_fix_phase1_tests.xml`.

## Evidence obtained

- Baseline: 178 passed, 0 failed/errors/skipped; JUnit SHA-256 `67267a5132d122663af7843e2f7615b46579511b2649f43150fb4a9a73ba8da4`.
- The exact historical transition interval is 446 solver steps (9,800 to 10,246).
- Interfacial energy is 14,854 at step 9,900 and 64,319 at step 10,000.
- Boundary-shear maximum changes from 52.7 at step 9,000 to `3.24e13` at step 10,000.
- Disconnected-grain fraction changes from 0 to 0.827 over the same saved-frame interval.
- No complete pre-avalanche restart exists; the final checkpoint cannot recover the onset.
- The archived reference is a current-geometry, orientation-derived line-source model. The historical production backend is an accumulated point-eigenstrain FFT surrogate with constant beta and is not presently entitled to an unqualified QIU-reference label.
- Diagnostics enabled versus disabled produced bitwise-identical `eta`,
  eigenstrain, active-phase mask, time, step, and final grain count.
- Diagnostic restart produces one record per accepted step without duplicates;
  the recorder can also be enabled at an uninstrumented legacy checkpoint.
- Post-instrumentation complete suite: 181 passed, 0 failed/errors/skipped in
  101.42 s; JUnit SHA-256
  `ca826ef620a1babf67f74473b64276066a3559ac1206cbd670bdc3d6626d14da`.
- Immutable HPC3 legacy replay plan
  `20260910T232437Z-nogit-f4803a` was submitted as Slurm job `55930486`.
- The historical Fourier projection fails the required equilibrium gate: a
  fixed random 17x19 field gives normalized residual about `0.69`.
- `FFTEigenstrainV2` independently solves displacement equilibrium with an
  explicit traction-free zero mode and explicit plane-stress/plane-strain
  constitutive choices. Ten corrected-kernel gates and one explicit
  legacy-failure regression pass, including a dense global displacement
  oracle and three energy directional derivatives.
- Archived two-reference coupling factors match the directly extracted
  pristine `functions_2ref.py` function for eight synthetic orientation pairs.
- The unique local donor/receiver sweep passes planar translation, shrinking
  circle, multiple-neighbor, disconnected-pair, cancellation, rigid-
  translation, sign-reversal, topology-persistence, and three-grid refinement
  tests. Exact discrete swept-area relative error is zero in the manufactured
  refinement cases; finite-difference work-conjugacy error is below `1e-4`.
- Focused operator/source/reference/integration suite: 64 passed.
- Historical tracking data show the source-gain mechanism directly: mean/max
  tracked domains per physical grain pair increased from 1.376/3 at step 9800
  to 4.341/34 at step 10246. The old unweighted one-source-per-domain map is
  therefore partition dependent and self-amplifies after fragmentation.
- `fft_eigenstrain_v2` now limits explicit external forcing by a configured
  maximum phase increment and rejects/halves any coupled trial that increases
  complete interfacial-plus-elastic energy. Variable-dt restart is exact.
- Forced-limit tests verify timestep reduction, raw-increment bounds, complete-
  energy monotonicity, diagnostic/output-cadence invariance, and a monotonically
  decreasing matched-time solution error across four tightened targets.
- The accepted 48x48, 18-grain, 400-step causal matrix is under
  `results/qiu_full_field_qualification_20260910/20260910T171500Z-reduced-matrix-r2`.
  Seven trajectories (`L0,L1,K,S,T,KS,KST`) completed; `F` is the nonduplicated
  alias of `KST`, and the distinct unimplemented reference-line `R` is recorded
  as excluded rather than impersonated.
- At matched physical time 8, all controls retain 16 grains. L0 has max stress
  0.1297, equilibrium residual 0.7062, and source-work error -0.00405; corrected
  `KS/KST` has max stress 0.02376, residual about `1e-14`, and source-work error
  about `-9.5e-9`. The timestep-only `T` is identical to L0 because its external
  limit correctly remains inactive in this resolved interval.
- Halving legacy dt changes its matched-time stress from 0.1297 to 0.1940 even
  though total energy agrees within 0.001%; this rejects timestep convergence
  of the accumulated legacy source. Operator-only and source-only controls each
  reduce stress, while only the source correction closes the work residual.
- No reduced trajectory crosses a morphology guard. Population-only guards in
  the 18-grain cell reflect ordinary loss of more than two grains per 100 steps
  and are not classified as avalanches. Maximum compactness remains 1.741.
- A first matrix attempt at `20260910T164900Z-reduced-matrix` is preserved and
  excluded: its absolute clipping threshold confused normal double-obstacle
  halo projection with a clipping spike. The corrected guard uses a 20-step
  running baseline plus multiplicative and additive spike margins.

## Decisions recorded

- Preserve the legacy implementation and historical output unchanged.
- Do not treat saved frames as restart states.
- Replay deterministically from the immutable initial state.
- Keep superseded diagnostic analyses with explicit exclusion manifests.
- Test the whole-grain-area/per-boundary source-overcount hypothesis rather than assuming it is the sole cause.
- If mathematical equivalence to the archived reference is not established, qualify the corrected backend as `FFT_EIGENSTRAIN_V2`, not `QIU_REFERENCE_V2`.
- Treat the non-equilibrated legacy Green projection as a demonstrated defect,
  but do not yet assign sole causal responsibility for the historical
  avalanche; source overcounting and explicit lag remain live hypotheses.
- Select `plane_strain` for the corrected FFT surrogate to preserve the old
  surrogate's constitutive choice. This is not a claim that the archived Qiu
  line kernel is a plane-strain eigenstrain model.

## Next automatic action

Commit the reduced-control machinery and results. Launch seed-5101 corrected
production and a factor-two tighter target only after a 384x384 short preflight
passes memory, checkpoint, diagnostic, and performance checks. Continue
monitoring the immutable legacy replay without consuming a second full worker
until that preflight is accepted.
