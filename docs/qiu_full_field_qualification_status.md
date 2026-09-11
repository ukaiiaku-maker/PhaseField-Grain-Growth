# QIU full-field qualification status

Last update: 2026-09-10, 384x384 corrected preflight accepted; legacy replay running; corrected production bundle next.

## Current source state

- Branch: `codex/qiu-full-field-qualification-v1`
- Audited base: `9f66c8d7a5a266687284d8da35aefbc6062808f7`
- Phase-0 commit: `66dce379ad6d9a0391ad76f4b653626891867d8a`
- Phase-1 instrumentation commit: `147141b54227e7c80c457e9d80c1a2ca3475bf61`
- Elastic-operator commit: `471412c`.
- Reference-identity commit: `13b43a6`.
- Local work-conjugate source commit: `e0a782b`.
- Energy-checked integration commit: `f0493bd`.
- Worktree: `/private/tmp/qiu-full-field-qualification-v1`
- Reduced-matrix commit: `9edbd0f32625604e0d43334842a7ad024dbbe7a3`.
- Scalable production runner commits: `6eb5977`, `8bb7837`.
- Historical production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`
- Historical QIU run: canonical, read-only, all integration-manifest hashes verified
- Current scientific decision: the historical QIU remains an unresolved non-self-similar transient and its backend is conclusively a legacy FFT eigenstrain surrogate, not the archived current-geometry Qiu reference formulation. The selected qualification backend is therefore honestly named `FFT_EIGENSTRAIN_V2`; production behavior remains pending.

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
- Ran 384x384 corrected-backend preflights from the canonical seed-5101 state,
  including exact trajectory comparisons across the performance revisions.
- Profiled the full production path and replaced phase-strided source/feedback
  scans, per-grain morphology masks, repeated rigid-translation label maps,
  repeated label reconstruction, and redundant dense PF rollback/history copies.

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
- The accepted corrected preflight is
  `results/qiu_full_field_qualification_20260910/20260910T190000Z-384-preflight-r5`.
  It completed three accepted 384x384, 800-starting-grain steps with zero coupled
  rejections, finite fields, equilibrium residual below `2.3e-13`, decreasing
  complete energy, restart checkpoint, per-step Parquet diagnostics, and movie
  frames. Its final state, energies, stresses, clipping statistics, and grain
  count are identical to the earlier unoptimized corrected preflights.
- Peak local RSS fell from 7.21 GB to 5.29 GB after eliminating redundant dense
  state copies. The observed local wall time is not a clean throughput estimate:
  two unrelated four-core DDD jobs were concurrently saturating the workstation.
  This is recorded as an operational limitation, not a numerical result.
- The production runner refuses cross-revision checkpoint resumes and verifies
  the canonical initial-state SHA before allocating the simulation.
- Latest focused validation after the production-memory changes: 58 passed,
  zero failures/errors. The subsequent complete repository suite passed all
  215 tests in 77.39 s with zero failures/errors; JUnit SHA-256 is
  `11d99dd625689ed7ba4f4a3f6fd4b3a94703bda27eeefa3ca3b186510ad4776c`.
- The full corrected seed-5101 production run is immutable HPC3 plan
  `20260911T003946Z-nogit-c53868`, source commit `8bb7837`, Slurm job
  `55932457`, using account `SDILLON1_LAB`. Its source-bundle SHA-256 is
  `e9926fc38ba41502d4deefadc57b344826924bfd2e4c0d75f4299ee1a6b8d5b0`.
  It is running concurrently with legacy replay job `55930486`, exactly at the
  two-full-worker ceiling. A delayed acknowledgement was reconciled against
  Slurm before retry; unused plan `20260911T004157Z-nogit-7dd0b1` remains only
  `PREPARED` and was never submitted.
- A fixed-physical-domain reduced spatial study at dx=1 and dx=0.5 is under
  `20260911T004500Z-spatial-convergence`. At matched t=4 it gives identical
  N=17 and population grain size, 0.78% interfacial/total-energy-density
  differences, 2.24% stress-p95 difference, and 1.40% mean-compactness
  difference. Elastic energy density differs by 3.73% (a failed 2% gate) while
  remaining about 4e-6 in absolute density; this is retained as a narrow
  spatial-convergence limitation rather than hidden.

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

Commit the spatial-convergence runner/result record and continue monitoring both
immutable HPC jobs. Retrieve and checksum terminal results locally. Start the
factor-two target refinement only after a worker slot is free.
