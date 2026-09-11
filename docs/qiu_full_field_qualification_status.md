# QIU full-field qualification status

Last update: 2026-09-11 05:02 PDT; legacy and corrected seed-5101 full-size jobs running on HPC3; corrected recovery audited through atomic checkpoint 1500; closure-specific source/work continuation remains first in the prepared queue.

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
- Spatial-convergence commit: `c9813b5`.
- Qualification-movie renderer commit: `847eb06`.
- Reduced long-time target study/analysis commits: `41859ec`, `f04cf8f`.
- Live report and handoff commit: `3b06bd5`.
- Paired-seed preparation commit: `4ef9390`.
- Three-grid spatial-convergence commit: `ae99d34`.
- Staged-source provenance fix: `a236192`.
- Matched-time/progress qualification analysis: `78da5ea`.
- Actual factor-two refinement runner: `a173624`.
- Bounded legacy-continuation/guard commit: `abde6e0`.
- Curated transition-contact-sheet commit: `167c30c`.
- Closure-specific source/work instrumentation commit: `bebd53b`.
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
- After adding only qualification analysis/movie tools and the configurable
  reduced-cell population guard, the expanded complete suite passed 216/216 in
  76.82 s; current JUnit SHA-256 is
  `33cc6789dc661d66bc58ac76fcde7e003ef27daf536a38ebd81dda9046e82b8a`.
- The full corrected seed-5101 production run is immutable HPC3 plan
  `20260911T003946Z-nogit-c53868`, source commit `8bb7837`, Slurm job
  `55932457`, using account `SDILLON1_LAB`. Its source-bundle SHA-256 is
  `e9926fc38ba41502d4deefadc57b344826924bfd2e4c0d75f4299ee1a6b8d5b0`.
  It is running concurrently with legacy replay job `55930486`, exactly at the
  two-full-worker ceiling. A delayed acknowledgement was reconciled against
  Slurm before retry; unused plan `20260911T004157Z-nogit-7dd0b1` remains only
  `PREPARED` and was never submitted.
- A nonintrusive 18:30 PDT live audit found legacy replay at step 1500/t=60
  and corrected seed 5101 at step 96/t=3.84/N=760. The corrected trajectory
  remains finite with zero coupled rejections; its latest equilibrium residual
  is `8.65e-14`, source-work error is `-3.49e-8`, and external dt limit is
  11.10 versus used dt 0.04. The live outputs continue to grow on node-local
  Slurm scratch and neither job has a terminal marker yet.
- The next atomically closed corrected milestone contains 208 steps through
  t=8.32/N=705. It has zero complete-energy increases and coupled rejections;
  all guard-critical fields are finite; maximum equilibrium residual is
  `2.16e-12`; maximum absolute source-work error is `1.33e-7`; maximum
  compactness is 1.795; and no disconnected grain appears. The largest
  100-step loss is 52 (about 7% of its preceding population, below the 10%
  guard), while clipping stays in the narrow 0.263--0.273 running baseline.
  Durable nonterminal audit: `hpc_live/20260911T024100Z-corrected-step208-audit.json`.
- The later corrected recovery snapshot contains 25 closed Parquet parts
  through step 400/t=16/N=617 and compact movie frames at steps 0, 200, and
  400. It retains zero complete-energy increases and coupled rejections, finite
  guard-critical fields, maximum equilibrium residual `3.42e-12`, maximum
  absolute source-work error `1.33e-7`, maximum compactness 1.795, and zero
  disconnected grains. The minimum external limit is 6.25 versus used dt 0.04;
  clipping stays within 0.254--0.273. Archive and audit are durable under
  `hpc_live/20260911T040500Z-corrected`, SHA-256
  `a06753ec606e29860b655624cc62195a572733c3c3e3648daa54ce7f82876bbd`,
  and remain explicitly nonterminal/nonpoolable.
- The corrected trajectory subsequently reached its first regular restart
  checkpoint at step 500/t=20/N=587. Its new durable archive contains all 500
  contiguous scalar rows in 32 closed Parquet parts, the restart state, compact
  frames through step 400, and the first cadence-dense full field at step 500.
  The audit still finds zero complete-energy increases and coupled rejections,
  finite critical fields, maximum equilibrium residual `3.42e-12`, maximum
  absolute source-work error `1.33e-7`, maximum compactness 1.795, and zero
  disconnected grains. Its external limit remains inactive (minimum 6.13
  versus used dt 0.04). Durable path:
  `hpc_live/20260911T044539Z-corrected`; matching local/remote archive SHA-256
  `9ec285d7213768c446ccc34b5996f352b510160c3f049b4cc5cf43c1b50d40fe`.
  This snapshot is explicitly nonterminal and nonpoolable.
- Corrected seed 5101 later reached step 1000/t=40/N=460. Its durable recovery
  archive contains all 1000 contiguous scalar rows in 64 closed parts, six
  compact frames through step 1000, and dense full fields at steps 500 and
  1000. The formal audit still finds zero complete-energy increases and coupled
  rejections, no diagnostic capture, finite critical fields, maximum
  equilibrium residual `3.42e-12`, maximum absolute source-work error
  `1.42e-7`, maximum compactness 1.795, and zero disconnected grains. The
  minimum external limit is 4.28 versus used dt 0.04. Durable path:
  `hpc_live/20260911T082618Z-corrected`; matching local/remote archive SHA-256
  `f1fdd7f98987a2e8373e4ef6c6ad7932d788234e19e3b53cb28f9eec089b5bfc`.
  This snapshot remains explicitly nonterminal and nonpoolable.
- The recovered corrected step-1000 fields were passed through the final movie
  pipeline. The resulting seven-frame GIF, SHA-indexed frame table, JSON
  metadata, and contact sheet are retained beside the snapshot as
  `corrected-step1000-preview.*`. The selector correctly used neutral
  initial/midpoint/terminal roles because no diagnostic capture exists; visual
  inspection shows smooth coarsening and spatially resolved, finite stress and
  eigenstrain fields. This is a renderer/data-retention check, not a terminal
  production movie or scientific endpoint.
- Corrected seed 5101 has now reached atomic checkpoint step 1500/t=60/N=364.
  Its durable v2 recovery contains all 1,500 contiguous, unique scalar rows in
  96 closed Parquet parts, eight compact frames through step 1400, three dense
  fields through step 1500, and the exact tracker byte offsets from the restart
  record. The local and remote archive SHA-256 is
  `fbf755562df4f866eb80f7eadae6bcc4271d4009940b4de2cbad515e03dc17a5`.
  The audit still finds zero complete-energy increases, zero coupled
  rejections, no capture, finite critical fields, equilibrium residual at most
  `3.42e-12`, absolute source-work error at most `1.42e-7`, maximum
  compactness 1.954, and no disconnected grain. The minimum external limit is
  4.28 versus used dt 0.04. Durable path:
  `hpc_live/20260911T115520Z-corrected`; the snapshot remains explicitly
  nonterminal and nonpoolable.
- The first archive assembled at that milestone is preserved but excluded: its
  filename-based Parquet bound stopped at step 1480 while the checkpoint was
  1500. `superseded_archive_manifest.json` records the mismatch and points to
  the range-verified v2 archive; no trajectory or live output was affected.
- The legacy step-1500 checkpoint plus all then-available movie frames were
  copied without pausing the solver to local recovery archive
  `hpc_live/20260911T013800Z-legacy`, SHA-256
  `b8ae0b4205a1074dacef9cef7633442ea7b6a6ba743145662952bdb34a5e182a`.
- After `/pub` metadata service recovered, a newer complete recovery archive
  was copied and independently rehashed at legacy checkpoint step 2500/t=100,
  with movie frames through step 2400. Durable path:
  `hpc_live/20260911T030100Z-legacy`; archive SHA-256
  `7b2c7ddb7481e7ca63c20134242368bd8a709f3fa35451f2cc30bb97dfcbf5d8`.
  Its manifest remains `pooling_allowed=false` until terminal retrieval.
- At 21:20 PDT, the legacy replay reached a closed step-3500/t=140 checkpoint.
  The checkpoint and all 12 then-available evolution frames (latest step 3000)
  were copied without pausing the solver, rehashed on HPC3 and locally, and
  retained under `hpc_live/20260911T041940Z-legacy`. Archive SHA-256 is
  `d58b052f5f8ff46112deee886f0f669013024b531851e22bcaaf0aa5c6e2ec5b`;
  the recovery manifest remains explicitly nonterminal and nonpoolable.
- At 23:05 PDT, legacy reached step 5000/t=200. Its complete checkpoint and
  14 saved evolution frames through step 5000 were copied and independently
  rehashed on HPC3 and locally under `hpc_live/20260911T060150Z-legacy`.
  Matching archive SHA-256 is
  `26c26bb615846c6c91e4d0b99091bb53ad31f1af5d24b547e8fbc895b7f422ce`;
  this recovery copy remains explicitly nonterminal and nonpoolable.
- At 01:45 PDT, legacy reached step 7000/t=280. Its complete checkpoint and
  16 saved evolution frames through step 7000 were copied and independently
  rehashed on HPC3 and locally under `hpc_live/20260911T083849Z-legacy`.
  Matching archive SHA-256 is
  `89e04dcd7b0f95b2e950b52f001b6dcf8d922ffa148f1be100d4f6462703bf6e`;
  this recovery copy remains explicitly nonterminal and nonpoolable.
- At 03:01 PDT, legacy reached step 8000/t=320. A further complete checkpoint
  and all 17 saved evolution frames through step 8000 were copied, independently
  rehashed on HPC3 and locally, and retained under
  `hpc_live/20260911T100125Z-legacy`. Matching archive SHA-256 is
  `d0b093ec1d3e3a4f5656d9009c381609c61925fb70a9061fff76b25dbf1dbbe5`;
  this recovery copy remains explicitly nonterminal and nonpoolable.
- All 17 replay frames from step 0 through 8000 are byte-for-byte identical to
  their canonical historical counterparts. An independent NPZ load also found
  identical key sets and bitwise-identical arrays in every frame. The durable
  machine-readable comparison is
  `hpc_live/20260911T100125Z-legacy/legacy_replay_equivalence_through_step8000.json`,
  SHA-256 `7277391fdd27816802a7ec4c86ab2578988024531ee80583c890d89ea9aa962d`.
  After removing only the intentionally different run-ID column, all 22,416
  checkpoint-bounded grain-track lines and 89,187 boundary-track lines are also
  byte-identical to the canonical records through step 8000. This establishes
  exact deterministic historical reproduction before the critical window,
  while leaving terminal avalanche reproduction pending.
- The replay reached the atomic step-9000/t=360 checkpoint, and its saved
  step-9000 frame SHA-256
  `7b5f88ab88daf286a5c8bd1cdd76de18965969c4aaecea7653d0a7af949bbaa8`
  is byte-identical to the canonical historical frame. As predicted, the old
  absolute clipping marker fired at step 9001 on clipped fraction 0.2414, with
  zero extinctions, compactness mean/max 1.337/1.448, and all critical fields
  finite. That marker is excluded from physical-transition inference; the
  wrapper continues and is writing per-step scalars and fields.
- A concurrency-safe recovery assembled from the atomic step-9000 checkpoint
  plus explicitly closed step-9001 scalar/field files is durable under
  `hpc_live/20260911T110245Z-legacy`. Its matching local/remote archive SHA-256
  is `20861606468eb6be39ae4012d7872e31ffd1e092fa7a71893586050e13ab3cea`;
  the step-9001 audit SHA-256 is
  `56c46c69da1d53528f46eaa91292445e9fdabfe7c9281856d3289bef865a00bf`.
  This snapshot remains nonterminal and nonpoolable.
- The first 16 active-replay rows (steps 9001--9016) have finite morphology,
  stress, eigenstrain, and energy, with N fixed at 495 and no extinction. They
  also exposed that `MigrationClosureSimulation._update_physics` overrides the
  instrumented base method: its active source revision never seeds the
  per-source energy baseline or per-boundary recorder, leaving
  `source_elastic_energy_change`/`source_work_error` undefined and the boundary
  stream empty. This is a read-only evidence defect, not a trajectory change.
  The active replay remains authoritative for exact trajectory/dense fields but
  is insufficient by itself for the source/work gate.
- Commit `bebd53b` adds the missing closure-specific recorder hooks: it resets
  only the read-only per-step source accumulator, seeds the pre-source energy,
  and records the unchanged midpoint source event and old-stress work. A paired
  full-size step-9000-to-9002 smoke test is bitwise identical in all seven
  numerical checkpoint arrays and identical in all nondiagnostic state. The
  corrected stream has finite per-step source energy/work and 4,078 finite
  per-boundary rows; its large work mismatch is exposed rather than corrected.
  Evidence JSON SHA-256:
  `1a8b4b7d8a553e2984f5739170c773ef83bc633a28e10484c8f5aa0ab0763fa7`.
- The complete suite after this read-only correction passes 233/233 in 38.76 s
  with zero failures/errors/skips. JUnit SHA-256 is
  `c6d732402d0fab27fd9dfdfbbceae370be7e19fe6ec82dd7b005f08f441da2dd`.
- The active legacy source commit `147141b` predates the relative clipping-guard
  correction and therefore still treats ordinary double-obstacle clipping
  (about 0.24 at a two-step step-8000 continuation smoke test) as a trigger at
  the absolute 0.02 threshold. Its immutable wrapper explicitly passes
  `--continue-after-guard`, so the false marker near the first instrumented step
  will not terminate the replay. The marker itself is excluded, but the
  continued scalar trajectory and fields remain usable. That source revision
  saves a guard field on every subsequent step, producing more data than needed
  but remaining within the 100-GB scratch allocation at the observed roughly
  12-MB compressed field size.
- Commit `abde6e0` adds a validated continuation from the exact atomic step-8000
  checkpoint. It changes no legacy backend, timestep, or coupling parameter,
  uses the corrected warm-up/relative clipping guard, continues rather than
  terminating after a scientifically meaningful capture, and saves the guard
  field only once before returning to ten-step field cadence. A local 384x384
  two-step resume smoke test completed with no capture. A matched smoke using
  active source `147141b` triggered its legacy absolute clipping marker at step
  8001, but all seven checkpoint arrays and all nondiagnostic checkpoint state
  were bitwise/equivalently identical between revisions at step 8002. Evidence:
  `hpc_live/20260911T100125Z-legacy/legacy_fallback_equivalence_step8000_to8002.json`,
  SHA-256 `5c5543ecbc51201cecbfe4e01ffcbb2c6cc97c835e265d319a735c46b2e5204c`.
  The complete suite
  passes 231/231 in 39.25 s; JUnit SHA-256 is
  `7ec04ae5216a61d2a8887238241964f8a0d91672f22aacaacabe942950871c38`.
- Prepared plan `20260911T101424Z-nogit-6bb2f0` is superseded before submission
  because it corrected guard/storage semantics but did not yet add the missing
  closure-specific source/work recorder hooks. It remains `PREPARED` with no
  Slurm job and must never be submitted.
- Immutable source/work continuation plan `20260911T111913Z-nogit-06fc67` is
  prepared but not submitted while two workers remain active. It resumes the
  clean, canonical-matching atomic step-9000 checkpoint and asserts source
  `bebd53b8706301715c216dafe224f2e0c4180aaa`, source-bundle SHA-256
  `d024b52083f7f1c4042d4d8a99987cff8c5c17160499c8d2bbbc2d4d1f8cffb5`,
  clean recovery SHA-256
  `bb7d4bdbd23c336f388ee28872155aea0ef280835bc87588edc8bae0b8c76e20`,
  and immutable HPC input SHA-256
  `55e7d13a136b7876eb198ba382f285e34036847dcbb8008644d1b58ca4374108`.
- The active HPC application manifests say `UNCOMMITTED` because their scripts
  queried Git from the parent stage directory. This does not make their source
  ambiguous: both immutable wrappers assert the detached commit and verify the
  source-bundle hash. `active_hpc_provenance_attestations.json` records those
  independent identities. Commit `a236192` corrects the lookup for queued runs.
  The final staged jobs use commit `a173624`, shared verified bundle SHA-256
  `c0ba9161a8f519436deec75c7017f03614fa04de43328221996b9efd0dd545b3`.
- The staged seed-5101 refinement now halves both the external phase-increment
  target (0.02 to 0.01) and the actual configured PF timestep (0.04 to 0.02).
  Thus the required refinement is exercised even if the external limiter never
  binds. Seven focused runner/movie/analyzer tests pass.
- The factor-two refinement is now immutable prepared plan
  `20260911T110924Z-nogit-6b926f`, with no Slurm job ID and state `PREPARED`.
  Its HPC input archive SHA-256 is
  `dfe68bb909d310bab479d328b3e76f9d1235a6b08d59c499d526f18a78d87583`;
  it asserts source `a173624e4a0b3b3e5074166a2fe0229ed2f34398` and the already
  attested source-bundle SHA-256
  `c0ba9161a8f519436deec75c7017f03614fa04de43328221996b9efd0dd545b3`.
  It is ready for submission when either active worker reaches a verified
  terminal state, and has not been submitted early.
- The paired corrected seeds are also immutable prepared plans with no Slurm
  jobs: seed 5102 is `20260911T112155Z-nogit-6247ef` (HPC input SHA-256
  `8c7596aeff28ebd6a02919368f5a0dab2de86f4404f20922c2c7309273c72a9b`)
  and seed 5103 is `20260911T112156Z-nogit-439d11` (HPC input SHA-256
  `08fdd599d2816544a469839bde5388fdea9ceeff4dc68b71c73a6aec9d78d90c`).
  Each prepares its 1068-grain state and compacts to exactly 800 within the
  immutable job before starting the same attested `a173624` production model.
- The complete active/prepared/superseded queue is machine-readable at
  `hpc_execution_queue_20260911.json`, SHA-256
  `817f517efe52e7bcf471e17c09bea5d8d8a94a068aa7de5c3a544df5a20540fe`.
  It records the two-worker rule and explicitly marks both unused plans as
  never-submit entries.
- `production_revision_equivalence.json` records that corrected seed-5101
  commit `8bb7837` and the staged refinement/paired-seed commit `a173624` have
  byte-identical FFT mechanics, PF solver, kernels, kinematics, and production
  config. Runner/provenance controls and the designed dt/target/seed factors are
  the only execution-relevant differences. The attestation SHA-256 is
  `2dfccd599372601b1266fc7583e764f2e37ca35da4b113c15da28a08431dfed2`.
- The production analyzer now reports timestep agreement at matched physical
  time and matched grain-size progress, evaluates interpolated energy-history
  norms, checks avalanche-class consistency, honors explicit source
  attestations, and creates dedicated timestep, legacy/corrected, and seed
  comparison plots. Its three focused tests pass.
- The complete latest branch suite passes 218/218 tests in 73.43 s with zero
  failures/errors/skips. The durable JUnit is
  `results/validation/qiu_fix_provenance_analysis_tests.xml`, SHA-256
  `c43697f80a640194da43ce6cb40410f334e604ec6afd821a9fce67cb3df8f28b`.
- Movie commit `cc1461e` fixes the final rendering path to merge coarse frames
  with dense transition fields and deduplicate by step, so the legacy onset is
  not reduced to 200-step cadence. Rendering is streamed, uses a symmetric-log
  field scale for the extreme legacy dynamic range, and writes explicit
  pre/transition/post contact-sheet panel metadata. Two focused tests and a
  two-frame 384x384 GIF/contact-sheet smoke render pass; FFmpeg is unavailable,
  so GIF is the retained animation format unless it is installed later.
- Commit `167c30c` adds an explicit, metadata-recorded transition-step override
  for legacy contact sheets. It preserves the raw diagnostic-capture path,
  step, and reasons while centering pre/transition/post panels on the onset
  selected from the completed per-step scientific analysis. Five focused tests
  and a full-size seven-frame renderer smoke test pass.
- A step-500 live render exposed that sparse-cadence ordinary coarsening could
  be mislabeled as a transition by the contact-sheet selector. The selector now
  centers pre/transition/post panels only on the recorded diagnostic-capture
  step; a run without a capture uses neutral initial/midpoint/terminal roles.
  The corrected four-frame full-size smoke render and 10 focused movie,
  analyzer, and runner tests pass.
- The complete suite after that renderer correction passes 229/229 tests in
  65.14 s with zero failures, errors, or skips. Durable JUnit:
  `results/validation/qiu_fix_capture_selection_tests.xml`; SHA-256
  `29a634fcaddb3f2a341a434d803f3e9ea8e25547d72c453b70616c748505eab5`.
- Commit `92fcfc9` adds a reproducible HPC3 source-attestation command. It
  independently rehashes each input archive and source bundle, verifies the
  asserted commit is in the Git bundle, and has generated attestations for both
  active jobs under `hpc_live/source_attestations`.
- A terminal-decision validator now enforces all five full production roles,
  explicit backend identity/classification, full source SHA, poolable terminal
  records, nonempty gate groups, and a zero-failure/error/skip final suite
  before `qualification_decision.json` can be accepted.
- The complete current branch passes 227/227 tests in 61.95 s with zero
  failures/errors/skips. Durable JUnit:
  `results/validation/qiu_fix_current_branch_tests.xml`; SHA-256
  `2564dc2f130a26c4139a40343e06276f810e42db640f2ae4b3f13b3f764336a6`.
- GitHub CLI authentication for `ukaiiaku-maker` is valid and no PR currently
  exists for this branch; the final report/decision commit can therefore be
  handed off as a new unmerged pull request once terminal evidence is complete.
- The accepted fixed-physical-domain reduced spatial study is
  `20260911T013000Z-spatial-convergence-r2` at dx=1, 0.5, and 0.25. At matched
  t=4 every grid gives N=17 and the same physical population grain size. From
  dx=0.5 to 0.25, differences are 0.22% in interfacial/total-energy density,
  2.77% in stress p95, 0.14% in mean compactness, and 1.48% in compactness p95.
  Elastic-energy density improves from a 3.73% coarse-pair difference to 2.59%
  on the finest pair but still fails the requested 2% sensitivity threshold;
  it remains about 4e-6 in absolute density. The convergence trend and narrow
  failed gate are retained honestly. The earlier two-grid study is superseded.
- The long reduced target study at
  `20260911T011000Z-reduced-timestep-convergence-r2` carries both target=0.02
  and target=0.01 from N=18 to N=5 at step 2405/t=96.2. The trajectories are
  bitwise identical, all observable tolerances pass, equilibrium residual stays
  below `4.6e-13`, and no guard or energy rejection occurs. This is explicitly
  *not* counted as exercised target convergence: the minimum external limits
  are 6.91 and 3.45, both far above dt=0.04. Full production refinement remains
  mandatory. Superseded attempt `20260911T010000Z-reduced-timestep-convergence`
  is excluded because the production 10% population guard is inappropriate in
  an 18-grain cell and stopped ordinary two-grain loss at step 167.

## Decisions recorded

- Preserve the legacy implementation and historical output unchanged.
- Do not treat saved frames as restart states.
- Replay deterministically from the immutable initial state.
- Keep superseded diagnostic analyses with explicit exclusion manifests.
- Exclude every short full-size preflight/profile from production inference;
  all five now have root-level exclusion manifests, while r5 remains accepted
  solely as the engineering preflight for the immutable HPC production path.
- Test the whole-grain-area/per-boundary source-overcount hypothesis rather than assuming it is the sole cause.
- If mathematical equivalence to the archived reference is not established, qualify the corrected backend as `FFT_EIGENSTRAIN_V2`, not `QIU_REFERENCE_V2`.
- Treat the non-equilibrated legacy Green projection as a demonstrated defect,
  but do not yet assign sole causal responsibility for the historical
  avalanche; source overcounting and explicit lag remain live hypotheses.
- Select `plane_strain` for the corrected FFT surrogate to preserve the old
  surrogate's constitutive choice. This is not a claim that the archived Qiu
  line kernel is a plane-strain eigenstrain model.

## Next automatic action

Continue monitoring both immutable HPC jobs without adding a third full worker.
Retrieve and checksum each terminal result. If the legacy replay completes
through its endpoint, retain it for exact trajectory/dense-field evidence but
launch prepared source/work continuation `20260911T111913Z-nogit-06fc67` when
that slot is verified free. The factor-two refinement follows after the
source/work continuation or another active worker reaches terminal state.
