# Qiu anisotropy and support-repair status

Updated: 2026-09-16 21:27 PDT

- KKT audit v2 supersedes only the old Gate-1 decision. At exact common time
  `0.2505237792079889` (coarse step 3936, fine step 7872), all fine-timestep
  tolerance states are bitwise identical and the coarse `1e-10`/`1e-12`
  states are bitwise identical. Two independent 128-step `1e-10` replays are
  also bitwise identical and reproduce the historical target. The old mismatch
  combined unequal endpoint steps and unequal parent/platform lineage.
- Revised audit classification: `COMPACT_SUPPORT_GATE_PACKAGING_MISALIGNMENT`.
  Superseding gate classification:
  `COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED_V2`. Production tolerance is
  `1e-10`. The intentional local/HPC platform pair v2 is the next gate; native
  and anisotropic Qiu tasks remain blocked until it passes.

- Branch: `codex/qiu-anisotropy-support-repair-v1`; integration base
  `88fd5dddd3fc425d0b920373cd3698080a3e7452`.
- Model identities: `QIU_LEGACY_FORENSIC`, `FFT_EIGENSTRAIN_V2`,
  `QIU_SI_REFERENCE`, and `QIU_SI_REFERENCE_ANISO` remain distinct.
- Jobs inspected: `55930486`, `55932457`, `55948258`, `55950433`.
- Lineage decision: no exact continuation chain. The legacy second job is an
  overlapping diagnostic replay; the FFT jobs are separate timestep cases.
- Retrieval: terminal FFT job `55932457` was reconciled from a stale runner
  state and retrieved with runner checksum verification. Its complete
  post-simulation package classifies the trajectory as
  `FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE`: all recorded total-energy
  increments are negative, full-horizon radius-squared linear-fit R-squared is
  0.9963, the maximum late one-step population loss is two, and no disconnected
  grains occur. The dt/2 refinement job `55950433` is also complete and was
  retrieved twice with an identical seven-entry verified checksum set. It
  independently has the same no-late-avalanche classification. The paired
  comparison is complete: target time differs by 1.75%, fitted radius-squared
  slope by 1.56%, and shared-interval normalized RMSE is 0.135% for grain
  count and 0.103% for total energy.
- Native plan: `20260912T180419Z-nogit-a25f0c` is prepared and unsubmitted;
  no duplicate exists. The immutable driver does not emit phase-field or
  restart checkpoints, so numeric A0 promotion will use a separately named,
  output-only instrumented native control built from the checksummed source.
- Support decision: the exact-stencil candidate graph, symmetric pair-mobility
  descent, simplex obstacle solve, KKT audit, and strict energy backtracking
  are implemented. The complete local suite passes 225/225. No HPC
  qualification has yet been released.
- Qiu anisotropy decision: the Qiu SI port through commit `97e6a9d` was
  selectively integrated; reduced controls remain blocked by Gate 1.
- Historical restart closure: immutable run
  `20260915T181129Z-nogit-2b882f`, authoritative Slurm job `56040138`,
  completed from exact step 126 through step 128. Its restarted final field is
  bitwise identical to the continuous final field. The verified result archive
  SHA-256 is `f0b4d488e297da133a51acc81321ebf13a74433f1f9f0f9d38f877943a5c2e12`.
  A lost submit acknowledgement created orphan
  duplicate job `56040137`; it was detected at six minutes and cancelled.
  Both jobs used identical inputs, and only ledger-owned `56040138` is
  authoritative.
- Support qualification: Slurm job `56040378`, run
  `20260915T182600Z-nogit-0d1663`, was cancelled by the owner after 89 seconds
  when a deeper source audit found that energy assembly still scanned all
  global phases at every cell. It is ineligible for scientific classification.
  The repaired implementation now evaluates the discrete energy, derivative,
  pair descent, projection, and line-search trials on one frozen stencil-local
  candidate graph; the complete suite passes 225/225.
  The cancelled source archive SHA-256 was
  `cca8f546f1e310eaf8f73dfd7f30ac44d4d6f564caa9b4985c05f0beaac4a592`.
  Replacement run `20260915T184102Z-nogit-583359`, Slurm job `56040470`, ran
  commit `db1770b2fe098274ccd8ce3a611ea79932ebd3bd`; its source archive SHA-256
  is `7ce21b3a3ece25fa4f09e5efd4a49fe87fa8fdd988a2e42cbabaef7667f6ab85`.
  It was owner-cancelled after 12:19 when external full worker `56040506`
  appeared and raised the account above the mission's two-worker ceiling.
  The retrieved partial archive contains an exact step-1024 checkpoint at
  time 0.06517691816792172. At that point mean/p95/maximum exact support was
  1.553/3/6, candidate maximum was 7, KKT residual was 2.22e-16, and the full
  timestep was accepted with decreasing energy. Prepared continuation
  `20260915T185600Z-nogit-8a5672` will resume this checkpoint without replay.
  An initial remote-directory setup timed out before upload or `sbatch`; a
  complete Slurm, accounting, command-path, comment, and remote-file check
  found no submission before the successful retry. Prepared
  plan `20260915T181823Z-nogit-edb039` has an empty source archive, was caught
  before submission, and must never be launched. Prepared plan
  `20260915T181845Z-nogit-76f707` is also superseded because it would have
  reported valid scientific gate failures as infrastructure failures.
- Local Gate-1 fallback: the user explicitly authorized local execution of the
  compact-support continuation while the HPC3 worker ceiling remained closed.
  The pre-launch reconciliation at 2026-09-15 12:45 PDT found only protected
  HPC3 job `55950433` live. Accounting showed predecessor `56040470` cancelled,
  and found no job for prepared run/comment
  `20260915T185600Z-nogit-8a5672` or its expected job name. Its remote directory
  contained no files. The prepared record is now
  `HELD_NOT_SUBMITTED_LOCAL_FALLBACK`, with null job ID and submission time.
- Local run `local-compact-support-step1024-20260915T200027Z` completed all six
  frozen cases; no worker or tmux session remains. It ran under a durable
  tmux/caffeinate supervisor, PID `17159`, at low priority with one
  numerical thread. It executes the byte-identical source archive from commit
  `db1770b2fe098274ccd8ce3a611ea79932ebd3bd` (SHA-256
  `7ce21b3a3ece25fa4f09e5efd4a49fe87fa8fdd988a2e42cbabaef7667f6ab85`).
  The immutable start checkpoint is exact step 1024, time
  0.06517691816792172, SHA-256
  `4b7138b5d73d30c7cff8d6e9a3d1b5c1891b3202c4214ce86efdd3df99df580e`.
  Controlled process interruption/restart advanced exactly through step 1099
  and step 1179 checkpoints without replay from step zero. The resumed
  coarse-dt/KKT-1e-8 case then reached its preregistered horizon at step 3968,
  time 0.25256055790072346. Python 3.13 rejected a NumPy Boolean during summary
  serialization; the terminal field remained intact. Its summary was recovered
  from that exact checkpoint without advancing the solver, and a JSON-only
  NumPy-scalar compatibility shim was installed. All first-case numerical and
  scaling targets pass. All remaining siblings also reached their registered
  endpoints.
- Case lineage is explicit. Only `dt-6.364933414836189e-05_kkt-1e-08` resumes
  the step-1024 checkpoint. The other five frozen-script cases independently
  start from the common deterministic step-zero field. Labeling those five as
  step-1024 restarts would contradict the executable and their checkpoint
  metadata. Every preserved checkpoint now records its actual parent, dt, KKT
  tolerance, case-local step/time, and continuous or replay role.
- The matched-time version-2 audit supersedes only that classification. At
  `t*=0.2505237792079889` (coarse step 3936, fine step 7872), all three fine
  fields are bitwise identical and the `1e-10`/`1e-12` coarse fields are
  bitwise identical. The loose coarse lineage differs by only
  `5.0417530025313027e-17` RMS with no discrete or morphology difference.
  Independent same-tolerance replays are bitwise identical. The stop was
  endpoint/package and lineage misalignment, classified
  `COMPACT_SUPPORT_GATE_PACKAGING_MISALIGNMENT`; Gate 1 is superseded by
  `COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED_V2`. Production KKT is `1e-10`.
- Platform confirmation used the same exact step-7744 parent and ran 128 steps
  locally and on HPC3. Job `56099749`, run
  `20260917T020427Z-nogit-6f9ac5`, was retrieved twice with stable verified
  checksums. macOS and Linux fields differ by `7.809703850438097e-18` RMS and
  `1.9984014443252818e-15` maximum; active sets, candidate graphs, support
  events, timesteps, morphology, and labels agree. Classification is
  `PLATFORM_NUMERICALLY_EQUIVALENT`.
- The pristine archived native Qiu path is independent of the generic
  compact-support operator. A separate maximum-1,000-step preflight is prepared
  as run `20260917T021402Z-nogit-c1b280`. The full native plan remains blocked
  behind a verified `QIU_SI_NATIVE_PREFLIGHT_PASSED` result. The preflight is
  replacement Slurm job `56100055` (run
  `20260917T022751Z-nogit-0eb60c`). The first environment-only attempt,
  `56099932`, failed before a scientific step because the base environment
  lacked `matplotlib`; its diagnostics were retrieved. The replacement reached
  pristine native `find_gb` but failed Numba compilation at
  `functions_4ref_new.py:118` before the first accepted step. Its historical
  classification is superseded by
  `QIU_SI_NATIVE_PREFLIGHT_OPERATIONAL_COMPATIBILITY_FAILURE` because zero
  accepted steps occurred.
- The bounded compatibility audit tested exactly two historical environments.
  Both `anaconda/2022.05` and `anaconda/2021.11` run the untouched 10-step
  native smoke, but only through deprecated object-mode fallback in `find_gb`.
  The minimal homogeneous-array patch has function-module SHA-256
  `1d656039b8dca20bac8f056ad195fdc201775331e73e0f6f5da1e866622f1f78`.
  It matches an independent oracle on five synthetic fixtures and 64 actual
  initialization operations. Pristine and patched downstream GB order,
  references, beta, line density, stresses, and pair terms are exactly equal
  in the historical environment. Modern Numba proceeds past the repaired
  assembly and then rejects unchanged boundary-sorter constructs. Extending
  the patch would violate the mission's array-assembly-only constraint.
- Bounded V2 compiler-gate preflight run
  `20260917T041831Z-nogit-8512a3`, Slurm job `56102284`, tests the pristine
  historical path and is explicitly ineligible for baseline use. It completed
  at exit `0:0` and was retrieved twice with stable runner checksums. Finite
  exact-state checkpoints were recorded at accepted steps 0 and 1, with
  SHA-256 values `eb60b24a2f120a079d28755885a20a1316cf8e6cad9cf6c21d14a2053063c882`
  and `278381ff6b66cdd95f1dca7e77d2d7d64ec502547aa0d3b67fdc52d8c3d9ce7a`.
  `find_gb` and `sort_gb` have no nopython signatures, so the mandatory
  no-object-mode gate fails. Classification is
  `QIU_SI_NATIVE_ENVIRONMENT_INCOMPATIBLE`. The result archive SHA-256 is
  `1a3c03d4f1905579d273207c3f94dce300b4b07552c7c29083a59b693ef8f413`.
  The full native baseline plan `20260912T180419Z-nogit-a25f0c` remains
  unsubmitted; A0 correspondence and all anisotropic production remain
  blocked.
