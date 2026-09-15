# Qiu anisotropy and support-repair status

Updated: 2026-09-15 11:37 PDT

- Branch: `codex/qiu-anisotropy-support-repair-v1`; integration base
  `88fd5dddd3fc425d0b920373cd3698080a3e7452`.
- Model identities: `QIU_LEGACY_FORENSIC`, `FFT_EIGENSTRAIN_V2`,
  `QIU_SI_REFERENCE`, and `QIU_SI_REFERENCE_ANISO` remain distinct.
- Jobs inspected: `55930486`, `55932457`, `55948258`, `55950433`.
- Lineage decision: no exact continuation chain. The legacy second job is an
  overlapping diagnostic replay; the FFT jobs are separate timestep cases.
- Retrieval: terminal FFT job `55932457` was reconciled from a stale runner
  state and retrieved with runner checksum verification.
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
  `20260915T182600Z-nogit-0d1663`, is running the corrected source commit
  `f2537c0fae40fbe114c13d2b4bfaef3f2c34efe3`; its source archive SHA-256 is
  `cca8f546f1e310eaf8f73dfd7f30ac44d4d6f564caa9b4985c05f0beaac4a592`.
  An initial remote-directory setup timed out before upload or `sbatch`; a
  complete Slurm, accounting, command-path, comment, and remote-file check
  found no submission before the successful retry. Prepared
  plan `20260915T181823Z-nogit-edb039` has an empty source archive, was caught
  before submission, and must never be launched. Prepared plan
  `20260915T181845Z-nogit-76f707` is also superseded because it would have
  reported valid scientific gate failures as infrastructure failures.
- Next automatic action: monitor and retrieve support qualification job
  `56040378`, run its committed postprocessor, and enforce Gate 1 from the
  resulting scientific classification.
