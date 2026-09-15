# Qiu anisotropy and support-repair status

Updated: 2026-09-15 11:34 PDT

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
  `20260915T181129Z-nogit-2b882f`, authoritative Slurm job `56040138`, is
  running from exact step 126. A lost submit acknowledgement created orphan
  duplicate job `56040137`; it was detected at six minutes and cancelled.
  Both jobs used identical inputs, and only ledger-owned `56040138` is
  authoritative.
- Support qualification: corrected prepared plan
  `20260915T181845Z-nogit-76f707` contains source archive SHA-256
  `37f5211a5fef614c9fe03fd83bba8fe177599573294688980331ef41ef8057ec`.
  It is intentionally unsubmitted until the closure worker exits. Prepared
  plan `20260915T181823Z-nogit-edb039` has an empty source archive, was caught
  before submission, and must never be launched.
- Next automatic action: retrieve and verify the historical closure, then
  submit the corrected 192×192 support qualification at both timesteps and
  all three KKT tolerances without exceeding two scientific workers.
