# Qiu anisotropy and support-repair status

Updated: 2026-09-15 11:22 PDT

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
  no duplicate exists.
- Support decision: the exact-stencil candidate graph, symmetric pair-mobility
  descent, simplex obstacle solve, KKT audit, and strict energy backtracking
  are implemented. The complete local suite passes 225/225. No HPC
  qualification has yet been released.
- Qiu anisotropy decision: the Qiu SI port through commit `97e6a9d` was
  selectively integrated; reduced controls remain blocked by Gate 1.
- Next automatic action: close the missing two-step historical restart and
  prepare the 192×192 HPC3 support qualification at both timesteps and all
  three KKT tolerances.
