# Local compact-support Gate-1 report

Status: **COMPLETE — `COMPACT_SUPPORT_OPERATOR_KKT_FAILURE`**

The user authorized the reduced compact-support qualification to run locally
while the HPC3 two-worker ceiling prevented the prepared continuation from
being submitted. Native Qiu, Qiu anisotropic controls, the FFT timestep
refinement, and full Phase-1 production trajectories remain HPC3-only.

## Identity and reconciliation

- Local run: `local-compact-support-step1024-20260915T200027Z`
- Execution site: `local_macos`
- Logical trajectory: `compact-support-qualification-v1`
- Predecessor Slurm job: `56040470`
- Prepared HPC3 continuation: `20260915T185600Z-nogit-8a5672`
- Prepared state: `HELD_NOT_SUBMITTED_LOCAL_FALLBACK`
- Scientific source commit: `db1770b2fe098274ccd8ce3a611ea79932ebd3bd`
- Source archive SHA-256: `7ce21b3a3ece25fa4f09e5efd4a49fe87fa8fdd988a2e42cbabaef7667f6ab85`
- Start checkpoint SHA-256: `4b7138b5d73d30c7cff8d6e9a3d1b5c1891b3202c4214ce86efdd3df99df580e`
- Resolved configuration SHA-256: `8ffcd7f51d42ebc977808b295d1346862adf128c9b64275a31d59c26ba56784d`

The pre-launch `squeue`, `sacct`, and remote-directory reconciliation found no
Slurm identity or remote input for the prepared continuation. Only protected
job `55950433` was live. The predecessor was `CANCELLED`; its exact checkpoint
was retained and loaded at accepted step 1024 and physical time
0.06517691816792172. All 204 executable/configuration files in the source
archive were byte-identical to the recorded scientific commit.

## Local environment and execution

The run uses macOS 26.1 build 25B78 on an x86-64 Intel Core i9-9880H, with
Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3, Numba 0.61.0, and OpenBLAS 0.3.21.
The numerical package versions match the qualified HPC3 environment. The
supervisor runs under `caffeinate` at nice level 10 and constrains OpenMP,
OpenBLAS, MKL, NumExpr, Accelerate, and Numba to one thread.

The campaign completed at 2026-09-16T02:01:42Z after all six frozen cases
reached their preregistered endpoints. No local worker or tmux session remains.
The scientific entrypoint was `scripts/qualify_compact_support.py`; no physical
parameter, tolerance, timestep, stop criterion, or model implementation was
changed.

## Checkpoint and restart evidence

The immutable start state was loaded through the frozen solver with finite
fields and phase-sum error `2.22e-16`. A planned signal produced an atomic
step-1099 checkpoint (SHA-256 prefix `260b5f3bf688`) and a new process resumed
from it. A later controlled supervisor restart produced and resumed exact
step-1179 checkpoint (SHA-256 prefix `03f33cd2abbf`). The run is therefore
durable across both scientific-worker and supervisor process termination.

Every case reached the A0-based physical-time ceiling. Coarse cases stopped at
step 3968 and time 0.25256055790072346. Fine cases stopped at step 7872 or
7936 and time 0.2505237792079889 or 0.2525605579007369. Across all six cases,
maximum KKT residual was `4.44e-16`, maximum/p95 exact support was `7/3`, and
every per-case energy, KKT, support, topology, and scaling target passed.

Python 3.13 rejected the NumPy Boolean produced while serializing that case's
otherwise complete summary. The exact terminal checkpoint was already atomic.
The summary was reconstructed from its stored field and history without taking
another numerical step. A recorded compatibility shim now converts NumPy
scalar diagnostics to native JSON scalars; it does not change the solver,
state, parameters, timestep, tolerance, or stopping rule.

## Independent case lineage

The frozen entrypoint defines six sibling cases. The coarse-dt/KKT-1e-8 case
is the sole case whose existing checkpoint matches and resumes the exact
step-1024 configuration. The remaining five cases start independently from
the shared deterministic step-zero Voronoi field. Their actual checkpoint
parents are recorded as `deterministic_initial_state`, step and time zero;
they are not concatenated to the terminal state of an earlier sibling.

Each checkpoint record contains the campaign run ID, case ID, actual parent
identity, dt, KKT tolerance, case-local step/time, and continuous or replay
role. The principal KKT-1e-10 coarse and fine cases passed exact restart and
restart-time equivalence. The campaign produced 282 checkpoint records. The
terminal fine/KKT-1e-12 checkpoint is step 7872, time 0.2505237792079889,
SHA-256 `8d307aa05840a32ff46801e8f1f38eb623faa0a658d53e1ade0bc4d678172095`.

## Gate-1 classification

The preregistered cross-tolerance bitwise equality criterion failed. At coarse
dt, KKT tolerances `1e-10` and `1e-12` produced the same final field, while
`1e-8` differed. At fine dt, `1e-8` and `1e-12` matched, while `1e-10`
differed; that case also reached the next 64-step boundary. Thus the upstream
classification is `COMPACT_SUPPORT_OPERATOR_TOLERANCE_NONCONVERGED`. Under the
Gate-1 packaging policy, tolerance nonconvergence maps to the terminal
classification `COMPACT_SUPPORT_OPERATOR_KKT_FAILURE`.

This is a scientific gate failure despite the very small measured KKT
residuals and the fact that all individual release targets passed. The local
and HPC3 platform confirmations are not released. Native Qiu and anisotropic
Qiu controls remain blocked.

The machine-readable package is
`results/local_compact_support_gate1_20260915T200027Z`. It contains energy,
KKT, support, timestep, runtime, and morphology tables; checkpoint lineage and
checksums; and `gate1_histories` and `morphology_histories` figures in PNG and
PDF form.
