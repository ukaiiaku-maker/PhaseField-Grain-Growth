# Anisotropic Phase-1 status

Branch: `codex/anisotropic-cahn-hoffman-phase1-v1`.
Initial HEAD: `9f66c8d7a5a266687284d8da35aefbc6062808f7` (clean audited parent).
Worktree: `/private/tmp/pfgg-anisotropic-cahn-hoffman-v1`.
Dedicated state: `/Users/sdillon/HPC3/anisotropic-phase1-state`.
A new VS Code window was requested with `code --new-window`.

## Phase 0

The primary checkout was clean at the audited parent; its branch was
`exp/long-time-kinetics-900K-20260824`. The ten-commit log and worktree list
were inspected before writing. Historical production source exists at
`4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`.
Both initial-state hashes match the brief. All eight non-QIU manifest,
performance and energy hashes match the integration inventory (24 artifacts).
Large result checksums remain a compute-side task. No initial-state evolution
or scientific simulation has run locally.

`matched_isotropic_initial_state = true`. The ten authorized full trajectories
are eight matched regimes plus E-only and M-only B0. A2 is intended, not yet
selected. All normalization and constitutive ladder calculations belong on
HPC3. No environment installation has run locally or on a login node.

## Audit findings and decisions

The existing arclength tracker has nearest-unvisited fallback jumps, pixel-count
lengths, heuristic TJ proximity, and index-based history reuse. It does not
provide the geometric guarantees needed for variational capillarity. Do not
feed its points directly into an anisotropic divergence or freeze normalization
from those approximate lengths.

The PF kernel applies one scalar mobility field to all local pairs. Its rate
contains phase-weighted Laplacians and an obstacle term, plus projected external
forcing. A sharp-interface pressure substitution alone is not proof that its
discrete update equals the requested target. PF and TJ coupling remain gated.

QIU was read only at `9edbd0f`, with uncommitted diagnostics and qualification
files. Its durable status does not freeze the shared core. No QIU source, state,
job, or output was edited. Production gate is CLOSED.

Live cluster: SDILLON1 available balance 552 SU; QIU job 55930486 has 3 allocated
CPUs and 72 h limit (216 CPU-h maximum exposure). Standard is available with a
14-day limit. No reservations were reported. No production resource request is
selected before profiling. No alternate account is authorized.

Tests: pending. Simulations: none submitted.
Current classification: ANISOTROPIC_IMPLEMENTATION_UNRESOLVED.
Next automatic action: implement and verify independent energy/force primitives;
prepare bounded HPC3 mathematical qualification; keep production fail-closed.

## Phase 1: primitive and initial-network qualification

Primitive commit: `78d8560`; initial-network qualification source: `db89a77`.
New files: anisotropy law, exact polygon variation, strict edge-graph ordering,
conforming upper-envelope reconstruction, geometry tests and compute entrypoint.
Local checks: 10 tests passed initially; 12 tests passed in 0.85 s after adding
a static Wulff mesh-convergence check and a returning-TJ-loop regression.
These are small non-evolving unit/analytical checks.

HPC3 environment-verification/full-regression/initial-network job: `55932113`.
Runner ID: `20260911T002118Z-nogit-d279e0`. Source archive SHA-256:
`a99aebe9f1d59d8a1b11fe194aa51645bf729d89cb7ffcf9305399f8df4e25e3`.
Requested: standard / SDILLON1, 2 CPUs, 12 GiB, 1 h; maximum exposure 2 CPU-h.
The second CPU provides memory under the live 6 GiB/core limit; numerical
threads are fixed to one and this is not a parallel-speedup claim.
The first attempt timed out in mkdir before upload or sbatch. Slurm then showed
no matching job, remote directories existed, and the SAME prepared bundle was
submitted once. No duplicate trajectory was launched.

The job uses a read-only existing environment, verified inside the allocation.
It runs the complete existing suite and the first 10 new tests, then measures
A0/A1/A2 on the exact initial network. The later two local tests and returning-TJ
fix require the final-source HPC3 regression before release.

The preregistration JSON stores exact historical configs for all ten planned
trajectories. It is explicitly non-executable: final normalization, full
resource requests and resolved anisotropic configs are not invented.

QIU now has uncommitted changes in `pf/solver.py` and `simulation.py` as well
as its mechanical coupling and diagnostics. No cherry-pick or freeze is valid
yet. Production remains closed. Next action: inspect retrieved qualification,
resolve any geometry failures, then verify the updated source on HPC3.
