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
