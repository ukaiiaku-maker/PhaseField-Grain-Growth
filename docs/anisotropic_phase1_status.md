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

## Phase 2: HPC3 environment repair and manufactured checks

Current scientific source HEAD: `f4aa53f2725c7f9bb3523e1165fe80203d5f487d`.
Job 55932113 failed before science because the reused environment lacked pandas.
Its inputs verified on compute. Two fetches completed; the archive SHA-256 was
independently verified, but its failed application remains unqualified. Actual
CPU 1 s; elapsed 10 s; maximum RSS 209.25 MB; CPU efficiency 5.00%.

Prepared a separate versioned environment overlay on compute, retaining the
old numerical dependencies and adding pinned pandas/pytest/matplotlib. No QIU
or historical environment was modified. Job `55932211`, runner
`20260911T002932Z-nogit-2371f7`, was submitted with the same 2 CPU-h maximum
request. A bounded query confirmed RUNNING on `hpc3-14-03` after a runner status
timeout. That timeout caused no resubmission.

The second allocation runs the complete 190-test expected suite, then initial
reconstruction/dense constitutive checks, then isotropic and A2 sharp-network
normal gradient flow at two timesteps and anisotropic TJ relaxation. These
are mathematical network tests; they cannot satisfy the missing PF gates.

Added validation, scientific and handoff reports that explicitly mark all
unperformed tests and all ten unreleased production trajectories. Current
scientific decision remains UNRESOLVED. Next action: retrieve the known job,
verify its actual outputs and accounting, and diagnose any failed gate.

## Phase 3: corrected validation submission and shared-core snapshot

Environment setup in 55932211 succeeded; test collection then failed because
`PYTHONPATH` omitted the source root needed by the `scripts` namespace. All
scientific calculations were still unstarted. Its verified, twice-fetched
archive is `a0d61c5cf56a4406bd1ddbc6a49a1fac722924523cfafe37aced80db262bd01e`.
Actual CPU 19 s, elapsed 304 s, peak RSS 473.55 MB, CPU efficiency 3.12%.

Corrected wrapper, same scientific source f4aa53f, submitted as **55932614**,
runner `20260911T004050Z-nogit-3e2ea8`. Request: 1 CPU, 6 GiB, 1 h on
SDILLON1 / standard. A local collection-only check found 190 tests in 8.05 s;
it executed no tests. Ledger reconciliation ran twice: updated, then unchanged.

The submission was briefly stopped before upload/sbatch when a new QIU job
appeared. Its account was then verified as SDILLON1_LAB, independent of the
550-SU SDILLON1 balance. Only the 216 CPU-h legacy QIU request draws on our
selected account. Empty Slurm reconciliation preceded resuming the same bundle.
Another mkdir timeout was resolved with the existing runner's bounded 55-second
SSH directory helper; no installed runner or authentication configuration changed.
A query confirmed 55932614 RUNNING on hpc3-14-03. No duplicate was submitted.

QIU later became clean at `847eb06e52c7276bf42ce9ab2f1a73ad74d21127`.
The shared-file differences were reviewed read-only in
`docs/anisotropic_qiu_shared_core_review.md`. Full-field corrections are guarded;
the generic label/previous-state changes are memory optimizations. None was
imported into this branch. A new preproduction review is still required if QIU
changes again; its earlier dirty state is not used as a claim about this snapshot.

Added compact evidence/report assembly and an idempotent local ledger sync tool.
Current decision: unresolved PF implementation; wait for the already running
mathematical qualification result and diagnose its actual tests.

## Phase 4: real-network failure, phase pullback and corrected geometry

Job 55932614 passed all **190 tests in 559.20 s**, then rejected the exact
matched initial network with `ValueError: unclosed periodic boundary`. No
constitutive normalization or manufactured evolution was executed. Its verified
archive is `901bd3d767793a577db1c31769162ee8c9e9764960c0f9e3471d57be3fc1d677`.
Elapsed 600 s, CPU 72 s, CPU efficiency 12%, maximum RSS reported as 2.14 GB.
The failure is retained; it is not a successful qualification.

Commit `4d0e94a` adds the exact triangle-network energy derivative with respect
to phase values, including resolved TJ constraints. Three static tests check
finite differences and phase-sum closure. This is not an evolving PF kernel.
Commit `81c8785` replaces rounded-coordinate equality with actual periodic
minimum-image tolerance matching and adds coordinates to topology errors.
All 15 small mathematical tests passed locally in 1.02 s; no evolution ran
locally. The complete suite and scientific checks run only in HPC3 allocations.

Prepared regression-only run `20260911T005754Z-nogit-669634` was superseded
before submission. It has no Slurm ID and must not be submitted.
Job **55932951**, runner `20260911T010258Z-nogit-ae4ebc`, tests source `81c8785`;
archive SHA-256 `96424071a82b2675bc136fa4efbfd4f4dacbf3c024cd3908845c090ebad63501`.
It requests 1 CPU, 3 GiB, 1 h on SDILLON1 / standard (maximum 1 CPU-hour).
The verified dedicated environment is reused unchanged. Geometry and sharp
network evolution run independently after the full suite, retaining either
failure. Account balance before submission was 549 SU. Remote directory timeout
was reconciled and directories verified through the transfer alias; exactly one
Slurm job was submitted. Latest status at this phase: RUNNING.

The PF integration audit additionally identifies configured/accepted timestep
clock discrepancies and step-based cadence/termination contracts. No historical
or QIU source was edited. Current scientific decision remains
ANISOTROPIC_IMPLEMENTATION_UNRESOLVED. Next action: retrieve and verify 55932951,
then diagnose the actual geometry and dissipation results before further release.

## Phase 5: initial-network and constitutive gates pass; sharp-time refinement

Job 55932951 passed **193 tests in 326.70 s** (480 deprecation warnings).
The periodic vertex-matching correction resolved the real initial-network
closure failure: 43,540 vertices, 44,340 elementary edges, 2,400 TJ-ended pair
paths, zero free loops. All edges belong to the closed periodic graph. A2
edge/path energy partition error is 1.87e-16. Dense constitutive measurements
passed all four A2 gates: p95/p05 gamma 1.64454, M 2.70451, M*S 41.7338,
minimum dense stiffness 0.128296. Normalizations remain provisional pending
geometric convergence. A3 has no release basis.

The job nevertheless failed its independent sharp-network time tests. A2
energy decreased, but dissipation-balance errors at dt=.002/.001 were
4.2868%/2.1004%, above 1%. The dt=.01 TJ converged to residual 9.995e-8;
its half-dt companion exhausted an equal-step budget, which represented only
half as much physical time. These are retained failed tests, not waived gates.
Job accounting: elapsed 866 s, CPU 181 s, efficiency 20.90%, peak RSS 2.14 GB.
Archive SHA-256 verified:
`86e69f6d877bef320b26e82dcf5d7e5f3ca6212e475fc93e49e2a22789ac8168`.
An idempotent second fetch completed. Compact evidence and a provisional
constitutive-range plot are assembled at
`results/long_time_kinetics_900K_anisotropic_20260910/20260911T010258Z-initial-network-time-unresolved/`.

Commit `54189b6fc7767535018afc5112bbffb3f47b8831` refines loop timesteps to
.00025/.000125, retains the 1% balance threshold and adds a 1e-4 final
position/energy refinement threshold. Both TJ timesteps now receive at most
200 physical time units. Constant-sample correlations return undefined rather
than roundoff-induced spurious correlation; the old raw A0 rank value is
retained as superseded evidence. No constitutive parameter changed.

Prepared run `20260911T013122Z-nogit-e7f36e` has source archive SHA-256
`7eb06d15256487ef35df2a3b6cb85babbc60c25a2f9ddf1c50718b62b5a27b54`.
Request: 1 CPU, 3 GiB, 1 h, SDILLON1 / standard. Preflight balance 547 SU,
remote quota 6.42 GiB/1 TiB and 43.75k/8m inodes; local free storage 197 GiB.
Only the two independently owned QIU jobs were active. The first lightweight
remote-directory operation timed out before upload/sbatch; the prepared record
has no job ID, and a name-based Slurm/accounting query found no submission.
Next automatic action: resume this same immutable bundle after directory
reconciliation, retrieve its complete test/time-refinement evidence, and
retain unresolved PF coupling and reduced-polycrystal gates.

The directory helper also timed out at 55 s, but a bounded login-node `stat`
and single-directory `mkdir` subsequently succeeded. Retrying the same
prepared bundle assigned exactly one job: **55933300**. No duplicate run ID or
Slurm submission was created. Ledger sync completed. The scientific source is
54189b6; full tests and corrected sharp-time checks are running in that allocation.
