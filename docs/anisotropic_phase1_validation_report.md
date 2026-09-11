# Anisotropic Phase-1 validation report

Current classification: **ANISOTROPIC_IMPLEMENTATION_UNRESOLVED**.
This report distinguishes completed checks from missing campaign gates.

## Provenance and isolation

The audited base is `9f66c8d7a5a266687284d8da35aefbc6062808f7`.
Historical production source is `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`.
The initial NPZ SHA-256 is
`106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`;
metadata SHA-256 is
`ff66cb5f938b5b70fe9e611674956f1e45f30f0979714507a853750f61e1a9bb`.
Both match locally and in the first compute allocation. Twenty-four modest
historical non-QIU manifest, performance and energy files match the canonical
integration inventory. Large historical inventories have not been rehashed.

All changes are in `/private/tmp/pfgg-anisotropic-cahn-hoffman-v1` on
`codex/anisotropic-cahn-hoffman-phase1-v1`. The historical solver, kernels,
kinematics, simulation, tracker and QIU code are unchanged on this branch.
QIU has only been read; no QIU job-control action or artifact write was issued.
A byte-level comparison of all 75 pre-existing tracked source files against the
audited base passed. The new modules are standalone and are not imported by the
historical execution path.

## Completed lightweight mathematical checks

Fifteen tests passed locally in 1.02 seconds at scientific source `81c8785`.
These contain no evolving PF or
network trajectory. The tests verify:

- analytic first and second inclination derivatives against fourth-order finite differences;
- rotation objectivity, fourfold symmetry, grain exchange and normal reversal;
- exact A0 energy, derivatives, stiffness and mobility;
- negative polygon-energy gradients at all interior and endpoint coordinates, for three perturbation sizes;
- exact force telescoping and closed-loop force cancellation with nonzero pressure;
- isotropic circle pressure and straight-boundary zero interior force;
- finite TJ energy variation and isotropic 120-degree force closure;
- connected graph partition, disconnected same-pair components and periodic winding;
- proper endpoint forces for a loop returning to a TJ;
- three-phase upper-envelope reconstruction and periodic label exchange;
- static Wulff weighted-curvature convergence across 128, 256 and 512 vertices;
- explicit rejection of unresolved zero-length, branch and half-box geometry.
- exact triangle-energy pullback to phase values in two- and three-phase cells;
- zero phase-sum derivative and exact A0 geometric reference subtraction.

This proves the independent polygon first variation. It does **not** prove a
work-conjugate coupling of those forces to the historical diffuse PF kernel.

## HPC3 execution

| Job | Purpose | Source | Request | Result |
| --- | --- | --- | --- | --- |
| 55932113 | Existing environment verification, regression, initial geometry | db89a77 | 2 CPUs, 12 GiB, 1 h | Failed before science: missing pandas |
| 55932211 | Dedicated environment setup and regression | f4aa53f | 2 CPUs, 12 GiB, 1 h | Environment succeeded; collection failed because wrapper omitted source root from PYTHONPATH |
| 55932614 | Corrected wrapper, full regression and matched initial geometry | f4aa53f | 1 CPU, 6 GiB, 1 h | 190 tests passed in 559.20 s; initial network failed periodic closure |
| 55932951 | Tolerance-based vertex matching, pullback regression and independent sharp-network flows | 81c8785 | 1 CPU, 3 GiB, 1 h | 193 tests passed; initial network and A2 constitutive gates passed; sharp-network time tests failed |
| 55933300 | Refined loop timesteps and equal TJ physical horizons | 54189b6 | 1 CPU, 3 GiB, 1 h | Completed: 193 tests, initial geometry, constitutive and tested sharp-network time gates pass |

All use `SDILLON1` / `standard`, no GPUs and no high QOS. Numerical libraries
use one thread. Requests respect the live 6 GiB/core memory rule; no
parallel-speedup claim is made. The balance before the latest submission was
547 SU; the concurrent legacy QIU maximum request was 216 CPU-hours on SDILLON1.
Its other 720 CPU-hour request uses SDILLON1_LAB under separate QIU ownership.
The five submitted anisotropy jobs total 7 CPU-hours maximum requested exposure.
The superseded regression-v4 plan has no Slurm ID and incurred no allocation.

Job 55932113 elapsed 10 s, used 1 CPU-second, had 5.00% CPU efficiency and
209.25 MB maximum RSS. Its runner retrieval is correctly marked incomplete.
The archive itself independently matches its compute-side SHA-256:
`4e1f446f06dbd0248485a0c0721f23b1f5a04ef1f59c83b2ee6b9b5465b85b73`.
An idempotent second fetch was performed. Failure remains an environment failure.

Job 55932211 elapsed 304 s, used 19 CPU-seconds, had 3.12% efficiency and
473.55 MB maximum RSS. Its independently verified archive is
`a0d61c5cf56a4406bd1ddbc6a49a1fac722924523cfafe37aced80db262bd01e`.
The separate environment overlay passed pinned package verification; only
the wrapper import path failed. The original Conda environment was untouched.

Job 55932614 elapsed 600 s, used 72 CPU-seconds, had 12% efficiency and
reported 2.14 GB maximum RSS. Its independently verified archive is
`901bd3d767793a577db1c31769162ee8c9e9764960c0f9e3471d57be3fc1d677`.
The complete 190-test suite passed with 480 deprecation warnings. Reconstruction
then raised `ValueError: unclosed periodic boundary`; no constitutive ranges or
manufactured evolution were produced. Rounded-coordinate hashing was replaced
by minimum-image tolerance matching for the next immutable source.

Job 55932951 elapsed 866 s, used 181 CPU-seconds, had 20.90% efficiency and
reported 2.14 GB maximum RSS. Its twice-fetched archive independently verifies as
`86e69f6d877bef320b26e82dcf5d7e5f3ca6212e475fc93e49e2a22789ac8168`.
All 193 tests passed in 326.70 s (480 warnings). The matched initial graph closes
with 2,400 TJ-ended pair paths; A2 satisfies all four constitutive strength gates.
See the scientific report for the actual distribution table and provisional
normalizations. Fixed-diagonal reconstruction convergence remains untested.

The A0 sharp circle radius errors are 8.11e-7 and 4.06e-7 at dt=.002/.001;
its dissipation balance errors are below 5e-12. A2 energy decreases at both
timesteps, but maximum balance errors are 4.2868% and 2.1004%, failing 1%.
The dt=.01 anisotropic TJ reaches residual 9.995e-8. Its half-dt test fails
because its equal-step budget is only half the physical horizon. The failed
raw result is retained. Job 55933300 uses dt=.00025/.000125 for the loop and
equal 200-time-unit TJ horizons, with unchanged balance/residual tolerances.
It also tests final loop position and energy differences below 1e-4.

Job 55933300 completed with application/finalization exit 0 and a complete
finalization marker. All 193 tests passed in 1219.61 s after a documented
BeeGFS I/O stall in an existing plotting test. No cancellation or replacement
submission occurred. The final archive was verified and fetched twice:
`de36bf39e1d2fc6862f0a3166cb9187ab3c9c322938bb80b703cc64e103b0923`.
Elapsed 1742 s, CPU 211 s, efficiency 12.11%, reported peak RSS 2.14 GB.

At dt=.00025/.000125 the A2 sharp loop has maximum relative dissipation errors
0.00517020/0.00257837 (0.517%/0.258%), satisfying the unchanged 1% gate.
No timestep rejection occurred. Relative final position and energy differences
are 1.10e-6 and 2.15e-8, below 1e-4. A0 circle radius errors are 1.01e-7 and
5.07e-8. Both anisotropic TJ relaxations reach force residual below 1e-7 at
t=180.55/180.58, before their common t=200 ceiling; final junction positions
differ by 3.10e-10. The largest energy increase is 1.78e-15 (roundoff).
These are sharp polygon checks with unit normalization constants, not PF
convergence or the reduced polycrystal effect-size study.
The loop uses 128 vertices, initial radius 8 and horizon 0.2. The anchored
three-arm junction uses unit junction mobility, anchors at radius 3 and fixed
grain orientations (0.1, 0.4, 0.7) radians.

All five submitted jobs are terminal and retrieved; reconciliation reports no
unretrieved jobs. Complete accounting and immutable identities are in
[the HPC3 execution report](anisotropic_hpc_execution_report.md).

The first mkdir attempt timed out before transfer or submission. Empty Slurm
reconciliation preceded retry of the same prepared bundle. Exactly one Slurm ID
was assigned to that bundle. A later monitoring timeout did not trigger a
resubmission. No full scientific trajectory or continuation has been launched.

## Outstanding scientific gates

| Requested gate family | Evidence and remaining scope |
| --- | --- |
| A: symmetry/objectivity | Passed for the independent A2 constitutive/geometry primitives. |
| B: isotropic nesting | Exact A0 primitive values and geometric reference subtraction pass; all historical source files are identical. No coupled anisotropic PF A0 path exists. |
| C: Cahn–Hoffman geometry | Straight/circle, endpoint closure, finite TJ variation and static Wulff mesh tests pass. Dynamic PF geometries have not run. |
| D: regularity and strength | A2 dense convexity/stiffness and all initial-network constitutive range gates pass provisionally; continuum reconstruction convergence remains open. |
| E: PF coupling | Not implemented: no pair-mobility/full-force operator, activation-pressure mapping or qualified TJ correction. |
| F: topology | Periodic closure, graph partition, disconnected paths and static split-energy checks pass. Dynamic merges/reconnection/history transfer are not qualified. |
| G: numerics | Static Wulff spatial convergence is tested. Sharp-network time refinement is separate from the unrun PF grid, width, cadence, restart, extinction and no-resurrection checks. |
| H: regression | 193 repository tests pass on HPC3; 75 pre-existing source files are byte-identical to the audited base. This does not supply missing anisotropic PF integration tests. |

- Initial-network reconstruction convergence and frozen normalization acceptance.
- Acceptance of measured constitutive ranges after reconstruction convergence.
- Variational diffuse PF coupling, pair-specific full-force mobility and TJ correction.
- Identical capillary pressure in PF evolution and activation-work diagnostics.
- Anisotropic PF dissipation and all timestep, grid, width, cadence and restart checks.
- Reduced 200-grain PF response, attribution controls and >=10% effect-size gate.
- A fresh QIU clean shared-core release snapshot; the earlier review found no physical correction to import.
- Resource profiling, funded production exposure and all ten production trajectories.
- Full historical artifact audit, matched comparisons, analysis figures and movies.

The static Wulff and polygon tests must not be substituted for these missing
PF or campaign tests. A3 is not released while A2 numerical qualification is
unresolved. No selected production law or qualified campaign is claimed.
