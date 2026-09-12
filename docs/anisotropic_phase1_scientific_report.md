# Anisotropic Phase-1 scientific report

Classification: **ANISOTROPIC_PF_TJ_ENERGY_GATE_FAILED**.

The implemented constitutive candidate is the preregistered A2 strong rounded
fourfold law: g_min=0.65, lambda=0.85, support power=16, mobility exponent=2.
It is intended for evaluation, not yet selected for production. The inverse
energy–mobility relation is a synthetic hypothesis. Orientations remain fixed,
and the principal comparison is to use the exact isotropically prepared state
without anisotropic equilibration before time zero.

Independent Cahn–Hoffman line forces and their polygon energy gradients pass
lightweight mathematical tests. This result does not establish anisotropic PF
motion, triple-junction relaxation in PF, or a qualified production model.

HPC job 55932951 reconstructed the exact matched initial state successfully:
43,540 vertices, 44,340 elementary edges and 2,400 pair paths, all ending at
junctions (zero free closed loops). Total polygonal boundary length is
19,434.318749. Periodic graph closure and edge-to-path energy partition passed;
the A2 relative energy partition error is 1.87e-16. This is discrete topology
closure, not continuum geometry convergence.

The provisional initial-network A2 results are:

| Quantity | Minimum | p05 | Median | p95 | Maximum | p95/p05 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma | 0.665237 | 0.739337 | 1.012053 | 1.215869 | 1.236863 | 1.644540 |
| M | 2.417578 | 2.501784 | 3.610909 | 6.766106 | 8.357383 | 2.704513 |
| stiffness | 0.128296 | 0.129883 | 0.198034 | 4.883831 | 8.159629 | 37.601701 |
| M × stiffness | 0.473134 | 0.486438 | 0.701413 | 20.300934 | 68.193146 | 41.733815 |

The dense minimum stiffness across all actual pairs and the fundamental-zone
scan is 0.1282964233, above 0.08. All four constitutive gates pass for A2.
Provisional Cgamma=1.2060485247581734 and CM=0.924620319961451 give means
gamma=1 and M=4. Length-weighted Pearson correlation is -0.9730280192; the
documented weighted rank statistic is -1. A3 is not released. These constants
are not production-frozen because reconstruction convergence remains untested.

The independent sharp polygon test at dt=0.002/0.001 lowered A2 energy from
38.277715 to 36.990254/36.990245 at t=0.2, but its maximum dissipation-balance
errors were 4.2868%/2.1004%, failing the unchanged 1% gate. A2 here uses unit
normalization constants, not the provisional initial-network constants. The
anisotropic TJ reached residual 9.995e-8 at t=180.55 with dt=0.01. The half-dt
run hit a step limit after only half as much physical time and is a failed test,
not evidence of a converged refinement.

The corrected immutable run **55933300** passed the tested sharp-network time
gates. At dt=.00025/.000125, A2 balance errors decrease to 0.5170%/0.2578%
without timestep rejection. Relative final loop position and energy differences
are 1.10e-6 and 2.15e-8. Both TJ relaxations reach residual below 1e-7 by
t=180.55/180.58, and their final positions differ by 3.10e-10. Energy decreases
apart from TJ roundoff of 1.78e-15. This resolves the tested sharp-time failures
while leaving the coupled PF and continuum-geometry gates open.

No anisotropic counterpart of B0, G, T, GT, GTC_GB, GSC_GBTJ_Ks025,
GTSC_GBTJ_Ks025 or GTSC_GBTJ_Ks030_LONG has run. Neither full-domain B0
attribution control has run. All ten remain `PREPARED_NOT_RELEASED` in the
preregistration and ownership ledger. They are not censored or failed scientific
trajectories; there are no anisotropic production trajectories to classify yet.

Consequently growth exponents, kinetics, texture selection, morphology,
mechanism residence, sink partition, event response and the >=10% global effect
requirement are **not evaluated**. Whether anisotropy changes the minimum
progressive-slowing model is **undetermined**. No movies or paired production
figures exist, and no inference about abnormal growth is supported.

The next scientific decision depends on reconstruction convergence and a
validated work-conjugate PF implementation.
The principal risk identified by the audit is substituting a sharp geometric
pressure into a diffuse kernel without establishing the corresponding discrete
energy derivative and pair mobility. Merely completing HPC3 jobs cannot resolve
that risk or qualify the campaign.

## Diffuse operator result

The resumed implementation supplied a discrete anisotropic energy derivative,
complete pair mobility, activation-pressure reuse and accepted-time clock
propagation. Job 55949185 passed 200 tests and all focused diffuse checks except
A2 triple-junction energy descent. Its largest single-step energy increase was
2.6432486007842755.

The controlled job 55949331 reduced accepted dt by factors 1, 2, 4 and 8 at a
common physical horizon. The maximum energy increases were 2.64325, 2.76135,
3.66392 and 3.30394. The failure persisted at the finest dt and did not decrease
with refinement, even though every field stayed finite and the maximum phase-sum
error was 2.22e-16. This is a structural failure of the current diffuse
active-set/projection formulation. It cannot support scientific attribution or
production. The A2 normalization remains provisional, and the >=10% effect
gate remains unmeasured.
