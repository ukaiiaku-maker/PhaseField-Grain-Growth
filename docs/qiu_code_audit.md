# Qiu code audit

Audited artifact: Zenodo 15120372 `PF_Codes.zip`, exact local MD5 `6cd49ca72eba89210abb96e700342f12`. Pristine files are under ignored `.external/qiu/PF_Codes/`.

> **2026-09-10 qualification correction.** The archived implementation below
> is a current-geometry line-disconnection construction. The repository's
> historical `qiu_full_field` backend is instead an accumulated point-
> eigenstrain FFT surrogate. The earlier statement in item 6 that its stored
> array was the physical self-stress conflated stress with the negative
> eigenstrain-energy derivative. A formal equilibrium audit now rejects that
> legacy FFT projection (modewise residual about 0.69 in the fixed regression).
> Its old regression results remain historical evidence only. The separately
> named `FFTEigenstrainV2` uses physical stress, for which
> `delta E/delta epsilon* = -sigma`, and is not called a Qiu reference port.

1. **PF terms.** `update_PF` in each `functions_*ref.py` implements a Steinbach-style pairwise multiphase evolution. The interfacial driving combines a sinusoidal double-obstacle-like term with finite-difference Laplacians. The interface inclination is decomposed between neighboring crystallographic reference directions.
2. **Mobility.** Pair mobility is initialized uniformly (`pmobi`), with zero diagonal; the scripts scale it as a function of interface thickness and pass it into the explicit update.
3. **Shear coupling.** `beta(ORR1, ORR2)` returns two reference coupling factors from grain orientation/misorientation. `cal_inc` resolves a local inclination between adjacent references.
4. **References.** The polycrystal uses two reference directions (`ref_theta_i=pi/2`); idealized `[100]` and `[111]` cases use four and six directions (`pi/ref`). Adjacent reference components reconstruct the local tangent.
5. **Internal stress.** `stress_field_bulk_single` precomputes unit line-disconnection kernels. `stress_field_line` integrates boundary tangent increments weighted by the reference coupling factors and superposes shifted kernels. `stress_field_extend` propagates line values across the diffuse interface.
6. **Elastic driving.** `update_PF` resolves local stress onto reference planes, combines it with coupling factors, and adds the resulting signed elastic work to pairwise PF driving.
   The implementation sets `E_el = -E_elastic(...)`; the independent FFT eigenstrain closure must therefore expose the negative elastic-energy derivative as its physical self-stress. A positive source self-work is an unstable sign error, not Qiu-type feedback.
7. **Burgers character.** Burgers content is represented through scalar coupling factors times discretized boundary displacement on reference directions, not event-resolved vector modes.
8. **Density.** Disconnection density follows the discretized GB line/tangent construction; there is no independent nucleation population or first-passage clock.
9. **Boundary conditions.** Periodic wrap helpers are used for geometry and finite differences.
10. **Orientation.** Each grain receives a fixed scalar orientation loaded from seed files (polycrystal) or prescribed by the benchmark script.
11. **Orientation evolution.** Grain orientation is not dynamically evolved; apparent rotation is therefore not an output of the reference implementation.
12. **Integrator.** Explicit forward time stepping followed by `renorm`; despite some imported libraries, the driver does not use an adaptive integrator.
13. **Typical numerics.** Polycrystal: `250x250`, `dx=dy=10`, `dt=10`, interface width `eta=5 dx`, 1001 order parameters, up to 300,000 steps. Idealized 4-reference: `dx=dy=1`, `dt=0.1`, `eta=5 dx`.
14. **Bottlenecks.** Dense `(grain,x,y)` storage, repeated overlap discovery, nearest-neighbor GB sorting, per-segment real-space kernel superposition, stress-array allocation, and Python-level I/O dominate.
15. **2-D assumptions.** GBs are lines, disconnections are point/line-section sources with implicit out-of-plane thickness, and plane stress components `sigma11,sigma12,sigma22` are used.
16. **Missing event kinetics.** Coupling is continuous; there are no discrete `(h,b,Nv)` events, Arrhenius hazards, mode competition, persistent identities, Burgers residual at TJs, climb state machine, or atomic-to-PF displacement ledger.

## Regression policy

The archived cases do not share one readiness state. The `[100]` four-reference
driver is self-contained and has now been staged byte-for-byte for a native
500x500, 200,000-step reproduction. The `[111]` drivers import absent module
names (`functions_6ref_new2_1` and `functions_6ref_new2_2`); the continuation
also requires an absent `OP_t142500.npz`. The polycrystal driver requires the
absent published `PolycrystalSeeds_1000Grains_25_250_1.txt`. Its included
generator is unseeded and cannot reconstruct that realization exactly. These
are archive-level reproducibility limitations, not defects in the legacy
surrogate or in `FFT_EIGENSTRAIN_V2`.

The pristine archive and every staged source remain immutable and checksummed.
The exact audit and promotion gates are in
`docs/qiu_si_reference_reproduction_status.md`. Earlier independently
specified small-geometry regressions remain formula/implementation tests only;
they are not a completed `QIU_SI_REFERENCE` reproduction.

The 2026-09-10 identity tests extract the undecorated `beta` function directly
from the pristine archived abstract syntax tree and compare eight synthetic
orientation pairs, including both piecewise breakpoints and periodic endpoints,
against `two_reference_coupling_factors`. Separate tests cover reference-sector
selection, tensor-axis conversion, stress sign reversal, and grain-order
reversal. These tests validate the audited scalar/reference formulas; they do
not turn the FFT eigenstrain surrogate into the archived line-source model.

The independently implemented regression at commit
`e6b0d8ea52a3d025d49c876fc60359f542c29024` passes. Its matched geometries
are equilibrated before physical time. In the 24-grain polycrystal, the
full-field shear feedback lowers the all-sample velocity-curvature correlation
from 0.436 to 0.388 and increases reverse-curvature motion among the upper
quartile of resolved curvature and speed from 12.7% to 30.8%. The four-grain
case independently increases active reverse motion from 3.09% to 11.38%.
Finite eigenstrain, nonlocal stress, feedback, and continuous phase-field
divergence are recorded in `results/validation/qiu_regression_benchmarks.json`;
the immutable dense runs are listed there. Two rejected proxy/metric attempts
remain beside it as failure records.

## Qualification consequence

The identity gate selects Path B from the qualification directive. The archived
current-geometry line construction remains reference evidence and is not
impersonated by a point/eigenstrain model. New production calculations use the
explicit regime/backend name `FFT_EIGENSTRAIN_V2` and
`configs/production/fft_eigenstrain_v2_qualification_900K.yaml`. The legacy
`QiuFullFieldLegacy` code remains reproducible for historical replay only. This
naming decision is final even if the corrected surrogate eventually exhibits
an elastic transition of its own.
