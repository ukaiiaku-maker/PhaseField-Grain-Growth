# Qiu code audit

Audited artifact: Zenodo 15120372 `PF_Codes.zip`, exact local MD5 `6cd49ca72eba89210abb96e700342f12`. Pristine files are under ignored `.external/qiu/PF_Codes/`.

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

The original drivers depend on seed data not included in `PF_Codes.zip` and on legacy `numba`/`sparse` combinations. Immutable regression inputs are therefore the published scripts plus recorded hashes. Reproductions use independently specified small geometries with the audited equations. Curvature-only, idealized shear, and polycrystal shear cases are written to `results/validation/qiu_*`; they are not considered passed until their manifests contain quantitative observables and `validation_passed=true`.

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
