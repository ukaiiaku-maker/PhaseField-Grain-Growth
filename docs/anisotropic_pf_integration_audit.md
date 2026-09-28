# PF integration audit: unresolved contracts

This audit concerns the historical non-QIU path at the audited base. The new
geometry and pullback modules are independent; they do not change these paths.

## Capillary and kinetic operators

The accepted kernel is a nine-point Laplacian, obstacle reaction, local phase
projection and scalar spatial gate. The free-energy diagnostic instead uses
forward two-direction gradients and a continuum obstacle normalization. These
are related discretizations, not an identity proving that an arbitrary geometric
pressure insertion is the discrete energy gradient of the accepted update.

The resolved triangular pullback now gives the exact phase-value derivative of
its polygonal energy, including its TJ endpoints. It is not yet combined with a
verified mobility metric or the historical profile-relaxation operator. A valid
reference correction must subtract the identical isotropic discrete geometric
contribution, account for the phase-value/velocity conversion, and demonstrate
dissipation under refinement. A0 must still reproduce the historical update.

The scalar `solver.mobility_scale` cannot assign three different physical pair
mobilities at one multiphase pixel. A pair operator must remain antisymmetric
before projection and must apply the same mobility to the full permitted driving
force, including shear and event impulses. The current code has no such
anisotropic pair operator.

`_activation_rates` and activation-work records currently obtain scalar
`gb_energy * segment.curvature`; neither is a local anisotropic pressure field.
They must receive the same qualified pressure used by motion, with a documented
mapping from the conservative interface nodes into finite kinetic domains.

## Actual physical time versus configured timestep

Several historical paths consume `config.pf.time_step` directly:

- TJ-coupled mode-flux remaining time;
- serial climb advancement and event-window start time;
- activation and TJ release clocks;
- shear-memory migration;
- compatibility and disconnection event advancement.

The solver can already return a smaller `diag.dt` through its stability bound.
For the historical matched constants, its 0.045 bound is above the requested
0.04, so the original campaign avoids this discrepancy. Strong anisotropy may
reduce the accepted timestep below 0.04. Reusing these clocks would then consume
more kinetic time than physical time and contaminate mechanism comparisons.
An anisotropic runner must advance them using the accepted interval, preserving
initialization behavior and transactional rejection/restart state.

## Cadence and stopping

The historical run loop stops at `max_steps` and schedules diagnostics and
checkpoints by step modulus. For unforced B0 it updates the entity snapshot only
at output steps. Therefore adding a force that reads that cached snapshot would
make B0 dynamics depend on output cadence.

The anisotropic driving geometry must update independently of output. Physical
output targets must not silently change the numerical timestep and therefore the
trajectory. Requested cadence and actual sample times must both be recorded.
Schedule state, normalization state, accepted physical time and all kinetic
clocks must survive restart. The hard endpoint is N<=100 or t=4000; reaching
100,000 smaller steps is not completion.

These are mandatory implementation and test gates, not claims that the original
isotropic campaign used anisotropic timestep control or that its historical
results should be reinterpreted. No mechanism-bearing anisotropic run has been
launched with these unresolved contracts.
