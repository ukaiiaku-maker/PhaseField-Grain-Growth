# Production shear-force pathway audit

## Executed continuous migration path

The production intermediate-stiffness campaign uses
`MigrationClosureSimulation` with `migration_closure=gate_only`.  At each PF
step the following path is executed:

1. `MultiphaseFieldSolver.step()` evaluates the constrained pairwise
   double-obstacle PF update.  The diffuse free-energy functional supplies the
   capillary contribution.
2. `MigrationClosureSimulation._update_physics()` measures each persistent GB
   domain's curvature, local normal velocity, and signed displacement from the
   PF fields.
3. `LocalShearMemory.migrate(beta, displacement, dt)` updates
   `s <- s + beta*dx_n` (plus optional deterministic relaxation, which is not
   active in these production configurations).
4. `LocalShearMemory.internal_shear_stress` returns `tau_int=-Ks*s`.
5. `_boundary_resolved_shear()` adds any full-field resolved shear.  The present
   `local_memory` backend has no full-field term, so it returns `tau_int`.
6. `_boundary_force_balance()` forms `p_shear=beta*tau_int` and the separate
   transient event-release pressure.  `_update_physics()` writes only
   `p_shear+p_event` to the antisymmetric per-phase `driving_field`; capillarity
   is not added twice.
7. `MultiphaseFieldSolver` passes that field to `pairwise_obstacle_step()` on
   the next step.  The local mobility field independently applies compatibility
   and climb pinning.

Thus the sharp-interface diagnostic decomposition is

\[
p_{net}^{diag}=p_{cap}^{diag}+p_{chem}^{continuous}+p_{shear}+p_{event},
\qquad p_{cap}^{diag}=\gamma\kappa,
\qquad p_{shear}=\beta\tau_{int}=-\beta K_s s.
\]

`p_cap` is a sharp-interface estimate of capillarity already evaluated by the
diffuse PF functional.  In the corrected area-loss closure,
`p_chem^{continuous}=0`: defect chemistry changes transition-state work,
pinning state, and GB/TJ sink competition, but is not inserted as a continuous
normal pressure.  `p_event` is the finite release impulse associated with
`normal_release_remaining`; it is not relabeled as chemical pressure.

The diagnostic arrest ratio is

\[
\chi_s=-\frac{p_{shear}}{p_{cap}^{diag}+p_{chem}^{continuous}}.
\]

It is positive when shear opposes the local non-shear drive, approaches one at
mechanical force balance, and exceeds one after reversal.  It is undefined
when the non-shear denominator vanishes.  The event impulse is reported but
excluded from `chi_s`, because it is a stochastic release mechanism rather
than stored backstress.

## Transition-state path

Continuous backstress is distinct from the activation calculation.
`MigrationClosureSimulation._activation_rates()` constructs candidate
disconnection modes and evaluates

\[
W_m^\ddagger
=p_{cap}V_{n,m}^\ddagger
+\tau_{resolved,m}V_{\tau,m}^\ddagger
+\Delta\mu_v N_{v,m}^\ddagger
-\Delta E_{TJ,m},
\]

then uses `DeltaG_eff=max(0,DeltaG0-W)` in the Arrhenius rate.  The mode's
resolved shear includes the projection of `tau_int` on its Burgers direction;
therefore shear can change individual barriers and the selected mode even when
the scalar continuous backstress is similar.  Compatibility filtering can
remove easy modes while a domain is blocked, and the GB/TJ gates separately
change mobility and state residence.

The production mechanism audit must therefore distinguish four coupled routes:

- continuous cancellation of capillary migration by `beta*tau_int`;
- transition-state work `tau_resolved*V_tau^ddagger`;
- mode removal/selection under compatibility gating;
- longer G/T/C residence and altered GB/TJ sink competition.

No discontinuous branch in the source is keyed directly on `Ks`.  Stiffness
enters through `tau_int=-Ks*s`; any arrest transition must emerge from the
coupled force, event, and compatibility dynamics rather than a hard-coded
stiffness gate.
