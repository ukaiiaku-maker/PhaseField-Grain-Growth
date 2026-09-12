# Native anisotropic phase-field operator

The anisotropic solver is opt-in through `PFConfig.anisotropy_strength`. A null
strength and `A0_ISOTROPIC` both call the historical Qiu kernel without changing
its arithmetic. This is the exact trajectory-nesting contract used for A0.

For a non-A0 strength, the implemented discrete functional is evaluated over
locally supported grain pairs. At each forward-difference cell,

\[
E_h={4\over w}\sum_c \Delta x^2\sum_{i<j}
\gamma_{ij}(\theta_{ij})
\left[\eta_i\eta_j-{w^2\over \pi^2}
\nabla_h\eta_i\mathbin{\cdot}\nabla_h\eta_j\right].
\]

The pair normal is the direction of `grad(eta_i-eta_j)`. The code differentiates every
term analytically, including the inclination derivative `gamma_theta`; its
returned capillary potential is therefore the exact nodal derivative of the
reported discrete energy on a fixed active-set branch.

The phase exchange at a pixel is

\[
\dot\eta_i=\sum_{j\ne i} M_{ij}(\theta)
\left[{\pi^2\over4w\Delta x^2}(\mu_j-\mu_i)
+{f_i-f_j\over N}\right],
\qquad \dot\eta_j=-\dot\eta_i.
\]

Here `mu` is the derivative of the same `E_h`, and `f` is the existing external
work field. Thus one pair mobility multiplies the complete allowed pair drive.
Each exchange is antisymmetric. The energy contains all diffuse junction
contributions, so no separate Herring or endpoint force is added to PF evolution.

The last evaluated pair chemical-potential difference is also the capillary
pressure source for disconnection activation work. This removes the former
anisotropic inconsistency in which activation used `gb_energy * curvature`
while PF would have used a different thermodynamic derivative.

Pair terms with identically zero local density are skipped without changing the
functional. Force/energy finite differences keep the global grain activity fixed,
while topology and extinction are tested separately. The stability bound uses a
conservative dense upper bound on the selected law's reduced coefficient
`M(gamma+gamma'')`.

HPC3 qualification rejected this implementation. Job 55949185 found a
2.6432486007842755 one-step energy increase during diffuse A2 triple-junction
relaxation. Job 55949331 repeated that calculation at accepted timestep factors
1, 1/2, 1/4 and 1/8; the maximum increases were 2.64325, 2.76135, 3.66392 and
3.30394. Finiteness and phase-sum conservation held, but the energy defect did
not converge away. The current active-set/projection formulation is therefore
not a qualified discrete gradient flow and must not be used for reduced or
production trajectories.

The Gate-1 correction removes both structural inconsistencies found in that
implementation. First, the gradient coefficient is fixed for every pair; it no
longer contains the local phase count that changed the energy of surviving pairs
when a third phase left a junction stencil. Second, the update no longer clips
and renormalizes the completed trial state.

For the obstacle constraint, the solver first computes all pair exchanges at a
node and sums each phase's proposed outgoing flux. If `dt` times a donor's total
outgoing flux exceeds its available phase fraction, all exchanges from that
donor receive the same availability factor in `[0,1]`. Every limited exchange
remains antisymmetric and has the form

\[
q_{ij}=a_{ij}(\mu_j-\mu_i), \qquad a_{ij}\ge 0,
\]

so the effective pair graph is symmetric positive semidefinite and its unforced
capillary work is nonpositive. The only post-update adjustment is a bounded
roundoff correction below `1e-12`; a larger obstacle violation raises an error.
No stress cap, constitutive mobility reduction, energy backtracking or smaller
stability limit is used.

The anisotropic route likewise does not delete a merely small positive grain at
the configured historical extinction threshold. A grain becomes inactive only
after its phase field reaches exact zero through the conservative exchanges.
This avoids a second zero-and-renormalize energy jump. The A0 route retains the
historical extinction arithmetic for exact nesting.

Focused local Gate-1 checks pass. The original TJ horizon is strictly energy
decreasing at both `dt=6.36493e-5` and `dt=3.18247e-5`, with zero positive
increments and maximum phase-sum error `2.22e-16`.

The single consolidated HPC3 qualification accepted the correction. Job
55965236 passed all 203 repository tests and every A0, derivative, planar,
inclusion, TJ, timestep, restart, topology, finiteness, nonnegativity, and
phase-sum gate. The TJ had zero positive energy increments at both accepted
timesteps. The verified result archive is
`412130123ca3e692c13979c06920b2d7be2e54eb591bcecdfc29d896c460fd7c`.
Production remains closed until the preregistered reduced polycrystal pilot
passes.
