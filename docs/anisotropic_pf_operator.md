# Native anisotropic phase-field operator

The anisotropic solver is opt-in through `PFConfig.anisotropy_strength`. A null
strength and `A0_ISOTROPIC` both call the historical isotropic Phase-1 kernel
without changing its arithmetic. This is the exact trajectory-nesting contract
used for Phase-1 A0; it is not the Qiu SI four-reference nesting test.

For a non-A0 strength, the corrected discrete functional is evaluated over
locally supported grain pairs. Let

\[
\widehat\gamma_{ij}(p)=|p|\gamma_{ij}(p/|p|)
\]

be the positively one-homogeneous pair norm, evaluated directly from rotated
even p-norms rather than from an angle. At each forward-difference cell,

\[
E_h={4\over w}\sum_c \Delta x^2\sum_{i<j}
\left[\gamma_0\eta_i\eta_j+{w^2\over4\pi^2\gamma_0}
\left(\widehat\gamma_{ij}(\nabla_h\eta_i-\nabla_h\eta_j)^2
-\widehat\gamma_{ij}(\nabla_h\eta_i+\nabla_h\eta_j)^2\right)\right].
\]

The difference/sum polarization preserves exact pair extinction: if either
phase and its local gradient vanish, that pair contributes zero energy. It also
recovers the former isotropic term exactly because
`widehat_gamma(p) = gamma0 |p|` in that limit. The code differentiates the
squared norm analytically. Its product with the Cahn--Hoffman vector is zero at
zero pair gradient, so no undefined normal or inverse gradient magnitude enters
the energy derivative.

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

The reduced pilot subsequently rejected this operator. In job 55966549 the A2
energy-only polycrystal increased from `5172.43154` to `5220.17837`, with 14
observed positive transitions and a maximum increase of `24.70799`. The state
remained finite, conservative, nonnegative, morphologically stable, and exactly
restartable. A0 passed the same checks. Because the failing control used
isotropic mobility, the defect lies in the anisotropic energy/update path and
cannot be attributed to inverse-correlated mobility. The remaining pilot and
all production work are closed pending a new thermodynamically consistent
polycrystal operator correction.
