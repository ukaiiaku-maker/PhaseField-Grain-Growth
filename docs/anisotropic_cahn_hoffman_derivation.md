# Anisotropic network first variation and PF coupling gate

## Convention and references

All new geometry uses Cartesian (x,y). Legacy pixel records use (row,column)
and require an explicit coordinate conversion. For an oriented edge A to B,
t = (B-A)/|B-A|, n = (t_y,-t_x), and theta = atan2(n_y,n_x).
Thus t = R(pi/2)n. Orientations are fixed scalar radians.

The source references are Hoffman and Cahn, *Surface Science* 31 (1972),
368–388, [DOI](https://doi.org/10.1016/0039-6028(72)90268-3), and Cahn and
Hoffman, *Acta Metallurgica* 22 (1974), 1205–1214,
[DOI](https://doi.org/10.1016/0001-6160(74)90134-5). Their bibliographic identities
are confirmed by [Cahn's NIST publication list](https://www.ctcms.nist.gov/~cahn/publications.html).
Publisher full text was not accessible during this session. The following
discrete derivation is independent and is checked directly by finite variation;
it is not represented as a full-text literature review.
Direct publisher abstract-page retrieval also returned HTTP 403 for both papers.
An additional search found the [TMS selected-works contents](https://www.tms.org/pubs/Books/PDFs/01-416X/01-416X-0.pdf),
which lists the two reprints on pages 293–313 and 315–324. Full reprint retrieval
was unsuccessful; this does not remove the full-text review limitation.

## Constitutive derivatives

Let u = cos(psi)^p + sin(psi)^p, h = u^(1/p). Then

    h' = u^(1/p-1) [sin(psi)^(p-1) cos(psi) - cos(psi)^(p-1) sin(psi)]
    h+h'' = (p-1) [cos(psi) sin(psi)]^(p-2) u^(1/p-2).

For p=2, h=1, h'=h''=0 exactly. For even p>2, h+h'' is nonnegative,
including zero at crystal axes. The isotropic part of the rounded support
gives a strictly positive floor proportional to 1-lambda. Numerical dense
evaluation is still required, including all initial-network pairs.

The fourfold misorientation factor is g_min+(1-g_min)sin(2 Delta)^2.
The inclination support is averaged over both crystals. Ca uses periodic
4096-point angular quadrature. Cgamma and CM are initial-length-weighted
normalizations, returned as an immutable new law. Normalizations from an
unqualified polygonal reconstruction remain provisional; they must not silently
become the production constants. Inverse energy–mobility correlation is a
synthetic constitutive hypothesis, not an empirical universal law.
Reported distribution percentiles use cumulative edge-length midpoint weights.
The rank statistic is length-weighted Pearson correlation of average tied
sample ranks; this convention is recorded rather than presented as a unique
definition of weighted Spearman correlation.

## Exact discrete variation

For edge displacement d, energy is L gamma(theta). Its derivative with respect
to d is c = gamma t - gamma_theta n = R(pi/2)xi, where
xi = gamma n + gamma_theta t. Consequently the negative energy gradient at
an interior vertex k is c_out - c_in. Dividing by dual length gives the local
vector force density; projection onto the node normal gives signed pressure.
The sum of the integrated interior forces telescopes to c_last-c_first.

Open endpoints have negative gradients +c_first and -c_last. They are returned
separately and excluded from the interior divergence. Incident endpoint forces
sum at a TJ to the outward Herring resultant. Joining subsegments adds their
endpoint terms once; adding them again as separate TJ forcing double-counts.
Closed loops have zero resultant but nonzero local pressure. A counterclockwise
circle has outward normal and negative pressure -gamma/R, hence shrinks for
positive mobility. With kappa=-dtheta/ds,

    dc/ds = -(gamma+gamma_theta_theta)(dtheta/ds)n = kappa_gamma n.

For network gradient flow, dE/dt = -sum Ldual v_n^2/M minus junction
|VJ|^2/MJ, provided interior nodes move normally and TJ velocities follow the
full resultant. A discrete timestep needs its own dissipation/convergence test.

## Reconstruction and topology

The existing arclength tracker does not provide a variational network: its
nearest-unvisited fallback can jump, lengths count pixels, TJ assignments use
proximity, and ordinal segment IDs can transfer histories after reconnection.
It is therefore not used for new line forces.

The new reconstruction obtains the upper envelope of linearly interpolated
phases on a conforming periodic triangular grid. Pair equality lines are clipped
against all competing phases. Shared vertices form an explicit edge graph;
degree and pair checks reject unresolved junctions. Deterministic graph walks
partition every edge once and unwrap the torus with explicit winding vectors.
Disconnected same-pair components remain separate. Zero-length edges,
unmarked branches and ambiguous half-box edges are rejected.

The fixed triangle diagonal introduces a grid-direction preference. Normal,
curvature, energy, TJ geometry and morphology convergence must be established
before this reconstruction can supply production forces. No history transfer
is implemented or claimed by this geometry module.

## PF coupling remains an explicit release gate

### Work-conjugate triangular pullback

An additional independent module now differentiates the resolved triangular
network energy with respect to the phase values at its mesh vertices. At a
side crossing, let d=eta_i-eta_j and let u be the mesh-side direction. The
linearized constraint gives

    delta x = -u [sum_a bary_a delta d_a] / (grad d dot u).

For a TJ inside a triangle, use the two constraints eta_i-eta_j=0 and
eta_i-eta_k=0. If their spatial gradients form the rows of A, then

    delta x = -inverse(A) [sum_a bary_a delta(eta_i-eta_j)_a,
                           sum_a bary_a delta(eta_i-eta_k)_a].

Contracting each endpoint variation with its negative energy gradient gives
the phase derivative. The implementation accumulates these pair differences
antisymmetrically, so the sum over phases is zero at every mesh vertex.
Three small finite-variation tests pass, including resolved TJs and exact
isotropic reference subtraction. Mesh-vertex crossings and singular junction
constraints are explicitly rejected because their fixed-topology derivative
is unresolved. This pullback does not yet constitute a qualified PF update.

### Legacy update and remaining contract

The historical rate at phase i is

    M0 gate [gamma0 (lap_i sum_eta - eta_i sum_lap
                    + obstacle (count eta_i - sum_eta))
             + external_i - mean_local(external)].

The geometric sharp-interface pressure gamma0*kappa is not this discrete
rate. Adding (p_gamma-gamma0*kappa) without a work-conjugate projection from
normal displacement to eta does not prove that the total update has the desired
velocity or descends network energy. Likewise a scalar spatial mobility cannot
represent three different pair mobilities at a TJ.

For an already qualified pair representation with base velocity b_ij and
noncapillary pressure q_ij, the correction would have to satisfy

    delta_v_ij = gate_ij M_ij (p_gamma_ij+q_ij) - b_ij.

This identity is a target, not a validated implementation. A native variational
kernel or a work-conjugate reference correction must first establish: exact A0
nesting, conservative antisymmetric phase contributions, full-force mobility
scaling, shared activation-work pressure, TJ endpoint accounting, and energy
dissipation under grid/time/interface refinements. Production is prohibited
until those tests pass. Historical PF and QIU equations are currently untouched.
