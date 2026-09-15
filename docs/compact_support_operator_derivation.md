# Compact-support anisotropic operator

## Discrete state and energy

At every cell (x), the phase vector belongs to the Gibbs simplex

\[
K=\{\eta:\eta_i\geq 0,\ \sum_i\eta_i=1\}.
\]

The discrete free energy (F_h) is unchanged from the homogeneous-norm A2
operator. Its chemical potential μ is the exact nodal derivative already
verified by finite differences. In particular, the support repair neither
omits low-amplitude energy terms nor weakens the A2 law.

For the exact support (S_x=\{i:\eta_i(x)>0\}), define the candidate set

\[
C_x=S_x\cup S_{x-e_1}\cup S_{x+e_1}\cup S_{x-e_2}\cup S_{x+e_2}.
\]

This graph is a deterministic function of the checkpointed field. A phase can
therefore enter only from the spatial stencil. There is no amplitude cutoff
and no remote nucleation path.

## Pair-mobility descent direction

For every unordered pair in (C_x), the existing symmetric nonnegative pair
mobility (M_{ij}) defines the antisymmetric exchange

\[
v_i^{ij}=M_{ij}c_h(\mu_j-\mu_i),\qquad
v_j^{ij}=-v_i^{ij}.
\]

Thus (v=-L_M\mu), where (L_M) is a weighted graph Laplacian, and

\[
\mu^Tv=-\frac12\sum_{ij}M_{ij}c_h(\mu_i-\mu_j)^2\leq0,
\qquad \sum_i v_i=0.
\]

The pair mobility continues to multiply the complete generic pair drive.

## Double-obstacle proximal active set

For trial scale α, each cell solves

\[
\eta^{n+1}_x=\operatorname*{argmin}_{z\in K,\ z_i=0\ (i\notin C_x)}
\frac12\left\|z-(\eta^n_x+\alpha\Delta t\,v_x)\right\|_2^2.
\]

The solution is the exact active-set simplex projection

\[
z_i=\max(y_i-\lambda,0),\qquad
\sum_i z_i=1,
\]

where (y=\eta^n+\alpha\Delta t,v). With

\[
r_i=z_i-y_i+\lambda,
\]

the obstacle KKT conditions are

\[
z_i\geq0,\qquad r_i\geq0,\qquad z_ir_i=0.
\]

The recorded KKT residual is the maximum primal, dual, complementarity, and
phase-sum defect. Values `1e-8`, `1e-10`, and `1e-12` are the only admitted
solver tolerances. They govern qualification acceptance; they do not prune a
phase or delete an energy term.

For an unforced step, deterministic backtracking halves α until

\[
F_h(\eta^{n+1})\leq F_h(\eta^n)+64\epsilon_{mach}
\max(1,|F_h(\eta^n)|).
\]

The pair direction is a descent direction before the obstacle solve, the
projection enforces the variational inequality, and the final acceptance test
enforces discrete energy descent to roundoff. The result has exact zeros and
requires no post-update cutoff or global renormalization.

## Scaling and restart

Candidate construction enumerates exact support once, then every force and
pair operation uses the local padded graph. Pair work scales as

\[
\sum_x {|C_x|\choose2}
\]

rather than the square of the global grain count. Sorted phase indices fix
loop order. Since (C_x) is reconstructed solely from η and periodicity, a
checkpoint needs no hidden graph state and restart is deterministic.

This derivation qualifies the algorithmic structure. Release still requires
the prescribed 192×192 coarse/half-step HPC3 tolerance, restart, energy,
morphology, and scaling campaign.
