# Native anisotropic phase-field operator

The anisotropic solver is opt-in through `PFConfig.anisotropy_strength`. A null
strength and `A0_ISOTROPIC` both call the historical Qiu kernel without changing
its arithmetic. This is the exact trajectory-nesting contract used for A0.

For a non-A0 strength, the implemented discrete functional is evaluated over
locally supported grain pairs. At each forward-difference cell,

\[
E_h={4\over w}\sum_c \Delta x^2\sum_{i<j}
\gamma_{ij}(\theta_{ij})
\left[\eta_i\eta_j-{w^2\over \pi^2 N_c}
\nabla_h\eta_i\mathbin{\cdot}\nabla_h\eta_j\right].
\]

`N_c` is the fixed local active-set count for that discrete variation. The pair
normal is the direction of `grad(eta_i-eta_j)`. The code differentiates every
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
Each exchange is antisymmetric before clipping or simplex normalization. The
energy contains all diffuse junction contributions, so no separate Herring or
endpoint force is added to PF evolution.

The last evaluated pair chemical-potential difference is also the capillary
pressure source for disconnection activation work. This removes the former
anisotropic inconsistency in which activation used `gb_energy * curvature`
while PF would have used a different thermodynamic derivative.

The active-set threshold makes the functional piecewise smooth. Force/energy
finite differences must therefore keep the support fixed, while topology and
extinction are tested separately. The stability bound uses a conservative dense
upper bound on the selected law's reduced coefficient `M(gamma+gamma'')`.

This document describes an implementation awaiting diffuse HPC3 qualification;
it is not a production validation claim.
