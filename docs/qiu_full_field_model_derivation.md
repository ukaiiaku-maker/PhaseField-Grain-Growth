# Full-field model identity and periodic elastic derivation

## Model identity

The historical `QIU` trajectory used an accumulated, Eulerian, point-deposited
eigenstrain field followed by a Fourier projection. The archived Qiu code does
not use that state variable. It reconstructs two orientation-dependent line
disconnection densities from the current grain-boundary network on every
update, superposes periodic real-space edge-dislocation stress kernels, and
projects that stress locally into the phase-field equation. Accordingly:

- `QiuFullFieldLegacy` means the preserved historical surrogate and remains
  available through the old `QiuFullField` import alias.
- `FFTEigenstrainV2` means the corrected periodic microelasticity problem
  derived below. It is not represented as a faithful Qiu implementation.

The archived kernels use the familiar two-dimensional edge-dislocation stress
shape with a fitted/scaled coefficient `G`; they do not expose an independently
identifiable Poisson ratio or a selectable macroscopic plane condition. The
FFT surrogate therefore makes its own `plane_stress` versus `plane_strain`
choice explicit. `plane_strain` preserves the constitutive choice made by the
historical FFT surrogate; it does not establish reference-model equivalence.

## Coordinates and constitutive law

Tensor index 0 is row/y and tensor index 1 is column/x everywhere in the FFT
backend. Thus `wave = (k_y, k_x)`, matching NumPy array axes. For an isotropic
two-dimensional law,

\[
 \sigma_{ij}=2\mu e_{ij}+\lambda e_{kk}\delta_{ij},\qquad
 e=\epsilon^c-\epsilon^*.
\]

The in-plane Lamé coefficient is

\[
 \lambda_{\rm strain}=\frac{2\mu\nu}{1-2\nu},\qquad
 \lambda_{\rm stress}=\frac{2\mu\nu}{1-\nu}.
\]

## Independent Fourier solution

For every nonzero wavevector,

\[
 \epsilon^c_{ij}=\frac{i}{2}(k_i u_j+k_j u_i),
 \quad k_j C_{ijkl}(\epsilon^c_{kl}-\epsilon^*_{kl})=0.
\]

Define

\[
 b_i=k_jC_{ijkl}\epsilon^*_{kl},\qquad
 A_{il}=\mu |k|^2\delta_{il}+(\lambda+\mu)k_i k_l.
\]

Then `i A u = b`, so `u = -i A^{-1} b`, with

\[
 A^{-1}=\frac{P_T}{\mu |k|^2}
          +\frac{P_L}{(\lambda+2\mu)|k|^2}.
\]

Substitution gives the compatible strain and the equilibrium stress directly.
At `k=0`, the production closure is traction-free mean strain:
`epsilon^c(0)=epsilon*(0)`, hence `sigma(0)=0`. A separately named
`fixed_mean_strain` option is retained for manufactured tests and future
explicit loading, but is not the selected production condition.

The stored elastic energy is

\[
 E_{el}=\frac12\int(\epsilon^c-\epsilon^*):C:
                 (\epsilon^c-\epsilon^*)\,dA\ge0,
\]

and its eigenstrain derivative is

\[
 \delta E_{el}=-\int\sigma:\delta\epsilon^*\,dA.
\]

`stress` in `FFTEigenstrainV2` is the physical Cauchy stress in these equations,
not its negative derivative. This sign convention is essential when the
source/force work-conjugate mapping is constructed.

## Gate result

The preserved historical projection is not an equilibrium Green operator. For
a fixed random symmetric eigenstrain on a 17 by 19 periodic grid its normalized
modewise equilibrium residual is approximately `0.69`, far above `1e-10`.
This is a fundamental mathematical failure, independent of whether the late
avalanche is also amplified by source overcounting or explicit coupling lag.

The corrected operator is tested for both plane conditions against:

- modewise equilibrium and stress symmetry;
- linearity, sign reversal, and translation;
- zero stress for uniform and compatible eigenstrains;
- nonnegative energy and three random energy directional derivatives;
- a dense global periodic displacement minimization on small odd grids; and
- a manufactured single Fourier mode.

All tests use the acceptance tolerances in the qualification directive. The
dense oracle constructs a global spectral strain-displacement matrix and solves
the displacement minimization by dense least squares; it does not reuse the
closed-form Green operator.

## Archived two-reference identity gate

For the polycrystal code, the directed misorientation is `theta=OR1-OR2` and
the archived `beta` returns two ordered coupling factors. Its three branches
are split at 70 and 150 degrees and are both scaled by 0.0333. The second factor
uses the complementary directed angle `sign(theta) pi-theta`. The local
inclination is rotated by the mean grain orientation and determines which of
the two orthogonal reference angles receives which ordered factor.

The production repository stores tensor axes as `(y,x)`, whereas the archived
functions name their first grid index `x`. The conversion is explicit in
`qiu_reference_geometry.py`; resolved shear maps `sigma11` to repository `xx`
(`stress[1,1]`) and `sigma22` to repository `yy` (`stress[0,0]`). Tests compare
the new factor function directly with the archived function extracted from its
AST, avoiding imports from the archive's unavailable legacy `sparse` stack.
