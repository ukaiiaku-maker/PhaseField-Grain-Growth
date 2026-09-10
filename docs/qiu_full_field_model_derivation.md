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

## Unique local interface sweep for `FFT_EIGENSTRAIN_V2`

At each pixel let `Delta eta_p=eta_p^{n+1}-eta_p^n`. Negative increments are
donors and positive increments are receivers. The unique local transfer from
donor `i` to receiver `j` is

\[
 q_{i\to j}=\frac{(-\Delta\eta_i)_+(\Delta\eta_j)_+}
                   {\sum_k(\Delta\eta_k)_+}.
\]

It follows exactly that each donor and receiver marginal is counted once and
that the absolute swept area is
`sum_x sum_j (Delta eta_j)_+ dx^2`. The construction does not iterate over
tracked boundary domains, so adding neighbors, splitting one pair into several
disconnected segments, or changing the tracking-domain length cannot multiply
the source.

For the canonical ordered pair `i<j`, the normal points from `i` to `j` and is
computed from the centered gradient of `eta_j-eta_i`; `t=(-n_x,n_y)` in stored
`(y,x)` order. With the orientation-derived directed coupling `beta_ij`,

\[
 B_{ij}(x)=\beta_{ij}\,\operatorname{sym}(t\otimes n),\qquad
 \Delta\epsilon^*(x)=\sum_{i<j}q_{ij}^{signed}(x)B_{ij}(x).
\]

The same tensor gives the feedback. Since

\[
 \frac{dE_{el}}{dq_{ij}}=-\sigma:B_{ij},
\]

the phase potentials are assigned symmetrically so that
`f_j-f_i=sigma:B_ij`, with `f=-delta E/delta eta`. This eliminates the old
midpoint-sample/whole-segment broadcast gain. A central finite-difference
energy perturbation agrees with this predicted work below `1e-4` relative
error. Reversing transfer reverses the source exactly; stationary fields add
none. Exact periodic rigid translations within the configured search radius
are detected and add no plastic source. Accumulated eigenstrain persists when
a phase or boundary disappears; topology changes never delete mechanical
history.

## Coupled time integration

Only `fft_eigenstrain_v2` changes timestep behavior. Before every trial, the
maximum centered external phase rate is evaluated on exactly the same local
phase support as the production obstacle kernel:

\[
r_{ext,max}=\max_{p,x}|M_0m(x)(f_p-\bar f_{\mathcal A})|,
\qquad
\Delta t_{ext}=\frac{\Delta\eta_{target}}{r_{ext,max}}.
\]

The first trial uses the minimum of the configured timestep, intrinsic
capillary bound, external bound, and any requested terminal-time bound. After
the PF trial, the unique sweep is mapped, eigenstrain is updated, and elastic
equilibrium is resolved. The complete interfacial-plus-elastic energy is then
compared with its pre-step value. An increasing trial is rejected, all PF and
mechanical arrays and clocks are restored, and the timestep is halved. Only an
accepted state updates entity bookkeeping and output. The configured default
`external_delta_eta_target=0.02` is an accuracy bound, not a stress cap; stress
and source remain uncapped.

Tests force the external limit active, prove exact checkpoint/restart with
variable timesteps, verify complete-energy monotonicity, bound raw pre-clipping
increments, and show monotonically decreasing matched-time solution error as
the target is tightened through four levels. Output cadence and read-only
diagnostics remain trajectory invariant.

### Historical overcount evidence

The legacy implementation computes one whole-grain area/perimeter displacement
inside every tracked GB domain, then deposits one unweighted point source for
each domain. It therefore changes if the same physical boundary is partitioned
differently. In the immutable historical run the mean/max domains per grain
pair were 1.374/3 at step 9000 and 1.376/3 at step 9800, but rose to 2.228/13
at step 10000, 3.619/33 at step 10200, and 4.341/34 at step 10246. Fragmentation
therefore multiplies further source deposition and supplies a concrete positive
feedback channel. The local-transfer map removes that dependency by
construction; it does not assume this channel is the sole avalanche cause.
