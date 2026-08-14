# Implemented physics equations

All scalar energies in event kinetics are in eV, temperature is kelvin, lengths are in the configured simulation unit, and PF energy/mobility units are self-consistent nondimensional units unless a configuration supplies a physical conversion. `physics_dimension` changes geometric event closure only; the spatial solver is two-dimensional.

## Multiphase field

For order parameters \(0\leq\eta_i\leq1\), \(\sum_i\eta_i=1\),

\[
F_{\rm int}=\frac12\sum_i\int_\Omega\left[\frac{\kappa_\eta}{2}|\nabla\eta_i|^2+W\eta_i^2(1-\eta_i)^2\right]dA,
\quad \kappa_\eta=3\gamma w,\quad W=6\gamma/w.
\]

The isolated equilibrium interface has energy \(\gamma\), profile \(\eta=[1+\tanh(x/w)]/2\), and \(\int(\eta')^2dx=1/(3w)\). The chemical potentials and constrained Allen–Cahn dynamics are

\[
\mu_i=2W\eta_i(1-\eta_i)(1-2\eta_i)-\kappa_\eta\nabla^2\eta_i,
\qquad
\dot\eta_i=-L(\mu_i-\bar\mu)+L f_i,
\]

where \(\bar\mu\) is the local active-phase mean, \(\sum_i f_i=0\), and \(L=M_0/(3w)\). This normalization yields the sharp-interface limit \(v_n=M_0\gamma\kappa\). A Euclidean simplex projection enforces bounds and filling after each explicit step. Periodic or zero-normal-gradient boundaries are supported.

## Disconnection modes and driving

A mode is \(m=(\mathbf b_m,h_m,N_{v,m},\Delta G^\ddagger_{0,m},\nu_m,N_{s,m})\), with \(\beta_m=|\mathbf b_m|/h_m\). Its reduced activation work and rate are

\[
W_m^\ddagger=p_nV^\ddagger_{n,m}+\tau_mV^\ddagger_{\tau,m}+N^\ddagger_{v,m}\Delta\mu_v,
\quad
p_n=\Gamma\kappa+\psi+p_q,
\]

\[
\tau_m=\mathbf t_m\cdot\boldsymbol\sigma\mathbf n,
\qquad
r_m=N_{s,m}\nu_m\exp[-\max(0,\Delta G^\ddagger_{0,m}-W_m^\ddagger)/(k_BT)].
\]

The `max(0,...)` is a specified barrierless transition, not a numerical patch: above it the rate is continuously attempt limited at \(N_s\nu\). Rates are formed in log space. Curvature couples to signed \(h\) through normal activation volume; shear couples to \(\mathbf b\). All admissible parallel modes compete using \(r_\Sigma=\sum_m r_m\), \(P(m|event)=r_m/r_\Sigma\). Necessary combinations use added mean residence times, \(r_{\rm combo}^{-1}=\sum_mr_m^{-1}\).

The synthetic isotropic spectrum uses finite \(|\mathbf b|\) shells, directions, and signed \(h\):

\[
\Delta G_0^\ddagger=G_{\rm core}+C_b(|b|/b_0)^{p_b}+C_h(|h|/h_0)^2+\xi_m.
\]

This is explicitly a surrogate closure. Disorder \(\xi_m\) is drawn once when a quenched library is constructed.

## First passage and geometric encounters

Each stochastic clock draws \(E\sim\mathrm{Exp}(1)\) and fires when

\[
H(t)=\int_0^t r(t')dt'\geq E.
\]

Piecewise-linear rate integration and hazard-space event-time interpolation are used. Overshoot is retained, allowing multiple events in one PF step. A geometric clock analogously integrates \(dH_{\rm enc}=\lambda\,dQ\), with \(Q\) chosen as GB measure change, TJ path, swept area/volume, or \(|\beta dx_n|\).

Persistent \(K\)-hit completion is the sum of \(K\) exponential passages, so \(T_K\sim\Gamma(K,r)\), \(E[T_K]=K/r\), and \(CV=K^{-1/2}\). Packet-reset completion over hazard \(\Lambda\) obeys

\[
P(N\geq K)=1-e^{-\Lambda}\sum_{j=0}^{K-1}\Lambda^j/j!=P(K,\Lambda),
\]

where the last form is SciPy's regularized lower incomplete gamma `gammainc(K,Lambda)`.

## Shear memory and nonlocal mechanics

The reduced entity state follows event-discretized

\[
ds=\beta\,dx_n-s\,dt/\tau_s-\sum_k\Delta s_k,
\quad E_s=\tfrac12K_ss^2,
\quad \tau_{\rm int}=-\partial E_s/\partial s=-K_ss,
\]

\[
v_n=M[\Gamma\kappa+\beta\tau_{\rm int}+\psi].
\]

The sign therefore follows the energy gradient. Reverse-curvature motion is possible only when the internal term opposes and exceeds capillarity. The `qiu_full_field` backend instead accumulates symmetric event eigenstrain, applies a periodic isotropic Fourier incompatibility projector, and computes plane-strain stress \(\sigma=2\mu\epsilon^{inc}+\lambda\operatorname{tr}(\epsilon^{inc})I\), with the zero wavevector removed. It is a nonlocal independently implemented surrogate, not a bitwise port of Qiu's line kernel.

At a TJ, event Burgers increments are conserved in a persistent residual \(\mathbf B_{TJ}\), with

\[
E_{TJ}=\tfrac12\mathbf B_{TJ}^T\mathbf K_{TJ}\mathbf B_{TJ}.
\]

Strict mode combinations require zero residual and compatible net step. Finite-residual configurations retain the state until a later event relaxes it.

## Atomic-to-PF event kinematics

An event records \(\Delta x_n=N_{disc}h_m\), \(\Delta\mathbf u_t=N_{disc}\mathbf b_m\), and \(\Delta N_v=N_{disc}N_{v,m}\). Hidden displacement accumulates until a configured PF release quota is reached. For unit out-of-plane thickness,

\[
\Delta\gamma_m=\frac{B_{\parallel,m}L_{swept}}{A_{RVE}},
\qquad
\Delta\epsilon_v=\frac{Q_m\Omega_{pd}}{V_{RVE}}.
\]

## Excess volume and climb

For GB measure change,

\[
dV_{ex}=\delta_V^{GB}dA_{GB},\qquad
dN_{req}=\frac{\delta_V^{GB}}{\Omega_{pd}}|dA_{GB}|,
\]

\[
q=N_{required}-N_{accommodated},\quad E_q=\tfrac12K_qq^2,\quad\Delta\mu=K_qq.
\]

The local exchange flux is

\[
J_{ex}=J_0[\exp(\alpha\Delta\mu/k_BT)-\exp(-(1-\alpha)\Delta\mu/k_BT)],
\]

with Onsager limit \(J_{ex}\simeq J_0\Delta\mu/(k_BT)\). Exponentials enter a documented finite attempt/flux limit before floating-point overflow. Transport uses

\[
D=D_0e^{-Q_D/k_BT},\qquad \tau_{tr}=C_{tr}\ell_{tr}^2/D.
\]

The stochastic production representation is the explicit serial state chain nucleation \(\rightarrow\) exchange \(\rightarrow\) transport \(\rightarrow\) quota completion. Its fixed-rate mean is \(r_{nuc}^{-1}+r_{ex}^{-1}+r_{tr}^{-1}\), never \((r_{nuc}+r_{ex}+r_{tr})^{-1}\).

## Analysis

In 2-D \(R_i=\sqrt{A_i/\pi}\), \(R_A=\sqrt{\langle A\rangle/\pi}\). Scaling fits scan/optimize \(n\) in

\[
R^n-R_0^n=K_n(t-t_0)
\]

and report residual autocorrelation, local slopes, and realization bootstrap intervals. Activation fits use per-event units only:

\[
K_n=K_0e^{-Q_{app}/k_BT},\qquad Q_{app}=-k_B\,d\ln K_n/d(1/T).
\]

The analytical comparator also implements intrinsic/drag \(dR/dt=K\Gamma(R)/R\), Class-B completion \(\Gamma=P(K,\Lambda(R))\), exchange crossover \(\Gamma=1/(1+R/R_x)\), and series/parallel activity composition.

