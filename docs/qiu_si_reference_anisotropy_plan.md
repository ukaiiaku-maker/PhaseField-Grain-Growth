# Qiu SI reference anisotropy plan

## Identity and scope

`QIU_SI_REFERENCE_ANISOTROPY` is a controlled extension campaign for the
pristine supplementary-information Qiu implementation. Its native control is
`QIU_SI_REFERENCE`; neither `FFT_EIGENSTRAIN_V2` nor
`QIU_LEGACY_FORENSIC` supplies its baseline. The existing large-polycrystal
`ANISOTROPIC_PHASE1` matrix remains separate.

The new port preserves native current geometry, initial order parameters,
orientations, crystallographic reference directions, beta functions,
line/disconnection reconstruction, elastic constants, stress kernels, elastic
PF force, boundary conditions, output definitions, and terminal physical time.
It introduces only pair-specific `gamma_ij(theta, phi_i, phi_j)` and/or
`M_ij(theta, phi_i, phi_j)`. The capillary contribution stays distinct from
native line-disconnection stress and the elastic force is neither scaled nor
duplicated.

## Baseline-first gate

The sole native baseline is the immutable archived [100] four-reference case:
500x500, 17 order parameters, `dx=dy=1`, `dt=0.1`, 200000 steps, four
references, and stress/GB output every 250 steps. The Qiu qualification session
owns its prepared HPC3 plan. This workstream will not copy or submit it.

Promotion requires checksummed matched-state records for beta, selected
reference directions, reference-resolved line density, `sigma11`, `sigma12`,
`sigma22`, signed elastic PF force, pre-renormalization delta eta, accepted
post-renormalization delta eta, every native energy quantity, and selected
morphology/stress frames at the archived plotting scale. Final-image agreement
is insufficient.

## Known input limits

- The [100] four-reference archive is self-contained and runnable.
- [111] drivers name missing modules and its continuation needs a missing
  `OP_t142500.npz`; both remain blocked.
- The published polycrystal realization is absent and the supplied generator
  is unseeded. Any new run must be named
  `QIU_SI_POLYCRYSTAL_NEW_REALIZATION`, never an exact reproduction.
- A fourfold A2 law will not be relabeled as a [111] reproduction.

## Candidate and normalization

The first [100] candidate is the preregistered rounded-fourfold A2 law:
`g_min=0.65`, `lambda=0.85`, support power 16, and mobility exponent 2.
Energy and mobility normalizations must be recomputed on the actual initial
four-reference boundary network. Every sampled pair/inclination must satisfy
`gamma > 0` and `gamma + gamma_theta_theta > 0`. Phase-1 800-grain
normalizations cannot be reused, and A3 remains unreleased.

## Control matrix

The controls are immutable once normalized:

| Case | Capillarity | Mobility |
|---|---|---|
| `QIU_SI_4REF_ISO_NATIVE` | exact archived SI | exact archived SI |
| `QIU_SI_4REF_A0_PORT` | port with anisotropy disabled | native scalar |
| `QIU_SI_4REF_ANISO_E` | A2 pair energy | native scalar |
| `QIU_SI_4REF_ANISO_M` | native isotropic energy | A2 pair mobility |
| `QIU_SI_4REF_ANISO_EM_INV` | A2 pair energy | A2 inverse-correlated mobility |

All cases share initial fields, orientations, references, beta, line sources,
elastic parameters, stress kernel, boundary conditions, outputs, and terminal
physical time. Cases will not be tuned separately.

## Coupling and correspondence gates

The port must derive capillarity from the discrete pair free energy used by the
PF update. Pair contributions are antisymmetric and conserve phase sum before
native renormalization. The discrete Cahn--Hoffman force is the derivative of
that same energy, each interior/TJ contribution occurs once, and pair mobility
multiplies the complete intended capillary-plus-elastic pair drive. No averaged
scalar TJ mobility or added weighted-curvature pressure is allowed.

Before a full anisotropic run, A0 must match native beta, references, line
density, stress, elastic force, pre/post-renormalization increments, and fields
at multiple checkpoints. Energy/force finite differences, endpoint/TJ closure,
planar equilibrium, Wulff relaxation, timestep and feasible grid refinement,
cadence invariance, restart, finite fields, clipping behavior, denominator
singularities, and accepted-time propagation must pass. Attribution compares
native, energy-only, mobility-only, and combined cases at matched physical
times.

### Current operator audit

The archived A0 update has now been transcribed as an independent controlled
port in `grain_growth_pf.pf.qiu_si`. Its nine-point stencil, beta function,
ordered-pair barrier term, native scalar mobility, and pre-renormalization
increment are tested directly. The implementation rejects non-antisymmetric
`eij` or elastic pair inputs and applies a pair mobility to the complete
capillary-plus-elastic-plus-barrier drive.

The A2 energy term is presently an explicit discrete variational correction
`E_A2 - E_A0`; it vanishes identically when anisotropy is disabled, counts
forward bonds once, and its force passes a finite-difference derivative test.
This does not yet close the full energy correspondence gate: the archived
native TJ capillary expression
`phi_j laplacian(phi_i) - phi_i laplacian(phi_j)` is not the same as the
chemical-potential difference of the port's pair energy when more than two
phases coexist. That distinction must be resolved or formally derived before
qualification. It is recorded here rather than silently replacing the native
TJ algebra.

The immutable control manifest is
`configs/production/qiu_si_4ref_controls.json`. Its normalizations and release
flags remain unset until the promoted native initial network is available.
`scripts/compare_qiu_si_matched_states.py` is prepared to compare the required
checksummed native/A0 state arrays without consuming images as numerical
evidence.
`scripts/freeze_qiu_si_normalization.py` refuses unpromoted or wrong-archive
network evidence, checks positive sampled gamma and stiffness, and writes a
new frozen manifest with the network hash and sampled distributions. The
checked-in template remains visibly incomplete until that evidence exists.

## HPC3 coordination and release

The two active `FFT_EIGENSTRAIN_V2` seed-5101 jobs are protected. No additional
full worker will be submitted while both slots are occupied. When a slot opens,
the owning Qiu session runs and retrieves the native baseline first; this
workstream then consumes it read-only, runs short A0 correspondence checks, and
only afterward runs reduced anisotropic controls. A 200000-step anisotropic
control requires every prior gate.

Separate output roots are mandatory:
`qiu_si_reference_native`, `qiu_si_reference_anisotropic`,
`anisotropic_phase1`, and `fft_eigenstrain_v2`. Reports must state that the
anisotropic result is a new extension, not a published Qiu anisotropic
reproduction. Current classification:
`QIU_SI_ANISOTROPIC_PF_UNRESOLVED`.
