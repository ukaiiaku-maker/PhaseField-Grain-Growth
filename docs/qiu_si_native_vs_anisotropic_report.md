# Qiu SI native versus anisotropic extension report

Status: **interim implementation report**. Native reproduction and all
trajectory sections remain pending. Current classification:
`QIU_SI_ANISOTROPY_UNRESOLVED`.

## Identity and provenance

This report compares `QIU_SI_REFERENCE` with the new
`QIU_SI_REFERENCE_ANISO` extension. It does not use `FFT_EIGENSTRAIN_V2` or
`QIU_LEGACY_FORENSIC` as a control. The extension is new work, not a
reproduction of a published anisotropic Qiu result.

| Source | SHA-256 |
|---|---|
| `PF_Codes.zip` | `2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90` |
| archived four-reference driver | `03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17` |
| archived four-reference functions | `3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf` |

The sole native plan is owned outside this worktree as
`20260912T180419Z-nogit-a25f0c`. It has not been duplicated.

## Unchanged native quantities

| Quantity | Native and extension value/source |
|---|---|
| geometry and initial fields | archived `[100]_4ref` initialization |
| domain | 500x500, periodic |
| order parameters | 17 |
| spatial and temporal steps | `dx=dy=1`, `dt=0.1` |
| physical horizon | 200000 accepted steps |
| orientations | 0, -22.6, and 28.1 degrees in the archived 17-entry assignment |
| crystallographic references | four, spacing `pi/4` |
| beta | archived orientation-pair function |
| line/disconnection density | archived current-geometry reconstruction |

## Pair-energy correspondence

For each ordered native pair, define `s=phi_i+phi_j` and
`d=phi_i-phi_j`. Holding `s` fixed under conservative pair exchange, the
archived nine-point capillary expression satisfies
`capillary_ij = -s delta(E_native_ij)/delta(d)`. The factor `s` is therefore
part of the native Onsager operator.

The A2 port adds `E_A2-E_A0` to this pair energy and adds its force with the
same factor `s`; pair-specific mobility then multiplies the complete
capillary, elastic, and antisymmetric barrier drive. Both energy/force
identities pass central finite differences locally. A0 still executes the
archived algebra directly, so nesting does not depend on cancellation in a
rewritten formula.
| stress and elastic force | archived `sigma11`, `sigma12`, `sigma22` kernels and signed pair force |
| barrier term | archived antisymmetric `eij`, `eee=20` |
| boundary conditions and cadence | periodic; GB/stress output every 250 steps |
| renormalization | archived clip-and-normalize operation after the recorded pre-increment |

## Declared extension quantities

| Quantity | A2 candidate |
|---|---|
| pair energy | rounded fourfold Cahn--Hoffman law |
| `g_min` | 0.65 |
| inclination weight | 0.85 |
| support power | 16 |
| pair mobility | inverse correlated, exponent 2 |
| energy/mobility normalization | pending promoted native initial network |

The control rows are `QIU_SI_4REF_ISO_NATIVE`, `QIU_SI_4REF_A0_PORT`,
`QIU_SI_4REF_ANISO_E`, `QIU_SI_4REF_ANISO_M`, and
`QIU_SI_4REF_ANISO_EM_INV`. Their machine-readable template is
`configs/production/qiu_si_4ref_controls.json`.

## Equation correspondence and local evidence

The A0 port preserves the archived nine-point Laplacian and ordered-pair
capillary, elastic, and barrier sum. An independent direct loop produces an
exactly equal pre-renormalization increment. Each implemented unordered pair
adds equal and opposite increments. The A2 mobility multiplies that pair's
complete drive.

The Cahn--Hoffman addition is currently the discrete derivative of
`E_A2-E_A0`, counted once per forward bond. Its finite-difference derivative
test passes. Native triple-junction capillarity has not yet been derived from
the same pair energy for more than two present phases, leaving the full
variational correspondence gate unresolved.

## Pending native and trajectory evidence

The following sections cannot be populated before native promotion and are
not inferred from morphology: beta/reference records; reference-resolved line
density; three stress components; signed elastic force; pre/post-renormalized
increments; native energy quantities; matched checkpoint fields; initial
gamma, mobility, and stiffness distributions; energy-only, mobility-only, and
combined trajectories; morphology/facet and stress/line maps; TJ residuals;
timestep/grid convergence; cadence invariance; and restart equivalence.
