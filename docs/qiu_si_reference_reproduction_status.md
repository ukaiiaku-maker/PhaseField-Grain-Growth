# QIU_SI_REFERENCE reproduction status

Status: archive identity verified; native benchmark reproduction not yet
complete. This workstream is independent of `FFT_EIGENSTRAIN_V2` and of the
closed `QIU_LEGACY_FORENSIC` investigation.

## Model identities

| Name | Meaning | Current scientific status |
|---|---|---|
| `QIU_SI_REFERENCE` | Pristine archived current-geometry, orientation-derived line-disconnection implementation | true Qiu baseline candidate; not yet reproduced |
| `FFT_EIGENSTRAIN_V2` | Corrected periodic accumulated-eigenstrain model developed in this repository | independent qualification running; not Qiu-equivalent |
| `QIU_LEGACY_FORENSIC` | Historical in-house point-eigenstrain surrogate formerly labeled QIU | deterministically reproduced, diagnosed, and rejected |

The legacy source singularity is closed. No further production allocation is
authorized for reproducing it. Its artifacts remain immutable forensic
evidence only.

## Archive and environment audit

The archived Zenodo payload `PF_Codes.zip` passes both recorded identities:

- SHA-256: `2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90`
- MD5: `6cd49ca72eba89210abb96e700342f12`

The new audit inventories every non-AppleDouble member and rejects unsafe ZIP
paths. The durable machine-readable record is
`results/qiu_si_reference_20260912/archive_audit.json`, SHA-256
`fbc92ee46264cea21729f7c1d82aea25577b4298d04ff5d1a960a8ed0f2daa58`.
It was produced by `scripts/audit_qiu_si_reference.py` using a dedicated,
ignored environment with Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3,
Matplotlib 3.10.0, Numba 0.61.0, and sparse 0.19.2. Import probes of the exact
four- and six-reference function files pass. The probes return archived
coupling values `(0.3996394353984736, -0.3996394353984736)` for the four-
reference `(0, -22.6 deg)` pair and `(0.58, 0.58)` for the six-reference
`(90, 57.8 deg)` pair.

After adding the audit/staging tool and its archive-path/import tests, the
complete repository suite passes **251/251** tests in 51.26 seconds with zero
failure, error, or skip. JUnit SHA-256:
`3e970b0a4220c5a9da20f7fe0c1da72ccfad77f142a1b1718a52d0ea393eadff`.

## Native benchmark readiness

| Archived case | Native parameters | State | Blocking evidence |
|---|---|---|---|
| `[100]` four-reference idealized | 500x500, 17 OPs, dx=1, dt=0.1, 200,000 steps | staged, not started | none in archive; full run remains computationally large |
| `[111]` six-reference idealized | 500x500, 17 OPs, dx=1, dt=0.1, 200,000 steps | compatibility repair required | driver imports absent `functions_6ref_new2_1`; archive contains `functions_6ref.py` |
| `[111]` continuation | same grid/timestep; manual continuation | blocked | driver imports absent `functions_6ref_new2_2` and requires absent `OP_t142500.npz` |
| polycrystal | 250x250, 1001 OPs, dx=10, dt=10, 300,000 steps | blocked for exact reproduction | required `PolycrystalSeeds_1000Grains_25_250_1.txt` is absent |

The polycrystal archive includes a random seed generator, but it initializes
neither Python's `random` state nor NumPy's state. A newly generated file is
therefore a new realization, not the missing published native input, and must
not be represented as exact reproduction.

The exact four-reference driver and support file are staged under ignored
`.external/qiu/staged/four_ref_native`. Their hashes are respectively
`03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17`
and `3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf`.
The stage manifest SHA-256 is
`c40fa6a2394d5a9db8628a829c1b65591ec1caff4672a35dbde702b8c6bcb8a7`.
No source line or native parameter has been changed.

## Promotion gate

No archived or independent result is presently called the Qiu baseline. The
four-reference run must first complete with exact source and environment
attestation, after which the following native observables must be checksummed
at matched states:

1. orientation-pair coupling factors `beta`;
2. reference-resolved line/disconnection density;
3. `sigma11`, `sigma12`, and `sigma22` fields using the archived plot scaling;
4. the signed elastic phase-field force;
5. the explicit pre-renormalization and accepted `Delta eta` arrays.

The six-reference compatibility aliases, if used, must remain declared repair
artifacts and be tested against the included function-file APIs before their
outputs are compared. The native polycrystal reproduction remains blocked
until the published seed realization is recovered or the comparison is
explicitly downgraded to a new-realization study. Transplantation to the
384x384 Phase-1 geometry is prohibited until this native gate passes.

## Resource scheduling

No third full worker is launched while corrected seed 5101 and its factor-two
refinement occupy the two protected HPC3 slots. The native four-reference
benchmark has an immutable, unsubmitted HPC3 plan
`20260912T180419Z-nogit-a25f0c`, input-archive SHA-256
`2783b5c019323ee0fb04742e35226c18938b5ddc7b704a3990ecc6025097dbcd`.
It requests 8 CPUs, 24 GB, 100 GB node-local scratch, and a 10-day ceiling;
the wrapper re-verifies the ZIP and both source files before executing the
unmodified archived driver, then archives all outputs for local retrieval.
Corrected seeds 5102
and 5103 are also held: they may start only after both running seed-5101 jobs
are terminal, locally retrieved, and have passed mechanics, work-conjugacy,
complete-energy, and timestep-comparison gates.
