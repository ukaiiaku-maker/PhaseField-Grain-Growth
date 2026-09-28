# Native Qiu SI independence and preflight audit

The pristine native baseline input is `PF_Codes.zip`, SHA-256
`2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90`.
The archived four-reference driver and its sole local module have SHA-256
`03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17`
and `3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf`.

An AST import and call-path inventory found only the archived
`functions_4ref_new` module and external numerical/plotting packages. Neither
the driver nor that module imports `grain_growth_pf`, `compact_support`, or the
generic KKT implementation. The native path is therefore independent of the
compact-support tolerance decision.

Prepared run `20260917T021402Z-nogit-c1b280` is a bounded preflight. It changes
the extracted copy's step ceiling from 200,000 to 1,000 and appends a final
state plus completion marker. It verifies archive and source hashes before
execution. Its artifacts state `baseline_eligible: false`; it is not the
promoted native baseline. The immutable 200,000-step plan remains
`20260912T180419Z-nogit-a25f0c` and has no Slurm identity.

After an exact identity and remote-directory reconciliation found no prior
submission, the first preflight attempt, job `56099932`, stopped before any
scientific step because the selected base environment lacked `matplotlib`.
Its partial diagnostics were retrieved and the failure was recorded as
operational. The repaired run `20260917T022751Z-nogit-0eb60c` uses the existing
checksummed analysis environment (including `matplotlib 3.10.1`) and was
submitted as job `56100055` after another identity reconciliation. The full
baseline remains dependency-blocked.

The replacement reached initialization and stress-kernel construction, then
failed before its first accepted step when Numba compiled pristine
`find_gb`. At `functions_4ref_new.py:118`, `np.concatenate` receives a tuple
whose first item is a Python list of arrays and whose second item is a 2-D
array; Numba 0.61.0 cannot type that signature in nopython mode. The retrieved
archive SHA-256 is
`b521e60c78111ecde5fe34798a412a838e3f2660e7e863574b9d13a1aa0bf65e`.
The result is `QIU_SI_NATIVE_PREFLIGHT_FAILED` with zero accepted steps. The
full native plan was not submitted, and no native or anisotropic downstream
task was released.
