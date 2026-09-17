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
submission, the orchestrator submitted the preflight as Slurm job `56099932`.
The full baseline remains dependency-blocked until the retrieved preflight is
checksum-verified and classified `QIU_SI_NATIVE_PREFLIGHT_PASSED`.
