# QIU shared-core snapshot review

Read-only review at QIU commit
`847eb06e52c7276bf42ce9ab2f1a73ad74d21127`, whose worktree was clean at inspection.
The earlier dirty-state observations in the status history remain valid for
those earlier times; they do not describe this later snapshot.

The QIU status records immutable corrected production job 55932457 (account
SDILLON1_LAB, source 8bb7837), legacy replay 55930486 (SDILLON1), and an unused
PREPARED plan rather than an ambiguous duplicate. Both known jobs were observed
running. Neither was changed by this session.

| Shared area | Observed change | Anisotropic decision |
| --- | --- | --- |
| PF numerical kernel | Add opt-in trial statistics and maximum external-rate evaluation; retain the original accepted-step kernel | No accepted non-QIU equation correction to import. External-force limiting will need a dedicated anisotropy derivation and tests. |
| PF solver | Opt-in diagnostics, external-drive timestep utility, label caching | Diagnostics are optional. Existing uncached labels remain correct. The old path is retained while no anisotropic PF coupling exists. |
| Simulation | FFT-v2/source selection and coupled energy rejection; forensic recording | Guarded full-field functionality is outside this campaign. No FFT backend is imported. |
| Simulation memory | Hold previous out-of-place eta by reference instead of copying | Generic memory optimization, not a missing physical correction. Existing copies are retained; anisotropy resource profiling must measure their cost. |
| Kinematics / segment tracking | No changes in these files relative to the audited base | No generic correction available to resolve the anisotropic reconstruction or history-matching gates. |
| Checkpoint/restart | Added full-field and forensic state, cross-revision checks in QIU qualification runner | Do not import QIU state into a non-QIU campaign. Anisotropic normals, normalizations, physical-time schedules and trial state still need their own restart contract. |

No full QIU merge or cherry-pick was made. This review does not freeze another
session's future work. Recheck its HEAD, durable status and job ownership before
any production release, and repeat affected qualification after a necessary
integration. The current anisotropic production gate is closed independently
because PF coupling and reduced validation remain unresolved.

At 2026-09-11 01:24 UTC, a further read-only snapshot found HEAD
`4ef93905937b7bfec6f87a5497fe7411d5924acf` with an uncommitted change confined
to `scripts/run_qiu_spatial_convergence.py`. Commits since the reviewed snapshot
changed QIU diagnostics, qualification scripts and documents, with no further
shared PF kernel/solver/simulation change. The worktree is no longer globally
clean, so this snapshot is not a production release freeze. The primary audited
checkout remains clean at `9f66c8d7a5a266687284d8da35aefbc6062808f7`.

The final read-only snapshot (2026-09-11 02:08 UTC) is clean at
`c68312799b4e41ba404cc0dd1ec3865bdcc9f9b8`. There are no source changes since
4ef9390, so the earlier shared-core relevance decisions are unchanged. Both
known QIU jobs remain running with their original IDs/accounts. No anisotropy
job remains active. No QIU changes were imported; the missing PF implementation
and reduced validation independently prevent production release.
