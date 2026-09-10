# QIU full-field qualification status

Last update: 2026-09-10, Phase 0 complete.

## Current source state

- Branch: `codex/qiu-full-field-qualification-v1`
- Base/current commit before Phase-0 commit: `9f66c8d7a5a266687284d8da35aefbc6062808f7`
- Worktree: `/private/tmp/qiu-full-field-qualification-v1`
- Historical production source: `4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`
- Historical QIU run: canonical, read-only, all integration-manifest hashes verified
- Current scientific decision: historical QIU remains an unresolved non-self-similar transient; source identity is provisionally `legacy FFT eigenstrain surrogate`, pending formal mathematical gates

## Commands completed

- Recorded `git status --short`, branch, HEAD, ten-commit log, and worktree list.
- Verified initial NPZ SHA-256 and every canonical QIU essential-file hash plus the exact 21-frame aggregate.
- Ran `PYTHONPATH=src:scripts /opt/anaconda3/bin/python -m pytest -vv --junitxml=results/validation/qiu_fix_baseline_tests.xml`.
- Ran `scripts/qiu_phase0_baseline.py` against the immutable historical QIU run.
- Inspected the archived `functions_2ref.py` pathways for `beta`, line-source construction, stress kernels, field extension, and PF feedback.
- Visually inspected the historical transition maps and terminal contact sheet.

## Evidence obtained

- Baseline: 178 passed, 0 failed/errors/skipped; JUnit SHA-256 `67267a5132d122663af7843e2f7615b46579511b2649f43150fb4a9a73ba8da4`.
- The exact historical transition interval is 446 solver steps (9,800 to 10,246).
- Interfacial energy is 14,854 at step 9,900 and 64,319 at step 10,000.
- Boundary-shear maximum changes from 52.7 at step 9,000 to `3.24e13` at step 10,000.
- Disconnected-grain fraction changes from 0 to 0.827 over the same saved-frame interval.
- No complete pre-avalanche restart exists; the final checkpoint cannot recover the onset.
- The archived reference is a current-geometry, orientation-derived line-source model. The historical production backend is an accumulated point-eigenstrain FFT surrogate with constant beta and is not presently entitled to an unqualified QIU-reference label.

## Decisions recorded

- Preserve the legacy implementation and historical output unchanged.
- Do not treat saved frames as restart states.
- Replay deterministically from the immutable initial state.
- Keep superseded diagnostic analyses with explicit exclusion manifests.
- Test the whole-grain-area/per-boundary source-overcount hypothesis rather than assuming it is the sole cause.
- If mathematical equivalence to the archived reference is not established, qualify the corrected backend as `FFT_EIGENSTRAIN_V2`, not `QIU_REFERENCE_V2`.

## Next automatic action

Implement read-only per-step diagnostics and the capture guard, prove diagnostic on/off trajectory invariance on a short run, and execute a deterministic legacy replay through the first guard trigger. Then independently derive and test the periodic elastic equilibrium operator.
