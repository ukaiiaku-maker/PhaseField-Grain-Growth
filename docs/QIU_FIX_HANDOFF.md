# QIU full-field qualification handoff

This campaign is active. Do not start another full-size worker while both jobs
below are running, do not resume the canonical historical directory, and do not
pool live snapshots or superseded preflights as terminal evidence.

## Source

- Worktree: `/private/tmp/qiu-full-field-qualification-v1`
- Branch: `codex/qiu-full-field-qualification-v1`
- Remote branch: `origin/codex/qiu-full-field-qualification-v1`
- Selected backend/config:
  `FFT_EIGENSTRAIN_V2` / `configs/production/fft_eigenstrain_v2_qualification_900K.yaml`
- Production physics source: `8bb7837677e6ca36fc9a952daf6aef601190b550`
- Later commits contain analysis/status/provenance tooling; they do not change
  the active production trajectory. Commit `a236192` fixes repository-root Git
  discovery for future staged jobs.

## Active HPC3 ownership

Project directory:
`/Users/sdillon/HPC3/qiu-qualification-20260910`

| Purpose | HPC3 run | Slurm | Source | State at handoff |
|---|---|---:|---|---|
| Legacy forensic replay | `20260910T232437Z-nogit-f4803a` | 55930486 | `147141b` | running |
| Corrected seed 5101 | `20260911T003946Z-nogit-c53868` | 55932457 | `8bb7837` | running |

Use `hpc3 status`, then `hpc3 fetch` only when terminal. Verify the runner's
archive checksum and `output/all-files.sha256`, then copy the unpacked result
under the durable local qualification root. Never clean remote results
automatically.

Unused run `20260911T004157Z-nogit-7dd0b1` is only a local `PREPARED` plan and
must not be submitted. It was created while reconciling a delayed submission
acknowledgement; Slurm confirmed the original job before any duplicate launch.

## Durable data

- Local qualification root:
  `/Users/sdillon/PF-graingrowth/results/qiu_full_field_qualification_20260910`
- Working-copy root:
  `/private/tmp/qiu-full-field-qualification-v1/results/qiu_full_field_qualification_20260910`
- Full JUnit:
  `/Users/sdillon/PF-graingrowth/results/validation/qiu_fix_production_preflight_tests.xml`
- First corrected live snapshot:
  `hpc_live/20260911T012000Z-corrected` under the durable root.
- Latest corrected restart snapshot:
  `hpc_live/20260911T044539Z-corrected`, step 500/t=20/N=587, archive
  SHA-256 `9ec285d7213768c446ccc34b5996f352b510160c3f049b4cc5cf43c1b50d40fe`.

Live snapshots contain only atomically closed Parquet parts and are marked
`pooling_allowed=false`. The active jobs retain checkpoint cadence <=500 and
compact movie frames sufficient to reconstruct label, stress, and eigenstrain
evolution.

The active wrappers verified their exact detached commits and bundle hashes,
but their internal application manifests report `UNCOMMITTED` because the
pre-fix scripts queried Git from the parent staging directory. Independent
machine-readable attestations are retained at
`active_hpc_provenance_attestations.json` under the durable result root. The
legacy step-1500 recovery archive, including checkpoint and evolution frames,
is at `hpc_live/20260911T013800Z-legacy` with SHA-256
`b8ae0b4205a1074dacef9cef7633442ea7b6a6ba743145662952bdb34a5e182a`.
The newest legacy recovery archive is
`hpc_live/20260911T041940Z-legacy`: checkpoint step 3500/t=140, 12 retained
evolution frames through step 3000, and SHA-256
`d58b052f5f8ff46112deee886f0f669013024b531851e22bcaaf0aa5c6e2ec5b`.
It is nonterminal and must not be pooled as production evidence.

The refinement and seed-5102/5103 job directories are staged but unplanned and
unsubmitted. All three now contain the verified `a173624` bundle with SHA-256
`c0ba9161a8f519436deec75c7017f03614fa04de43328221996b9efd0dd545b3`.
The refinement changes both the external-increment target from 0.02 to 0.01
and the actual configured PF timestep from 0.04 to 0.02; it therefore remains
a factor-two test even if the external limiter is inactive.
`production_revision_equivalence.json` in the durable root independently
records the byte-identical mechanics/PF/config hashes across active corrected
commit `8bb7837` and staged commit `a173624`.

Recreate or verify the active-job attestations with
`scripts/attest_qiu_hpc3_source.py`; its generated records are under
`hpc_live/source_attestations` in the durable root. Before issuing the final
decision, run `scripts/validate_qiu_qualification_decision.py` against
`qualification_decision.json`; it refuses missing terminal roles, nonterminal
or nonpoolable production records, nonzero final test outcomes, ambiguous
source SHAs, and missing gate groups.

## Next execution order

1. Retrieve and audit the first terminal active job; retain remote data.
2. When one slot frees, launch seed-5101 target=0.01 through the full critical
   interval from the exact same initial state and physics revision.
3. Prepare/hash two additional 384x384 initial states with the same
   1068-to-800 equilibration/compaction protocol, then run the corrected model.
4. Run `scripts/analyze_qiu_full_field_qualification.py` on terminal paths.
5. Render legacy, corrected seed-5101, and at least one additional seed with
   `scripts/render_qiu_qualification_movie.py`; retain each `.frames.csv`.
6. Issue `qualification_decision.json`, update the live report and validation
   files, rerun the complete suite, commit, push, and open (but do not merge) a
   pull request.

The final avalanche conclusion remains deliberately unset until the full
corrected/refined evidence exists.

The renderer at commit `cc1461e` merges compact cadence frames with dense
diagnostic fields, prefers the dense archive at duplicate steps, streams data
to bound memory, and writes transition-resolved morphology/stress/eigenstrain
contact sheets. This is required for the legacy step-9000--10246 onset rather
than rendering only the coarse cadence frames. FFmpeg is not installed on the
workstation, so the verified fallback product is an animated GIF with the same
CSV frame/hash index and JSON metadata.
