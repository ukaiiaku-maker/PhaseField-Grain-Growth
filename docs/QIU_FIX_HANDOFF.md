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
- Later commits contain analysis/status tooling; they do not change the active
  production trajectory.

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

Live snapshots contain only atomically closed Parquet parts and are marked
`pooling_allowed=false`. The active jobs retain checkpoint cadence <=500 and
compact movie frames sufficient to reconstruct label, stress, and eigenstrain
evolution.

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
