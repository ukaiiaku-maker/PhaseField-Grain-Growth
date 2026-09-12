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
  discovery for future staged jobs. Commit `6e7d9eb` provides the reusable,
  checkpoint-counter-based live-snapshot utility used at step 2000.

## Active HPC3 ownership

Project directory:
`/Users/sdillon/HPC3/qiu-qualification-20260910`

| Purpose | HPC3 run | Slurm | Source | State at handoff |
|---|---|---:|---|---|
| Corrected seed 5101 | `20260911T003946Z-nogit-c53868` | 55932457 | `8bb7837` | running; checkpoint 5500, latest closed step 5740/t=229.6/N=133 |
| Refined seed 5101 | `20260911T110924Z-nogit-6b926f` | 55950433 | `a173624` | running; checkpoint 1500, latest closed step 1820/t=36.4/N=480 |

Use `hpc3 status`, then `hpc3 fetch` only when terminal. Verify the runner's
archive checksum and `output/all-files.sha256`, then copy the unpacked result
under the durable local qualification root. Never clean remote results
automatically.

Both legacy jobs are terminal, locally retrieved, and remain preserved remotely.
The clean source/work continuation is
`20260911T111913Z-nogit-06fc67` / Slurm `55948258`, source `bebd53b`; its
726,135,940-byte terminal archive SHA-256 is
`f0c613ac63b358a92b1add3944202ef446b71a9f8863f2aec88a19fbd0899a23`.
The durable extraction is
`hpc_terminal/20260911T111913Z-sourcework-legacy-seed5101`. It reaches the
same step-9984/N=100 endpoint, all seven numerical checkpoint arrays are
bitwise identical to the first replay, and its 984 scalar plus 1,918,433
per-boundary rows close the source/work audit. The retrieval, source/work, and
trajectory-equivalence sidecars are colocated with that extraction. Neither
legacy replay is poolable as a qualified production result.

Unused run `20260911T004157Z-nogit-7dd0b1` is only a local `PREPARED` plan and
must not be submitted. It was created while reconciling a delayed submission
acknowledgement; Slurm confirmed the original job before any duplicate launch.

## Durable data

- Local qualification root:
  `/Users/sdillon/PF-graingrowth/results/qiu_full_field_qualification_20260910`
- Working-copy root:
  `/private/tmp/qiu-full-field-qualification-v1/results/qiu_full_field_qualification_20260910`
- Latest full JUnit:
  `/Users/sdillon/PF-graingrowth/results/validation/qiu_si_reference_split_tests_20260912.xml`,
  251/251 passed in 51.26 s, SHA-256
  `3e970b0a4220c5a9da20f7fe0c1da72ccfad77f142a1b1718a52d0ea393eadff`.
- First corrected live snapshot:
  `hpc_live/20260911T012000Z-corrected` under the durable root.
- Latest corrected restart snapshot:
  `hpc_live/20260911T152821Z-corrected`, step 2000/t=80/N=317, archive SHA-256
  `2a8e3d1ab77404c5efb0aef778496fb2f9ea86bb498480f53b7d00845ad57812`.
  It contains all 2,000 contiguous scalar rows in the 128 parts closed by the
  checkpoint recorder, 11 compact frames through step 2000, and dense fields
  through step 2000. The earlier step-1500 v2 recovery remains valid; its
  first filename-bounded archive remains excluded because it stopped at step
  1480.

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
`hpc_live/20260911T142720Z-legacy`: atomic checkpoint step 9500/t=380/N=495,
all 500 contiguous scalar rows from steps 9001--9500, 18 compact frames, and
101 bounded dense fields sufficient for ten-step movie reconstruction. Its
matching local/remote SHA-256 is
`a79f99d741883784edd8a897fbbce8e088a24f43d90582b0c5bfa676b673a4ba`.
It is nonterminal and must not be pooled as production evidence.
The terminal analyzer now records raw-capture provenance separately from
scientific transition acceptance and independently resolves the first
extinction, 100-step burst, morphology, energy, and nonfinite-field steps.
Through step 9500, the only such onset is a first total-energy increase at step
9002; all population and morphology indicators remain absent while the
source/stress/eigenstrain precursor grows.
All 18 common frames through step 9000 are byte-identical to the canonical
historical files, with a second array-level bitwise comparison also passing.
The checkpoint-bounded grain and boundary tracks through step 8000 are likewise
byte-identical after dropping only the intentionally different run-ID column.
The step-8000 equivalence JSON has SHA-256
`7277391fdd27816802a7ec4c86ab2578988024531ee80583c890d89ea9aa962d`.
The step-9001 audit explicitly classifies its zero-extinction clipping marker as
oversensitive and has SHA-256
`56c46c69da1d53528f46eaa91292445e9fdabfe7c9281856d3289bef865a00bf`.

The corrected step-1000 recovery has also been exercised through the final
renderer. Its seven-frame preview animation, SHA-indexed frame CSV, metadata,
and neutral initial/midpoint/terminal contact sheet are retained as
`hpc_live/20260911T082618Z-corrected/corrected-step1000-preview.*`. These prove
movie reconstruction from the retained cadence but are not terminal products.

The active legacy bundle predates the relative clipping-guard correction and
is expected to mark ordinary ~0.24 double-obstacle clipping near step 9001.
Its immutable wrapper does request continuation after a capture, so it should
continue to the historical endpoint; exclude only that oversensitive first
marker, not the continued trajectory. The old recorder will save one compressed
guard field per subsequent step, so monitor its 100-GB scratch allocation.
The first 16 dense-window rows show that this same active source also misses
closure-specific source/work hooks: those two columns are undefined and the
per-boundary diagnostic stream is empty. Its exact trajectory/fields remain
usable, but a second continuation is required for this gate.

Prepared plan `20260911T101424Z-nogit-6bb2f0` is superseded before submission
and must never be launched. Its replacement is immutable plan
`20260911T111913Z-nogit-06fc67`, which resumes the clean atomic step-9000 state
using source `bebd53b8706301715c216dafe224f2e0c4180aaa`, source-bundle SHA-256
`d024b52083f7f1c4042d4d8a99987cff8c5c17160499c8d2bbbc2d4d1f8cffb5`,
and clean recovery SHA-256
`bb7d4bdbd23c336f388ee28872155aea0ef280835bc87588edc8bae0b8c76e20`.
Its immutable HPC input SHA-256 is
`55e7d13a136b7876eb198ba382f285e34036847dcbb8008644d1b58ca4374108`.
A paired full-size two-step smoke is bitwise trajectory-identical to active
source `147141b`, while producing finite source energy/work and 4,078 finite
per-boundary records. Evidence JSON SHA-256:
`1a8b4b7d8a553e2984f5739170c773ef83bc633a28e10484c8f5aa0ab0763fa7`.

The refinement is prepared as immutable, unsubmitted plan
`20260911T110924Z-nogit-6b926f`; its HPC input archive SHA-256 is
`dfe68bb909d310bab479d328b3e76f9d1235a6b08d59c499d526f18a78d87583`.
Seed 5102 is immutable prepared plan `20260911T112155Z-nogit-6247ef`
(HPC input SHA-256
`8c7596aeff28ebd6a02919368f5a0dab2de86f4404f20922c2c7309273c72a9b`),
and seed 5103 is immutable prepared plan `20260911T112156Z-nogit-439d11`
(HPC input SHA-256
`08fdd599d2816544a469839bde5388fdea9ceeff4dc68b71c73a6aec9d78d90c`).
Both have no Slurm job and remain unsubmitted. All three corrected-job directories
contain the verified `a173624` source bundle with SHA-256
`c0ba9161a8f519436deec75c7017f03614fa04de43328221996b9efd0dd545b3`.
The authoritative machine-readable execution queue is
`hpc_execution_queue_20260911.json` under the durable qualification root,
SHA-256 `37d95dd44d77797fb9c6257d90f48774918b71a12f28fe60d94da08f50195801`.
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

1. Continue monitoring corrected and refined seed 5101 without exceeding two
   full-size campaign workers; retain remote data and verified snapshots.
2. Do not launch corrected seed 5102 or 5103 merely because a slot opens. Wait
   until both seed-5101 trajectories are terminal and locally retrieved, then
   require passing mechanics, work-conjugacy, complete-energy, and timestep
   comparison gates. Launch seed 5102 first only if that combined gate passes;
   seed 5103 follows only while the same gates remain satisfied.
3. Run `scripts/analyze_qiu_full_field_qualification.py` on terminal paths.
4. Render legacy, corrected seed-5101, and at least one additional seed with
   `scripts/render_qiu_qualification_movie.py`; retain each `.frames.csv`.
5. Issue `qualification_decision.json`, update the live report and validation
   files, rerun the complete suite, commit, push, and open (but do not merge) a
   pull request.

The final avalanche conclusion remains deliberately unset until the full
corrected/refined evidence exists.

The legacy trajectory is not a Qiu baseline. Its deterministic singularity is
closed forensic evidence and receives no further production allocation. The
true archived-code workstream is `QIU_SI_REFERENCE`; its source/input audit,
native-case blockers, staging hashes, and promotion gate are recorded in
`docs/qiu_si_reference_reproduction_status.md`. It must remain separate from
the independent corrected `FFT_EIGENSTRAIN_V2` qualification.

The renderer at commit `cc1461e` merges compact cadence frames with dense
diagnostic fields, prefers the dense archive at duplicate steps, streams data
to bound memory, and writes transition-resolved morphology/stress/eigenstrain
contact sheets. This is required for the legacy step-9000--10246 onset rather
than rendering only the coarse cadence frames. FFmpeg is not installed on the
workstation, so the verified fallback product is an animated GIF with the same
CSV frame/hash index and JSON metadata.

The later contact-sheet correction selects pre/transition/post panels from the
recorded `diagnostic_capture.json` step. In runs without a capture it uses
initial/midpoint/terminal roles, avoiding false event labels from ordinary
population loss between sparse frames.
The complete post-correction suite passes 229/229 tests; durable JUnit SHA-256
is `29a634fcaddb3f2a341a434d803f3e9ea8e25547d72c453b70616c748505eab5`.
For the legacy replay, use the renderer's `--transition-step` override added at
commit `167c30c`; this keeps the oversensitive raw step-9001 capture provenance
but centers the final contact sheet on the onset selected from per-step data.
