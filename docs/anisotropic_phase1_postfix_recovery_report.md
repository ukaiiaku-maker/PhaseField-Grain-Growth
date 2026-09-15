# Phase-1 anisotropic postfix recovery

## Recovery identity and scope

This recovery uses branch `codex/aniso-phase1-postfix-recovery-v1`, created in
an isolated worktree directly from `1316cc89dbabfb41cb883b0d4a4c74738cc2bef6`.
The executable scientific equations are unchanged from that commit. Recovery
edits are limited to checkpoint/finalization code, focused tests, and reports.
The later Qiu/SI source is absent and outside this recovery.

## Job 55968528 archive audit

The canonical local record is
`/Users/sdillon/HPC3/anisotropic-phase1-state/hpc3-results/pfgg-anisotropic-cahn-hoffman-v1/20260912T212636Z-nogit-6d3d72`.
The partial archive was fetched twice with matching SHA-256
`f4b97f553ec5e6923bb205eb95700ecec63a7fb60b74dd6c3cd94edf749980f9`.
Slurm recorded `TIMEOUT`, 04:00:25 elapsed against 04:00:00 requested, two
allocated CPUs, 03:58:45 CPU time, 49.65% CPU efficiency, and 1.19 GiB peak
memory. The application wrapper had not returned; emergency finalization
created a partial archive. The final marker therefore records incomplete
finalization and cannot confer scientific completion.

The archived identities are:

- scientific source commit: `1316cc89dbabfb41cb883b0d4a4c74738cc2bef6`
- submitted source archive: `57197ed07f845a4823a745403dcb227466864ba32b04db21091667db59e42f6e`
- runner input bundle: `8d922df923de1bf2d57a305e16e309830837a73c7fbad9e495604e4adf3c6c25`
- captured transition state: `75f189987c95a94645b628548078ecb9d37c460ada921b81c8ad233658558a0c`
- environment manifest: `714dd16b2f3624b877352abbbe3e75a576883507e8b286dfd0691eb70d30438a`
- resolved job manifest: `54f5fe3d606a37f9c1a3252625a5212e41e306e9691b6113f91ba70c9b3bff65`
- job wrapper: `81271fe933acc28b72b0f79d6ce1c18d81a17990309a0732caadfad8d8ac179c`

The archive contains a complete `pytest.log` (`206 passed in 85.67s`) and a
passing manufactured `qualification.json`. Those results are reusable upstream
evidence because the scientific source is identical. It contains only progress
text for the polycrystal calculation and an empty `output/postfix/` directory.
There are no phase fields, active masks, orientations, physical times, accepted
step records, anisotropy/configuration identities, pair-support state,
timestep histories, checkpoint schemas, energy traces, final states, case
summaries, or restart comparisons.

| Requested record | Preserved progress | Classification | Recovery action |
|---|---:|---|---|
| A0 continuous | 64/64 | `LOG_ONLY` | rerun from initial state |
| A0 midpoint restart | replay not independently serialized | `LOG_ONLY` | rerun from new midpoint checkpoint |
| A2 energy-only continuous | 64/64 | `LOG_ONLY` | rerun from initial state |
| A2 midpoint restart | replay not independently serialized | `LOG_ONLY` | rerun from new midpoint checkpoint |
| half-dt A2 continuous | 96/128 | `LOG_ONLY` | rerun all 128 steps |
| half-dt A2 midpoint restart | no state | `LOG_ONLY` | rerun from new midpoint checkpoint |

No exact scientific state from job 55968528 is reusable or restartable.

## Finalization repair and preregistered decision

Each recovery case writes atomic, fsynced checkpoints at steps divisible by 16
and at stage boundaries. The checkpoint embeds the complete field, active mask,
orientations, mobility scale, exact step/time, solver configuration, energy and
support histories, timestep records, schema, source, input, configuration, and
environment identities. Pair support is a pure derivation of the field and
active mask and has no mutable solver cache. A separate exact midpoint and
continuous final state support the restart replay. `CASE_INCOMPLETE` remains
until the case summary, energy CSV, timestep CSV, and `CASE_COMPLETE` marker are
durable. `USR1` requests an exact checkpoint and never creates a completion
marker.

The timestep result is preregistered as converged when the matched-horizon
coarse/half-dt final-energy relative difference is at most 1% and each of
boundary density, grain-area CV, and area-weighted mean-radius differs by at
most 10%. Classification priority is energy failure, restart failure, timestep
nonconvergence, complete qualification, then operational incompleteness.

## Sizing and preflight

Job 55968528 did not timestamp its individual progress records, so it cannot
support separate defensible measurements for coarse A2, pre/post-saturation
half-dt A2, or restart replay. The only exact rate is the aggregate bound:
03:58:45 CPU time and 04:00:25 elapsed while completing the tests,
manufactured/actual-state work, A0 continuous/restart, A2 continuous/restart,
and 96 half-dt continuous steps. Assigning that full elapsed time to the 224
logged continuous accepted steps gives a conservative upper bound of 64.4
elapsed seconds per logged step; it is not represented as a case-specific
measurement. The recovery records monotonic per-stage timings so this omission
will not recur.

The recovery repeats 256 continuous steps plus 128 restart steps, omits the
already passed full test/manufactured work and actual-state audit, and adds
per-step diagnostics plus atomic I/O. Its conservative predicted runtime is
eight hours; twelve hours supplies the required 50% margin. One CPU and 3 GiB
limit exposure to 12 CPU-hours. On 2026-09-14 the `standard` partition was UP
with a 14-day maximum, the SDILLON1 balance was 477 SUs, the remote project root
existed, and no job with the recovery name/comment appeared in `squeue`. One
unrelated FFT job, 55950433, was running. No Qiu recovery or production job was
submitted.

The executable recovery commit is
`ad698bc8ac3b6d26d04d07bb241b85e5f0562f97`. The immutable runner ID is
`20260914T202421Z-nogit-f0fd7d`, submitted once as Slurm job `56017421` with
job name `aniso-postfix-recovery-v1`. Its identities are:

- recovery source archive: `f9507613c56b664b0b9d57953cbd9ad27bba741494e8bc71e9b8e98288197cf9`
- runner archive bundle: `49400c14508e937f7335f9f2c65a8a969e72f8c96f3633f7c375239ee518ccae`
- logical input manifest: `8f614b6b07fc97e0257e46cebca2d0ce29e25eeb99b6d5cfefbf7c7480a3c2f9`
- resolved recovery config: `51daaf8fd254c5a320c93512d4fe356aeed0101b91bd4f3fdfc43634a0e41558`
- environment: `714dd16b2f3624b877352abbbe3e75a576883507e8b286dfd0691eb70d30438a`
- initial phase field: `2986f2bf744107c46aeedf05d00849a5f35db8c3ba85f70f7a7645b7922b2e61`
- orientations: `217578bc628ee08179cf59123c21a9af7b31ecdfc4eb8ee6a0816446221f3141`

The first remote-directory operation timed out before upload. Reconciliation
found no matching Slurm job and no remote bundle, so the same immutable
prepared run was retried and submitted; no second run ID was created.

## Completion result

Job 56017421 terminated after 11:57:27 with Slurm state `FAILED` and exit
`1:0`. This is the recovery program's deliberate incomplete exit after the
checkpoint request, not a numerical exception. The application exit is 1 and
finalization exit is 0. It used 11:31:49 CPU time on one CPU (96.43% CPU
efficiency) and 981068 KiB peak RSS. The requested resources were one CPU,
3 GiB, and 12 hours; actual exposure was 11.96 CPU-hours.

The result archive was fetched and independently compared with the remote
archive and runner checksum. All three SHA-256 values are
`b1d3eec6b0b693641958c78fcc73bf2bd7a83bc827eff861fdd661cfd116fc42`.
The final marker SHA-256 is
`0146d43016914e5846bc0606ccdc61d794b2b4aa52ced80c0c8156f6e556282e`.
The marker records `complete:false`, `application_exit:1`, and
`finalization_exit:0`.

The recovery recomputed every requested trajectory from the authoritative
deterministic initial state because job 55968528 supplied no exact state. It
resumed no state from that job. A0 and coarse A2 each completed both the
64-step continuous path and exact midpoint replay. Fine A2 completed its
128-step continuous path, then replayed from the midpoint through exact step
126 of 128. Its checkpoint contains the complete field, active mask,
orientations, mobility scale, physical time 0.004009908051346792, 62 replay
timesteps, the complete 129-point continuous energy trace, 129 support records,
and all submitted identities. It correctly retains `CASE_INCOMPLETE`; no fine
case summary or `CASE_COMPLETE` exists.

| Case | Energy, initial to final | Maximum increment | Positive steps | Restart result |
|---|---:|---:|---:|---|
| A0, 64 steps | 5210.112257121886 to 5191.6061749583605 | -0.1458073707026415 | 0 | exact; both field hashes `1d53555e5a2a345ca5bf04e22f1daf6ffbad6a67d663e2fe369819189e868850` |
| A2 energy-only, 64 steps | 5236.012485284901 to 5149.6346852197785 | -0.9399010361212277 | 0 | exact; both field hashes `705c5a32329aa697b51273a473c6a3143e1c12c35ff126bf6f915a4ed89a36f3` |
| A2 energy-only, half dt, 128 continuous steps | 5236.012485284901 to 5115.519299781972 | -0.47039203997974255 | 0 | incomplete at replay step 126/128; continuous final field `a12a42877337e1c2f71d52c1a6415ee424ace386f7abf3325e9070208d25bc62`, preserved replay field `52a8b65b8b7cdbd662d2dce7279a78664e812b2c0b7b33b4aa5ad4331ca66101` |

Every recorded state is finite and nonnegative, with maximum phase-sum error
4.44e-16 for A0 and 2.22e-16 for both A2 cases. The positively homogeneous
zero-gradient branch remained finite throughout. Coarse A2 reached 70 maximum
active phases per cell, 2,415 maximum active pairs per cell, 49,763,328
cell-local supported-pair instances, and 230,742 zero-gradient supported-pair
instances. Fine A2 grew to 184 active phases per cell, 16,836 active pairs per
cell, 494,818,360 supported-pair instances, and 3,581,251 zero-gradient
instances. Its smallest nonzero pair-gradient magnitude was 7.89e-31; all
evaluations remained finite. This growth was numerically finite but not
computationally manageable within the conservatively requested 12-hour job.

The completed coarse and fine continuous states satisfy the preregistered
timestep comparison: final-energy relative difference 0.0066690, boundary
density difference 0.0017495, grain-area CV difference 0.0006918, and
area-weighted radius difference 0.0001340. These are all below their 1% and 10%
thresholds. This does not qualify the model because the exact fine restart
comparison is missing two accepted steps.

Case summary SHA-256 values are
`91e3108d65c634df6eefa134fbee7bcb23d091791f508bec7dd07704744282d1`
for A0 and
`c283826a03a6c876643db4846466e377de651837285cd94c734038502cfefcd7`
for coarse A2. Fine A2 has no admissible case summary; its checkpoint JSON and
NPZ hashes are
`32963dd9bf37f5666790876e67c84ae851a9106b4d70ad3422a807e8a07ea6a3`
and `dfa0d876ea9d3106463b510ea6f5bf063f82dd11d39552918f5171329b205cfd`.

The generated submission included `#SBATCH --signal=B:USR1@900`, but Slurm
delivered that signal to the batch shell while it was waiting for `run.sh`; the
shell trap did not promptly invoke `checkpoint.sh`. An overlapping step on the
already assigned node invoked that uploaded checkpoint command once. The
scientific process then completed its current accepted step, atomically wrote
restart step 126, and exited before the hard limit. The wrapper copied the full
incremental output and finalized the archive successfully. The external runner
renderer was subsequently corrected to emit `#SBATCH --signal=USR1@900`, which
the installed HPC3 `sbatch` manual defines as signaling all job steps. Its three
focused renderer tests pass. This runner correction is outside the isolated
scientific-source branch.

The required final classification is exactly:

**`A2_POSTFIX_OPERATIONALLY_INCOMPLETE`**

The recorded energy, convergence, finiteness, and low-gradient evidence is
favorable, but the missing two replay steps prohibit an exact fine restart
decision. Accordingly the reduced A2 mobility-only and combined controls are
not released. No A3, production, or Qiu/SI job was submitted or modified.

The verified runner record is
`/Users/sdillon/HPC3/anisotropic-phase1-state/hpc3-results/pfgg-anisotropic-cahn-hoffman-v1/20260914T202421Z-nogit-f0fd7d`.
The independently extracted archive is
`/Users/sdillon/HPC3/anisotropic-phase1-state/extracted-20260915T0825Z`, and
the case records are under its `output/postfix` directory.
