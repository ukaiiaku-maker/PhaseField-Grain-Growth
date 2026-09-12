# QIU full-field qualification: Phase-0 baseline audit

> **Naming correction (2026-09-12).** This historical filename is retained for
> provenance, but the audited trajectory is a Qiu-labeled in-house surrogate,
> not a Qiu baseline. The only true-baseline candidate is the still-unreproduced
> archived implementation now named `QIU_SI_REFERENCE`.

## Scope and immutable provenance

This audit is read-only with respect to the canonical Phase-1 campaign and its
historical Qiu-labeled surrogate result. The qualification worktree was created from audited HEAD
`9f66c8d7a5a266687284d8da35aefbc6062808f7` on branch
`codex/qiu-full-field-qualification-v1`. The historical production source is
`4761ef957715ba2faa84f015a0e4f4c4cd21c7aa`.

Canonical QIU run:

`/Users/sdillon/PF-graingrowth/results/long_time_kinetics_900K_20260824/20260825T012949Z-6d83ee5c82/QIU-T900-s5101`

The initial-state SHA-256 is
`106819770289d156b0ff6a35d356a9aa172258f428325ec50c5f62f380ea56f6`,
as required. Every essential QIU artifact matches `integration_manifest.json`:

| Artifact | SHA-256 |
|---|---|
| `manifest.json` | `8b5d93811657698d546a97080da6ae5459602671e862224cb1e06a0aee6ad50c` |
| `checkpoint.json` | `310e6a6810300f8056cf7aa792a87fcfc265b30917e17d725082d8b631746e37` |
| `checkpoint.npz` | `93006d857d525e4d6066e2612974af30f3af044b19f448150b0b2826f24d6255` |
| `grain_tracks.csv` | `babc3d35f8a22be4affacf7763a8012d5b985386bc2e25565b76bc50d9c7399f` |
| `boundary_tracks.csv` | `5a924273913540028ea22edc48b55ccf4bc6c4266f4d6e85d848f4345e226e36` |
| `mechanism_state.csv` | `83a9a2915f5d60e8446eeb9f67048ad8152659fe4a863e738f66562c2e1ce39b` |
| `energy.json` | `31e4ec35f31d102c41a5095396a2bc32e492f0e0696974c137eee781a5b34516` |
| 21 `frames/*.npz`, aggregate name/content digest | `416798f3b552a2b97f6cf289f4901f28244ca0740dac8bcee9db4a01d5b7bad1` |

The aggregate is calculated exactly as the integration tool does: update the
SHA-256 stream with each sorted filename followed by the raw bytes of that
file's SHA-256 digest, with no separator.

## Baseline software validation

The command

```text
PYTHONPATH=src:scripts /opt/anaconda3/bin/python -m pytest -vv \
  --junitxml=results/validation/qiu_fix_baseline_tests.xml
```

passed **178/178 tests** with zero failures, errors, or skips in 120.92 seconds
on Python 3.13.5. JUnit SHA-256:
`67267a5132d122663af7843e2f7615b46579511b2649f43150fb4a9a73ba8da4`.
The longer time than the earlier warm-tree audit includes fresh worktree/JIT
compilation and is not a scientific discrepancy.

## Historical transition

The historical change is from 494 grains at step 9,800 to 99 grains at step
10,246: **446 solver steps**, not 246. The following metrics use the campaign's
compactness definition, `P/sqrt(4*pi*A)`. Connected components and aspect ratios
are evaluated directly from periodic saved label fields.

| Step | N | mean/p95/max compactness | mean/p95/max aspect ratio | disconnected grains | max components/grain | RMS/max boundary shear | interfacial energy |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9,000 | 495 | 1.337 / 1.400 / 1.460 | 1.280 / 1.560 / 1.738 | 0 | 1 | 0.589 / 52.7 | 14,857.7 |
| 10,000 | 451 | 3.316 / 4.957 / 6.542 | 3.556 / 7.600 / 27.783 | 373 | 11 | 1.516e12 / 3.240e13 | 64,319.1 |
| 10,200 | 135 | 4.178 / 6.788 / 8.086 | 5.605 / 15.063 / 40.532 | 122 | 20 | 3.167e12 / 3.452e13 | 43,781.5 |
| 10,246 | 99 | 4.366 / 7.350 / 7.650 | 6.052 / 14.822 / 60.473 | 90 | 24 | 3.210e12 / 4.891e13 | nearest saved energy 43,781.5 |

At step 9,800 the mean/p95/max compactness remains
1.337/1.400/1.487 and the interfacial energy is 14,838.8. The energy remains
14,854.0 at step 9,900, then rises by a factor of 4.33 by step 10,000. The
grain count falls by 43 from step 9,800 to 10,000, by 316 in the next saved
200-step interval, and by another 36 before termination. The largest spatially
connected cluster inferred from the last pre-loss centroids contains 282 of
the 316 grains lost between saved steps 10,000 and 10,200.

Visual inspection of the historical movie, transition label maps, and terminal
contact sheet agrees with the quantitative record: the polygonal network at
step 9,000 becomes diagonal, multiply disconnected streaks by step 10,000.
This precedes the largest population collapse.

The first quantity identifiable at historical cadence is simultaneous stress
and energy blow-up. The maximum saved boundary shear changes from 52.7 at step
9,000 to `3.24e13` at step 10,000, while interfacial energy increases rather
than dissipates. The historical records do not contain full elastic energy,
per-step clipping, or per-step extinction data, so the precise first departing
solver step requires deterministic dense replay.

## Restart audit

There is no scientifically complete pre-avalanche restart. The only checkpoint
is terminal at step 10,246. It contains `eta`, `eigenstrain`, `active_phases`,
mobility, driving, orientations, previous entity state, and embedded JSON, but
cannot be used to replay the onset. The saved frames contain labels and derived
fields, not `eta`, accumulated eigenstrain, active-phase state, or complete
solver metadata. A frame is therefore not promoted to a restart.

The required forensic run is an exact deterministic replay from the unchanged
initial-state NPZ. Ordinary output can be used through a verified checkpoint
near step 9,000; per-step scalar diagnostics and ten-step field snapshots begin
there. Because the historical QIU path has no stochastic event/compatibility
modules, any divergence before instrumentation is enabled is itself evidence
of hidden nondeterminism or a trajectory-changing diagnostic.

## Initial model-identity and source-mapping findings

Read-only source comparison already rejects the claim that the historical
backend is a faithful port of the archived Qiu implementation:

- the archived code reconstructs a current-geometry line-disconnection density
  from the full boundary line and crystallographic reference decomposition;
- it derives two coupling factors from grain orientations and periodic
  misorientation, rather than using one constant `easy_beta=0.35`;
- the historical Python backend accumulates a history-dependent symmetric
  eigenstrain at one Eulerian midpoint per boundary and samples/broadcasts one
  pair force over the boundary;
- whole-grain area changes are converted to a normal displacement separately
  for every neighboring boundary domain, creating a concrete source-overcount
  hypothesis that must be tested with manufactured geometries;
- the historical FFT code describes its constitutive response as plane strain,
  whereas the archived Qiu code works with in-plane stress kernels and does not
  establish equivalence to that accumulated-eigenstrain closure.

At Phase 0, before the Fourier equilibrium, energy-gradient, source-integral,
and discrete-work gates were run, the historical model remained an unqualified
FFT eigenstrain surrogate and its Phase-1 transient was unresolved. The later
terminal source/work audit now rejects the surrogate and closes its initiating
numerical singularity; see the live forensic report.

## Durable artifacts

The accepted compact Phase-0 results are under
`results/qiu_full_field_qualification_20260910/20260910T160500Z-phase0-baseline-r4/`.
Earlier r1-r3 analysis attempts are retained with exclusion manifests; they
were superseded only because their compactness estimator did not match the
campaign definition. No historical artifact was changed.
