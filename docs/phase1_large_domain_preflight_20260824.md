# Phase-1 large-domain preflight

## 432 by 432 preferred-domain attempt

The preferred 432 by 432, approximately 1000-grain initialization was started
from clean source SHA `427a29a` using `/opt/anaconda3/bin/python`, one process,
and the frozen 900 K physics.  Streaming Voronoi assignment and the exact
simplex fast path kept peak resident memory near 6.3 GB on the 32 GB host, so
memory was not the limiting resource.

The initializer remained inside the PF equilibration loop and had not reached
its first 100-step progress checkpoint after more than nine minutes.  Even the
optimistic upper bound of 100 steps in nine minutes projects to more than 6.25
days per 100,000-step trajectory.  Nine cases with at most two workers would
therefore require more than 28 days before entity tracking, event processing,
checkpoint compression, or rendering overhead.  This meets the directive's
explicit multi-week-wall-time fallback condition.  The diagnostic was stopped
with `KeyboardInterrupt`; it created neither a campaign root nor a reusable
partial initial-condition archive.

## One authorized fallback

Phase 1 consequently uses 384 by 384 with approximately 800 equilibrated
grains.  Its initial area per grain is

\[
384^2/800=184.32,
\]

which exactly matches the old 192 by 192, 200-grain value.  The precursor count
is initialized at 1068, preserving the validated 267/200 starting ratio; the
generator still determines the actual compacted state rather than assuming
that the ratio alone fixes it.  Grid spacing, interface width, event-domain
length, event density per physical GB length, excess-volume density, PF
mobility, time stepping, and every mechanism parameter remain unchanged.

The 384 by 384 B0 and `GTSC_GBTJ_Ks025` 100-200-step benchmark remains the
mandatory gate before long production.

The first 384 by 384 initialization pass also exposed a read-only scaling
cost: the solver evaluated full interfacial energy after every equilibration
step even though initialization does not consume that diagnostic, and the
long-run loop would retain 100,000 energy records.  Phase 1 now evaluates the
same diagnostic every 100 steps while leaving the PF update untouched.
Regression coverage proves bitwise-identical trajectories with dense and
sparse energy diagnostics.  Common-state equilibration omits the unused
energy calculation entirely.  The repaired full suite passes 174 tests.

The next 384 by 384 diagnostic reached its first compaction checkpoint at step
100 with 1038 active grains.  It then exposed that the legacy compiled PF
kernel still visited every pixel for every phase, despite compact local phase
support.  Phase 1 now constructs conservative row/column support masks that
include the exact one-cell cardinal halo and skips only phase rectangles that
cannot enter the existing stencil.  A regression compares the optimized path
with the same compiled kernel supplied all-true masks and requires bitwise
equal output for both periodic and non-periodic boundaries.  A 180-phase,
128 by 128 local benchmark improved from 0.0706 to 0.0158 seconds per step
(4.48-fold) with identical filling-constraint residuals.  The large-domain
preflight is restarted from the resulting clean SHA; no interrupted
initialization artifact is reused.

## Completed 384 by 384 production gate

The formal two-case preflight used source SHA
`f75e44c40f59b57174aa988c51f13a4e7b870fcf` and the common equilibrated state
`results/initial_conditions/seed-5101-5ba2f47025f7d22a.npz`.  Equilibration
finished at step 506 with exactly 800 active grains.  Both cases then completed
200 solver steps from that identical state:

| Case | Solver seconds/step | Peak RSS | Run size at 200 steps |
| --- | ---: | ---: | ---: |
| B0 | 1.019 | 5.13 GB | 22.7 MB |
| `GTSC_GBTJ_Ks025` | 7.260 | 8.05 GB | 349.4 MB |

Two concurrent worst-case workers require about 16.1 GB of resident memory,
comfortably below the 32 GB host capacity.  The coupled case's apparently
large 200-step directory is dominated by the rotating restart pair:
127.9 MB for `checkpoint.npz` and 208.2 MB for `checkpoint.json`.  Treating
that pair correctly as fixed rather than linearly accumulated gives a
conservative 100,000-step projection near 6.98 GB, not 174.7 GB.  B0 projects
to roughly 1.17 GB.  The production launcher now records fixed restart and
accumulating output separately.

The preflight supports two workers.  Ordinary scalar/frame output is reduced
from every 100 steps to every 200 steps to preserve disk headroom; the complete
event ledger, deterministic 10% unit-stride event traces, 500-step rotating
checkpoint, 1000-step movie cadence, and one-percent growth-progress frames
remain unchanged.  Mobile cases are expected to reach 100 grains before the
100,000-step safety ceiling, while censored cases will be evaluated from their
late-time kinetic trends before any extension is considered.
