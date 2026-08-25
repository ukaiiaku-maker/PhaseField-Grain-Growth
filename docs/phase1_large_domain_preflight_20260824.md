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
