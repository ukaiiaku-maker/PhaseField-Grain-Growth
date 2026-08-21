# Expanded jerky-mechanism analysis

## Decision

Do not launch the 36-run temperature/seed campaign yet.  The expanded analysis
strengthens the conclusion that C is the minimum one-mechanism surrogate for
GTSC at 900 K, but it also shows that the current GTSC intermittency is mainly a
serial-climb renewal process.  It does not show substantial systematic stress
buildup before release at the saved boundary-track cadence.

No new phase-field simulations were launched for this analysis.

## Event-burst causality

The operational large-burst definition is the strictly positive top 5% of
grain-radius-rate intervals.  A release is linked to every grain named in its
ledger row and is considered preceding over one saved-frame interval.  The null
randomly reassigns the same number of release-associated intervals within each
grain, using 200 deterministic shuffles.  This controls both grain lifetime and
between-grain differences in event propensity.

| Case | P(large) | P(large given release) | Actual risk ratio | Shuffled ratio | Excess | Top-5% release-fraction excess |
|---|---:|---:|---:|---:|---:|---:|
| G | 0.0500 | 0.0563 | 1.126 | 0.991 | 0.135 | 0.0155 |
| T | 0.0500 | 0.0571 | 1.143 | 1.087 | 0.055 | 0.0019 |
| C | 0.0500 | 0.0781 | 1.561 | 0.997 | 0.564 | 0.1125 |
| GT | 0.0500 | 0.0603 | 1.206 | 0.989 | 0.218 | 0.0280 |
| GTSC | 0.0500 | 0.0773 | 1.546 | 1.004 | 0.542 | 0.1121 |
| B0, S, QIU | 0.0500 | no releases | -- | -- | 0 | 0 |

The maximum release-leading growth-rate cross-correlation excess over shuffled
counts is 0.294 for G, 0.207 for T, 0.169 for C, and 0.128 for GTSC; it is zero
for B0, S, and QIU.  Thus CV and stationarity alone are not event-specific, but
C and GTSC have nearly identical, nonzero release-burst coupling.  G is causal
but weaker by the top-motion statistic, and T is only weakly resolved at this
output cadence.

The slow-climb discriminating case gives the largest causal risk-ratio excess,
1.156, versus 0.414 for reference climb and 0.233 for fast climb.  Longer serial
waiting therefore concentrates motion more tightly after quota completion.

## Pre-release work and effective barriers

For every reconstructable GB pin episode, saved curvature, resolved shear, and
free-volume deficit were converted to capillary, shear, and vacancy work using
the release event's activation volumes.  The exact activation-work row was then
appended at release.  Short episodes that begin and end between boundary-output
frames cannot be reconstructed.

The paired median effective-barrier change from the first saved pinned frame to
release is small:

| Case/event | Reconstructed events | Median wait | Median Delta effective barrier | 10th--90th percentile |
|---|---:|---:|---:|---:|
| G compatibility | 6,297 | 1.189 | -0.00069 eV | -0.00968 to +0.00508 eV |
| GT compatibility | 6,045 | 1.165 | -0.00064 eV | -0.00926 to +0.00458 eV |
| C climb completion | 16,562 | 1.554 | -0.00049 eV | -0.00874 to +0.00566 eV |
| GTSC compatibility | 1,751 | 1.296 | -0.00036 eV | -0.00716 to +0.00361 eV |
| GTSC climb completion | 15,162 | 1.532 | -0.00097 eV | -0.00915 to +0.00450 eV |

At 900 K, kBT is about 0.0776 eV.  These median changes are less than 1.3% of
kBT, and their distributions include both increasing and decreasing barriers.
The current saved data therefore do not support a picture in which most pins
steadily accumulate capillary/shear/vacancy work until release.  G is better
described by its persistent multihit clock, and C by its serial climb cycle.
The GS compatibility distribution is broader (-0.0380 to +0.0199 eV over its
10th--90th percentiles), so shear can still matter for a minority of G releases.

## Arrhenius work tails

Signed multipliers are exp(W/kBT); change fractions use exp(abs(W)/kBT), so
either acceleration or suppression counts.

- In GS, the shear multiplier is 1.011, 1.331, 1.563, and 2.399 at the 50th,
  90th, 95th, and 99th percentiles.  Shear changes 2.13% of releases by more
  than 2x, 0.219% by more than 5x, and 0.031% by more than 10x.
- In GTS, the corresponding shear multipliers are 1.018, 1.249, 1.451, and
  2.257; 1.61% exceed a 2x magnitude change.
- In GTSC, the shear multiplier reaches only 1.059 at p99 and no event changes
  by 2x.  Conditional on nonzero vacancy work, the C multiplier reaches 1.103
  at p99 and likewise never changes a rate by 2x.

The work tails confirm that S can modify G/T hazards when C is absent, but S and
vacancy chemical work are small hazard corrections in the present GTSC
parameterization.  C dominates through its serial gate, not through a large
vacancy-work term in the GB Arrhenius barrier.

## Full 2^4 factorial effects

These are descriptive high-minus-low contrasts from one realization, not
inferential estimates with seed-level uncertainty.

- Jerkiness CV: C +1.252, S +0.642, and SxC -0.597.
- Stationary fraction: C +0.0802, S +0.0676, and SxC -0.0651.
- Linear growth coefficient: C -0.02137, S -0.01868, and SxC +0.01846.
- Event Fano factor: C +10.84, G +6.87, and GxC -7.35.
- Top-5% release-fraction excess: C +0.0934, TxC -0.0112, and G +0.0094.

The C and S effects are strongly non-additive, exactly as suggested by the raw
table.  Once C imposes the slow serial gate, it masks most incremental S effects.

## Common-topology windows

C and GTSC remain close when compared at the same microstructural state rather
than at equal elapsed duration.

| Window | Case | dR/dt | Stationary fraction | Pinned fraction | Release Fano |
|---|---|---:|---:|---:|---:|
| N=190 to 175 | C | 0.007746 | 0.9490 | 0.5124 | 1.081 |
| N=190 to 175 | GTSC | 0.006831 | 0.9511 | 0.4920 | 1.110 |
| N=175 to 160 | C | 0.005183 | 0.9733 | 0.4838 | 0.977 |
| N=175 to 160 | GTSC | 0.004876 | 0.9737 | 0.4652 | 0.900 |
| N=190 to 160 | C | 0.005486 | 0.9634 | 0.4946 | 1.148 |
| N=190 to 160 | GTSC | 0.005241 | 0.9630 | 0.4756 | 1.118 |

This rules out unequal endpoint topology as the explanation for C approximately
equaling GTSC at 900 K.

## Waiting-time survival and hazard

Fits include right censoring.  Exponential has one fitted parameter; gamma and
Weibull have two, and models are ranked by AIC.

- G_LOW, G_REF, and G_HIGH strongly prefer gamma waiting times, with shapes
  2.82, 2.43, and 2.21.  Their exponential delta-AIC values are 1,227, 978,
  and 659.  This is consistent with persistent multihit kinetics rather than a
  single Poisson barrier.
- C_FAST prefers gamma with shape 1.22.  C_REF and C_SLOW prefer Weibull with
  shapes 0.934 and 0.905; exponential delta-AIC is 60 and 46.  The decreasing
  effective hazard is more consistent with mixed domain populations/history
  than a homogeneous serial Erlang clock.
- TJ waiting was not tracked from encounter.  Lower-bound fits from the first
  persisted TJ hit to release reject an exponential for T_REF and T_HIGH by
  delta-AIC about 28 and 27, but gamma and Weibull are close.  T_LOW contains
  many same-step hit/release sequences, so its fitted sub-frame shape is not a
  reliable physical lifetime.

## Mechanism occupancy reconstruction

Only a generic GB `blocked` flag was persisted.  Completion events identify G
or C involvement; episodes with neither completion remain unresolved.  When
both completions occur, their times split multiple- and single-gate portions.
This is event-classified blocked-domain time, not an exact internal-state trace.

For GTSC, recorded GB blocked-domain time is reconstructed as:

- C-limited: 85.87%
- G-limited: 5.21%
- simultaneous G and C: 1.97%
- unresolved: 6.95%

GTC is almost identical: 85.45% C-limited, 5.18% G-limited, 2.15% multiple,
and 7.22% unresolved.  GTSC also has 291 resolved TJ cycles and 394.4 units of
post-first-hit TJ waiting time, but the absence of time-resolved TJ tracks means
that quantity cannot be placed in the same exact occupancy denominator.

The defensible conclusion is that C controls most persisted GTSC GB arrest at
this parameterization; G and T are active but secondary, while S modifies
hazards rather than forming a separate gate.

## Formal minimum-model ranking

The standardized observable vector uses linear K, stationary fraction, top-5%
motion concentration, event Fano, median pin duration, and top-5% release
coupling.  Against GTSC, C has RMS discrepancy 0.133 and is the Pareto-optimal
one-mechanism model.  SC improves this to 0.095 with two mechanisms, and TSC to
0.067 with three.  G, T, and S are all dominated among one-mechanism models.

Thus C is the current minimum model for reproducing GTSC observables, while SC
and TSC quantify the modest improvements available at added complexity.

## Observability limits and next decision gate

The existing files cannot supply exact four-way time occupancy because they do
not persist reason-specific G/T/C pending flags or TJ tracks.  They also cannot
resolve stress accumulation in pin episodes shorter than the four-step boundary
output cadence.  Spatial/arclength branching analysis was not included in this
highest-priority pass.

Before a temperature ensemble, persist reason-specific gate states, TJ blocked
state, instantaneous component work/rate at ordinary output frames, and stable
spatial/arclength event coordinates.  The smallest validation is then two
diagnostic replays, C and GTSC at 900 K/seed 5101, to verify exact occupancy and
event-triggered barrier histories without changing physics.

If those replays confirm the reconstruction, replace the initial 36-run launch
with a staged extreme-temperature screen: B0/C/GTSC at the lowest and highest
planned temperatures and three paired seeds (18 runs).  Fill the two interior
temperatures only if C continues to track GTSC or if an emerging divergence
needs localization.

## Provenance

- Factorial raw data:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c`
- Discriminating raw data:
  `results/jerky_discriminating_900K/20260821T080200Z-124f03cd0d`
- Expanded outputs:
  `results/production_summaries/jerky_mechanism_expanded/`
- Analysis implementation: `08e47a0d74a98198740fb2291146d9dea9b06139`
- Simulation implementation: factorial `027b91c`; discriminating `2265d32`
- Branch: `fix/jerky-mechanism-integrity-20260820`
- Validation: 124 tests passed on 2026-08-21 using `/opt/anaconda3/bin/python`
- Factorial bookkeeping: 16 G/T/S/C subsets including B0, plus QIU, equals
  17 completed cases; there was no missing eighteenth case.

