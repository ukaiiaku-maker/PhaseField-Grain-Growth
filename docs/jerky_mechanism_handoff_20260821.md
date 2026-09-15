# Jerky grain-growth mechanism handoff (2026-08-21)

## Decision

At the present common parameter set, **C (serial climb/free-volume gating) is the
provisional minimum model that reproduces the long-horizon GTSC response**.  C
and GTSC differ by 1.9% in end grain count, 1.7% in jerkiness CV, less than
0.1% in effective stationary fraction, 3.8% in pinned fraction, and 0% in
median pin duration.  Their Fano factors differ by 13.4%.

This is not yet a manuscript-level conclusion that climb alone is sufficient.
The discriminating screen shows that C is strongly activation-scale controlled,
and its mean chemical activation work is much smaller than capillary work.  The
current hierarchy is therefore dominated by the serial climb gate rather than
by a large chemical transition-state correction.  G/T/S effects can be masked
under that dominant C scale.

## Integrity repairs

- Triple-junction releases now use signed capillary, shear, and vacancy driving
  work in their effective barriers and write `tj_compatibility_release` rows to
  `activation_work.csv`.
- The corrected TJ gate ignores the legacy pixel radius and is applied once in
  physical length units.
- A separate `compatibility_pending` state prevents climb-only arrests from
  falling through the GB compatibility-release path.  C now contains only climb
  completions; it no longer secretly contains G.
- Gate-only climb completion still accommodates its quota once.
- Frame and analysis artifacts now report boundary shear RMS, stored shear
  energy, nonzero GB coverage, free-volume coverage, TJ-specific work, effective
  stationarity, motion concentration, pin durations, and shear-release audits.

## Main mechanism results at 900 K

| case | N end | CV | stationary | Fano | pinned | median pin | interpretation |
|---|---:|---:|---:|---:|---:|---:|---|
| B0 | 100 | 1.191 | 0.751 | — | 0 | — | smooth reference |
| G | 100 | 1.429 | 0.837 | 13.59 | 0.126 | 1.28 | robust GB-local event intermittency |
| T | 100 | 1.329 | 0.808 | 5.25 | 0 | — | TJ waiting control, less burst-count overdispersion |
| GT | 100 | 1.522 | 0.862 | 14.23 | 0.117 | 1.28 | slower/more stationary than G or T; no decisive series-fit advantage |
| GS | 120 | 2.042 | 0.921 | 20.83 | 0.052 | 1.28 | shear materially modifies G hazards |
| GTS | 123 | 2.125 | 0.922 | 16.25 | 0.049 | 1.28 | T adds little beyond GS at this scale |
| C | 159 | 3.210 | 0.963 | 20.99 | 0.497 | 1.60 | dominant serial climb gate |
| GTC | 158 | 3.209 | 0.958 | 19.61 | 0.495 | 1.60 | close to C |
| GTSC | 162 | 3.267 | 0.963 | 18.51 | 0.479 | 1.60 | full reference; close to C/GTC |
| QIU | 150 | 3.127 | 0.956 | — | 0 | — | continuous-mechanics reference, not event limited |

Exact fitted exponents were not used as acceptance criteria.  G, T, and GT all
preferred the simple linear fit over the selected finite window; T was faster
and less Fano-intermittent than G.  GT's series fit was informative but not
better than the simple linear fit, so this one-seed matrix does not establish a
clean crossover.

## Work and internal-state interpretation

- In GS, TS, and GTS, mean absolute shear work is about 0.009--0.011 versus
  0.048--0.051 capillary work (roughly 18--23%): shear materially changes the
  hazard when C is absent.
- In GTSC, mean absolute shear work falls to 0.00116 versus 0.04665 capillary
  work (about 2.5%).  C suppresses the shear contribution enough that S is not
  required to reproduce the present full-model response.
- Overall free-volume work is about 0.00036 in GC/GTC/GTSC (less than 1% of
  mean capillary work); its effect is concentrated in genuine GB compatibility
  events.  C's dominant effect is serial gating, not large average chemical work.
- All smoke S cases had nonzero evolving shear, finite stored energy, and
  approximately 99--100% nonzero GB coverage.  Non-S cases correctly had zero
  stored shear.
- With the reference packet size 2, GB events zeroed the pre-event shear state
  in 99.8--100% of releases.  Packet sizes 1 and 0.5 reduced this to 96.8% and
  80.4% without materially changing the short GS kinetics.  The near-total
  reset is therefore a coarse packet-scale artifact, not a required physical
  outcome.

## Discriminating screens

- Raising G barriers increased stationary fraction from 0.768 to 0.881,
  pinned fraction from 0.065 to 0.217, and median pin duration from 0.48 to
  2.88.  G remains a robust intermittent activated skeleton.
- Raising T barriers increased stationary fraction from 0.746 to 0.852, but T
  retained much smaller event-count Fano factors (1.79--3.25) than G.
- GT followed the same monotonic arrest trend, with stationary fraction
  increasing from 0.760 to 0.911.
- Changing climb barriers from fast to reference to slow increased stationary
  fraction 0.896 -> 0.945 -> 0.952, pinned fraction 0.354 -> 0.516 -> 0.626,
  median pin duration 0.48 -> 1.44 -> 4.64, and Fano 35.7 -> 23.7 -> 43.4.
  The apparent minimum model is consequently sensitive to the climb scale.

## Morphology and resolution status

All 12 corrected smoke movies rendered successfully.  At matched grain count
near 175, spectral roughness ratios were within 2.1% of B0 and high-frequency
boundary power within 1.2%; there was no recurrence of progressive large-scale
waviness.  This supports the connected physical-arclength/local-pinning
architecture.

The long factorial did not save frame archives, so its morphology conclusion is
limited to boundary-track statistics.  Grid sensitivity of pin duration,
nonzero coverage, and local bowing remains unresolved.  Qiu's long-horizon
stress also grows much larger than in smoke and should remain a reference until
its quantitative reproduction and resolution behavior are checked separately.

## Recommended smallest next statistical campaign

Run **B0, C, and GTSC only**, at **four temperatures** (for example 800, 900,
1000, 1100 K) and **three matched seeds**: 36 runs total.  Use physical TJ and
GB correlation lengths, `gate_only`, the reference climb scale, and a reduced
packet size of 0.5 as the provisional shear-relaxation scale.  Add one 900 K
grid-refinement triplet for C and GTSC before treating pin-duration or morphology
statistics as resolution independent.

This is the smallest campaign that can test whether C reproduces GTSC across
temperature, estimate seed variability, and preserve a B0 kinetic baseline.
Do not launch the larger manuscript ensemble until this comparison either
confirms C or shows a temperature-dependent need for G/T/S.

## Provenance

- Corrected smoke:
  `results/jerky_integrity_smoke/20260820T220630Z-fffb17dfa0`
- Full factorial:
  `results/jerky_factorial_900K/20260820T233721Z-f5d010440c`
  (16 G/T/S/C subsets including B0, plus QIU: 17 completed cases total)
- Discriminating screen:
  `results/jerky_discriminating_900K/20260821T080200Z-124f03cd0d`
- Aggregate tables:
  `results/production_summaries/jerky_mechanism_integrity_smoke.csv`,
  `results/production_summaries/jerky_mechanism_factorial_900K.csv`, and
  `results/production_summaries/jerky_mechanism_discriminating_900K.csv`
- Smoke/factorial simulation code: `027b91c`
- Discriminating-screen simulation code: `2265d32`
- Branch: `fix/jerky-mechanism-integrity-20260820`
