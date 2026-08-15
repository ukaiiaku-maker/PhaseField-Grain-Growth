# Validation status

This file is generated/updated by the validation campaign. A production label is forbidden unless all prerequisite checks named in the run manifest pass.

| Group | Status | Evidence |
|---|---|---|
| Source audit | complete | `source_manifest.md`, `external_sources.md`, `qiu_code_audit.md` |
| Intrinsic PF | quantitative validation passed | `results/validation/numerical_validation.json` |
| Entity/geometry | fast validation passed | entity persistence/retirement and TJ path tests |
| Stochastic engine | quantitative validation passed | `results/validation/stochastic_validation.json` |
| Mode selection/TJ | fast validation passed | isotropy, discrete minimum Burgers, feasible-combination and TJ persistence tests |
| Shear/climb | fast validation passed | sign, balance, exchange, transport and serial-time tests |
| Production scaling | blocked by validation | campaign manifests |
| Activation campaign | blocked by scaling validation | campaign manifests |

Fast-suite result at the v0.7 milestone: **20 passed in 13.56 s** on Python 3.11.3, NumPy 2.4.6, SciPy 1.17.1. The subsequent integrated engine/elastic-backend additions raise this to **22 passed in 13.90 s**. A 33-regime, three-step matrix completed only as a wiring smoke test; its near-zero fitted coefficients are deliberately not scientific scaling results. These are prerequisite component validations, not a claim that mesh/ensemble production campaigns have completed.

Quantitative intrinsic campaign at commit `7afc7d40e68fad2a726479d84bc7f485db04ba75` passed: six circle cases over `dx={1,0.75,0.5}` and `M={0.1,0.2}` had `R²>0.99999`; sharp-interface slope errors were 0.27–1.18%; doubling mobility changed the slope by factors 2.013–2.015; mesh slope spread was 0.247%; interfacial energy decreased monotonically; measured grid surface-energy anisotropy was 2.60%; the equal-energy TJ angles were 121.56°, 115.71°, and 122.72° (maximum 120° error 4.29°).

Stochastic campaign at commit `7ba2e13b01e2779634793929dbe5efcdeed164d9` passed. The single-hit mean was 0.40366 for expected 0.4 (KS p=0.824), paired event times were invariant over an 8x timestep range to `3.8e-14`, the K=5 mean/CV were 1.9980/0.44658 versus 2.0/0.44721, packet completion was 0.37649 versus exact 0.37729, and parallel/serial means matched 1/7 and 0.7. Four-temperature mechanism-isolation fits recovered an input 0.65 eV barrier as 0.6516±0.0004 eV (single hit) and 0.6490±0.0013 eV (persistent K=5).

The first 65-run scaling pilot is preserved as a failed scientific result in `pilot_scaling_summary_failed_extinction.csv`. It exposed non-monotone grain count from regrowth of diffuse tails after nominal extinction. Production interpretation was rejected; the solver now permanently deactivates an order parameter once its maximum falls below the configured extinction threshold, and `test_extinct_grains_cannot_resurrect` guards the correction.

The five-seed corrected rerun confirmed monotone grain count (50 to 26–28 grains) and exact restart with the active mask. Its free-exponent result is not accepted: a 50-grain system enters finite-statistics slowdown before `n` and the post-transient window are jointly identifiable. The unreliable analysis is preserved; the next discriminating case uses 200 initial grains.

The subsequent five-seed, 200-order-parameter run also had monotone topology (mean 198.2 to 112.2 grains), but revealed a separate protocol error: the initial diffuse Voronoi relaxation was still included in physical time. Its two-rate trajectory is preserved and rejected for exponent inference. Equilibration now runs before time zero with mechanics/events disabled; time, hazard, and output clocks are then reset and the equilibrated state is written as `t=0`.
