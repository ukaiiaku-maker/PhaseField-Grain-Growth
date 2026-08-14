# Validation status

This file is generated/updated by the validation campaign. A production label is forbidden unless all prerequisite checks named in the run manifest pass.

| Group | Status | Evidence |
|---|---|---|
| Source audit | complete | `source_manifest.md`, `external_sources.md`, `qiu_code_audit.md` |
| Intrinsic PF | quantitative validation passed | `results/validation/numerical_validation.json` |
| Entity/geometry | fast validation passed | entity persistence/retirement and TJ path tests |
| Stochastic engine | fast validation passed | exponential, Erlang, packet-reset, time-dependent and parallel tests |
| Mode selection/TJ | fast validation passed | isotropy, discrete minimum Burgers, feasible-combination and TJ persistence tests |
| Shear/climb | fast validation passed | sign, balance, exchange, transport and serial-time tests |
| Production scaling | blocked by validation | campaign manifests |
| Activation campaign | blocked by scaling validation | campaign manifests |

Fast-suite result at the v0.7 milestone: **20 passed in 13.56 s** on Python 3.11.3, NumPy 2.4.6, SciPy 1.17.1. The subsequent integrated engine/elastic-backend additions raise this to **22 passed in 13.90 s**. A 33-regime, three-step matrix completed only as a wiring smoke test; its near-zero fitted coefficients are deliberately not scientific scaling results. These are prerequisite component validations, not a claim that mesh/ensemble production campaigns have completed.

Quantitative intrinsic campaign at commit `7afc7d40e68fad2a726479d84bc7f485db04ba75` passed: six circle cases over `dx={1,0.75,0.5}` and `M={0.1,0.2}` had `R²>0.99999`; sharp-interface slope errors were 0.27–1.18%; doubling mobility changed the slope by factors 2.013–2.015; mesh slope spread was 0.247%; interfacial energy decreased monotonically; measured grid surface-energy anisotropy was 2.60%; the equal-energy TJ angles were 121.56°, 115.71°, and 122.72° (maximum 120° error 4.29°).
