# Validation status

This file is generated/updated by the validation campaign. A production label is forbidden unless all prerequisite checks named in the run manifest pass.

| Group | Status | Evidence |
|---|---|---|
| Source audit | complete | `source_manifest.md`, `external_sources.md`, `qiu_code_audit.md` |
| Intrinsic PF | fast validation passed | planar stationarity; circle slope; mobility; energy; restart in `tests/pf` |
| Entity/geometry | fast validation passed | entity persistence/retirement and TJ path tests |
| Stochastic engine | fast validation passed | exponential, Erlang, packet-reset, time-dependent and parallel tests |
| Mode selection/TJ | fast validation passed | isotropy, discrete minimum Burgers, feasible-combination and TJ persistence tests |
| Shear/climb | fast validation passed | sign, balance, exchange, transport and serial-time tests |
| Production scaling | blocked by validation | campaign manifests |
| Activation campaign | blocked by scaling validation | campaign manifests |

Fast-suite result at the v0.7 milestone: **20 passed in 13.56 s** on Python 3.11.3, NumPy 2.4.6, SciPy 1.17.1. These are prerequisite component validations, not a claim that mesh/ensemble production campaigns have completed.
