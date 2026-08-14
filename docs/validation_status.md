# Validation status

This file is generated/updated by the validation campaign. A production label is forbidden unless all prerequisite checks named in the run manifest pass.

| Group | Status | Evidence |
|---|---|---|
| Source audit | complete | `source_manifest.md`, `external_sources.md`, `qiu_code_audit.md` |
| Intrinsic PF | pending run | `configs/validation/intrinsic.yaml` |
| Entity/geometry | pending run | unit and geometry tests |
| Stochastic engine | pending run | stochastic distribution tests |
| Mode selection/TJ | pending run | mode statistics and conservation tests |
| Shear/climb | pending run | mechanics and quota-conservation tests |
| Production scaling | blocked by validation | campaign manifests |
| Activation campaign | blocked by scaling validation | campaign manifests |

