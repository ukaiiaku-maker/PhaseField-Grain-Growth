# Source manifest

The files below were present at campaign startup on 2026-08-14. They are immutable inputs and are excluded from Git; this manifest makes the exact inputs auditable.

| File | SHA256 | Represents | Equations/implementation used |
|---|---|---|---|
| `PF_Codes.zip` | `2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90` | Qiu et al., *Why grain growth is not curvature flow*, Zenodo 15120372 | Qiu/Salvalaglio multi-order-parameter PF, two/four/six disconnection-reference shear coupling, line-disconnection stress construction; audit in `qiu_code_audit.md` |
| `README.md` | `3240f08fa64bd15757caaafc903b47da4be5f0334c3710da46d3aaaf758d5c86` | Qiu archive description | DOI/provenance only; left unchanged |
| `GrainGrowth_V6.docx` | `7468b3e8a8aab4354e90edb696b9a698fd4ad0aa71b7f0be4a8e563650f11d15` | *Non-Steady-State Grain Growth: A Hazard-Limited Framework* draft | Class A/B/C/D closures; geometric hazard classes; series/parallel composition; exact equations are transcribed in `physics_equations.md` and `analysis/analytical_models.py` |
| `Supporting_Information.docx` | `e4a52f9bbd530f2d39d3fc1bf6f15a20370fb2535fb8fe3086bdc962404d1ea7` | Supporting information to the hazard-limited framework | Lambert-W, Poisson multihit, Zener, vacancy, TJ, exchange and diffusion-transient derivations |
| `Archive.zip` | `8e6af6b7bdff6272a6de30662692a18d84a61220c98479296e509f4ddbdae76de` | Snapshot of the supplied manuscript and predecessor code collection | Redundant preservation copy; contents individually inventoried below |
| `agg_copula_map_demo_SI_MARGINALS_3D_COPULA_RTERM_FULLJKO_v2.m` | `c1e6ec6337421c5918d949c2939e101ca39c2ed81f5f19e87363b93e94067c07` | Sintering/nucleation-renewal aggregate model | Renewal statistics and correlated marginal sampling used as conceptual checks only |
| `classB_multihit_jerkiness_map.m` | `a7c3f1119f122c8e7f448a2c5808d196a479e1a8632f0cb7242d094cb64c24cd` | Work-limited multihit parameter sweep | Poisson-tail closure and jerkiness map regression targets |
| `sweep_classB_multihit_v5_mechanisms.m` | `6dd20b90ce38a296d6e30662692a18d84a61220c98479296e509f4ddbdae76de` | Class-B mechanism sweep | Mechanism-specific size scalings and work gating |
| `plot_classB_worklimited_multihit_v5.m` | `c5b38ad5d414c2afb18e6626f73f374251452d134cb01efb9f0fe8833988951c` | Class-B analysis/plotting | Plot definitions and diagnostic quantities |
| `graingrowth_polycrystal_TJ_hazard_gate_RENEWAL_patched_v4 copy.txt` | `6ecbd654201aceeba34fc763be2832ff991594cbb314239b252780611b222e4e` | Earlier PF/TJ hazard-gate prototype | Pixel-threshold pathology reference; entity-based implementation replaces pixel state |
| `potts_grain_growth_GS_dependent_avalanche.m` | `7e60e757bdf166bdfb3ec0a29ac0fdd075fbd94bd2ce26ea0bf99b1e69c8299d` | Earlier Potts avalanche model | Legacy random-obstacle and grain-size-dependent behavior as external controls only |
| `codefromShen_merged.ipynb` | `4d1138097840dfb197509e6082d01c0cbab6d07979936346aeab46cd35bb7910` | Previous PF notebook/archive | Numerical and initialization comparison; not imported as package code |
| `Instructions.md` | `a86ee1782d499c839b99627f7f398dfa2edacc918308c01ef7f23e3c745a7aa5` | Authoritative campaign directive | Requirements and acceptance criteria |

The `PF_Codes.zip` MD5 is `6cd49ca72eba89210abb96e700342f12`, exactly matching the Zenodo record. Original files are never used as mutable run outputs.

