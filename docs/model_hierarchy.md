# Model hierarchy

The code separates four layers that must not be conflated:

1. `pf`: constrained multiphase-field capillarity and optional mechanical driving.
2. `entities` and `encounters`: persistent grain/GB/TJ identities and geometric reaction-coordinate accumulation.
3. `disconnections`, `stochastic`, `climb`, and `mechanics`: mode admissibility, Arrhenius first passage, serial accommodation, and stored/nonlocal stress.
4. `analysis`: grain tracks, scaling, activation energy, Qiu metrics, and intermittency statistics.

Two compatibility models are supported: `explicit_modes` derives difficult states from finite mode feasibility, while `geometric_surrogate` creates them by cumulative hazard in GB-length change, TJ path, swept area, or slip. Two mechanics backends are supported: `local_memory` and `qiu_full_field` (a periodic FFT eigenstrain Green-operator implementation inspired by, but not copied from, Qiu’s line-disconnection solver).

Every stochastic clock uses cumulative hazard. Parallel alternative modes compete by summed rate; serial climb stages retain separate residence times.

