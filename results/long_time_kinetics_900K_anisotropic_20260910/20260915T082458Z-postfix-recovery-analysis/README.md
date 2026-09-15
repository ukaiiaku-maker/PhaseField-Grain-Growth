# Anisotropic postfix compact analysis

Scientific source: `1316cc89dbabfb41cb883b0d4a4c74738cc2bef6`. Verified result archive: `b1d3eec6b0b693641958c78fcc73bf2bd7a83bc827eff861fdd661cfd116fc42`.

The three continuous trajectories are complete. A0 and coarse A2 restart exactly; the half-dt restart ends at step 126/128. The formal classification remains `A2_POSTFIX_OPERATIONALLY_INCOMPLETE`.

Coarse A2 energy decreases from 5236.012485284901 to 5149.634685219778; half-dt A2 decreases to 5115.519299781972. Both have zero positive accepted-step increments. The matched-horizon energy difference is 0.6669%.

The main operational result is low-amplitude support proliferation: the half-dt final field has mean support 116.8 phases per cell at 1e-14, but 1.62 at amplitude 0.01. This drives the late-stage cost increase.

Figures: [energy](energy_evolution.png), [convergence and morphology](convergence_morphology.png), [microstructures](microstructure_fields.png), [grain areas](grain_area_statistics.png), [support growth](support_growth.png), [amplitude support](amplitude_support.png), [interface orientation](interface_orientation.png), and [runtime scaling](runtime_scaling.png).

Machine-readable results are in `analysis_summary.json`; plotted data are in the CSV files; hashes are in `artifact_manifest.json`.
