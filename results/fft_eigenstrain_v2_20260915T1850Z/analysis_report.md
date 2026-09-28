# FFT_EIGENSTRAIN_V2 post-simulation analysis

Classification: **FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE**.

Job `55932457` completed 8379 steps through physical time 335.16, reducing the population from 799 to 100 grains. Total energy decreased at every recorded step; the largest increment was -2.151e-01.

The full-horizon linear fit of population-radius squared has slope 4.032 and R-squared 0.9963. The largest one-step population loss was 3; in the latter half it was 2, with a maximum 100-step loss of 8. No disconnected grains were recorded. These diagnostics do not show the late non-self-similar collapse seen in the separate QIU_LEGACY_FORENSIC model.

This result is a single timestep calculation. The running dt/2 refinement remains a separate trajectory, so timestep convergence is not claimed here.
