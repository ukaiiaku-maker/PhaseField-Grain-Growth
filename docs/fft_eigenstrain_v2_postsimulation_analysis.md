# FFT_EIGENSTRAIN_V2 post-simulation analysis

The retrieved coarse-timestep seed-5101 trajectory (run
`20260911T003946Z-nogit-c53868`, Slurm job `55932457`) completed 8,379 steps
through physical time 335.16 and the preregistered 100-grain target. It is a
separate constitutive model and is not a native Qiu SI control.

The analysis classifies this trajectory as
`FFT_EIGENSTRAIN_V2_COMPLETED_NO_LATE_AVALANCHE`. Total energy decreased at
all 8,378 recorded increments. The maximum one-step population loss was three
grains; during the latter half it was two, and the maximum late 100-step loss
was eight. No disconnected grains were recorded. The full-interval linear fit
of population-radius squared gives slope 4.0325 and R-squared 0.9963.

The checksummed package in
`results/fft_eigenstrain_v2_20260915T1850Z` includes kinetic, energetic,
mechanical, topology, grain-distribution, and field-evolution figures in PNG
and PDF; full diagnostic and boundary CSV tables; a JSON summary; and an
artifact manifest tied to the verified result archive. The dt/2 calculation
is still running as a separate refinement trajectory, so this result does not
yet establish timestep convergence.
