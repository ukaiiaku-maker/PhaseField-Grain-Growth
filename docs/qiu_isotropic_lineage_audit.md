# Qiu isotropic lineage audit

The cited calculations do not form one exact continuation chain. Jobs
`55930486` and `55948258` are both `QIU_LEGACY_FORENSIC`. The second starts
from the exact step-9000 recovery and replays through step 9984 with added
source/work diagnostics. Its final numerical checkpoint arrays match the
first replay, but steps 9000–9984 overlap. It is therefore an exact diagnostic
replay, not a successor segment.

Jobs `55932457` and `55950433` are `FFT_EIGENSTRAIN_V2`. They start from the
same seed-5101 state but use accepted timesteps 0.04 and 0.02. They are base
and timestep-refinement calculations, not segments. Job `55932457` completed
at step 8379, physical time 335.16, and 100 grains; its verified retrieved
archive contains checkpoint SHA-256
`d190c50b68328629e3c28779df4e5dcc8f3f80670a54a507d74ee34ee3513665`.
Job `55950433` remains active and is not fetchable yet.

The sole native `[100]` four-reference plan
`20260912T180419Z-nogit-a25f0c` is still prepared and has no Slurm ID. Thus no
native result exists to promote. The authoritative isotropic baseline remains
the pristine `QIU_SI_REFERENCE`, while its current classification is
`QIU_ISOTROPIC_LINEAGE_INCOMPLETE`.

The machine audit records all available hashes and explicitly preserves these
separate identities. `QIU_LEGACY_FORENSIC` is classified duplicated because
of the controlled overlap. `FFT_EIGENSTRAIN_V2` is classified incompatible as
a stitch because its two trajectories differ in timestep. Neither is the
native Qiu SI baseline.
