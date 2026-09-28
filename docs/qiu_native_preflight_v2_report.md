# Native Qiu-SI preflight V2 report

## Decision

The V2 compiler gate is classified
`QIU_SI_NATIVE_ENVIRONMENT_INCOMPATIBLE`. The native baseline is not eligible
for release. The prepared full-baseline run
`20260912T180419Z-nogit-a25f0c` remains unsubmitted, and A0 correspondence and
anisotropic Qiu controls remain blocked.

## Executed identity

- Run: `20260917T041831Z-nogit-8512a3`
- Slurm job: `56102284`
- Terminal state: `COMPLETED`, exit `0:0`, wall time `00:01:21`
- Source variant: `pristine_archive`
- Archive SHA-256: `2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90`
- Driver SHA-256: `03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17`
- Function-module SHA-256: `3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf`
- Environment: Anaconda 2022.05, Python 3.9.12, NumPy 1.21.5,
  Numba 0.55.1, llvmlite 0.38.0

The run was a bounded compiler gate and was marked `baseline_eligible=false`
before submission. It used the exact native 500x500, 17-grain initialization.

## Numerical evidence

The native path produced finite, checksummed checkpoints at accepted steps 0
and 1. Phase sums remain normalized within `2.22e-16`, phase values remain in
`[0,1]`, and all recorded force and morphology diagnostics are finite.

| Accepted step | Physical time | Grain count | Boundary density | Mean active support | Maximum support | Checkpoint SHA-256 |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 0.0 | 17 | 0.027008 | 1.123520 | 3 | `eb60b24a2f120a079d28755885a20a1316cf8e6cad9cf6c21d14a2053063c882` |
| 1 | 0.1 | 17 | 0.025152 | 1.081024 | 3 | `278381ff6b66cdd95f1dca7e77d2d7d64ec502547aa0d3b67fdc52d8c3d9ce7a` |

Compilation evidence shows nopython signatures for the other instrumented
native kernels, while `find_gb` and `sort_gb` have none. The historical stack
therefore executes these functions through deprecated object-mode fallback.
The mandatory no-object-mode condition fails, so checkpoints 10, 100, and
1,000 and the step-100 restart check were intentionally not run.

The separately qualified `numba_compatibility_v1` array repair is exactly
equivalent on the actual initialization and downstream quantities. Modern
Numba then exposes additional incompatibilities in unchanged sorter logic.
Repairing those constructs would exceed the authorized array-assembly-only
scope, leaving no permitted source/environment combination for production.

## Retrieval and integrity

The result was fetched twice. HPC3 runner verification was stable both times.
The result archive SHA-256 is
`1a3c03d4f1905579d273207c3f94dce300b4b07552c7c29083a59b693ef8f413`;
the runner checksum-map identity stored by the orchestrator is
`d10cafbb7a99d1e7801b7d6f8f7b66c3699e0689eafc5026591d785a41de9dc8`.
Every artifact listed in `artifact_manifest.json`, including both checkpoints,
matches its recorded SHA-256.

The auxiliary `all-files.sha256` has a stale self-entry: it records the empty
file hash for itself before its contents were written. This is a packaging
defect only. That file is outside `artifact_manifest.json`; its defect does not
alter the verified runner archive, scientific artifacts, or classification.

Curated evidence is stored in
`results/qiu_native_preflight_v2_20260917T041831Z/`.
