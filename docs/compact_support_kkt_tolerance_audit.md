# Compact-support KKT-tolerance audit v2

Decision: **`COMPACT_SUPPORT_GATE_PACKAGING_MISALIGNMENT`**  
Superseding Gate-1 result: **`COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED_V2`**

The historical package and its classifications remain unchanged:

- packaged: `COMPACT_SUPPORT_OPERATOR_KKT_FAILURE`
- upstream: `COMPACT_SUPPORT_OPERATOR_TOLERANCE_NONCONVERGED`

Those labels resulted from a comparison that did not hold accepted step,
physical time, or parent lineage fixed. They do not describe a KKT residual,
feasibility, complementarity, termination, or determinism failure.

## Endpoint alignment

The exact common time is `0.2505237792079889`. The aligned states are coarse
step 3936 and fine step 7872. All fine checkpoints already existed at step
7872. The three missing coarse states were reconstructed from the nearest
earlier authoritative checkpoints using only 86, 128, and 59 accepted steps.
No trajectory was restarted from step zero and no phase-field interpolation
was used.

At fine dt, the `1e-8`, `1e-10`, and `1e-12` fields are bitwise identical.
At coarse dt, `1e-10` and `1e-12` are bitwise identical. The coarse `1e-8`
field differs from the tighter field by RMS `5.0417530e-17` and maximum
absolute difference `2.3425706e-14`, with zero dominant-label disagreement,
zero non-tie disagreement, and identical grain count and morphology scalars.
That case alone resumed a checkpoint produced on HPC3; the tighter siblings
started locally from the deterministic initial field. It is therefore not a
matched tolerance comparison.

The old fine mismatch was entirely an endpoint error: the `1e-10` final field
was taken at step 7936, while `1e-8` and `1e-12` stopped at step 7872. At the
common step 7872 all three are identical.

## Determinism and convergence

Two independent one-thread local processes started from the same exact
fine-dt, `1e-10`, step-7744 checkpoint and ran 128 accepted steps. Their
checkpoint SHA-256 values are identical, their fields are bitwise identical,
and their energy histories, timestep histories, support entries, and support
retirements match exactly. Both reproduce the historical step-7872 field hash
`03afcef0198aa0011dcbc81723770c11ffc1dd67553598185aa907b7f09aa91b`.

For every requested continuous field and scalar metric, the `1e-10` versus
`1e-12` distance is zero at both timesteps. Thus the tighter-distance
contraction criterion passes with a zero replay noise floor. All original
energy, KKT residual, feasibility, phase-sum, support, topology, and scaling
checks remain passing. A `1e-14` reference is unnecessary.

The first-divergence audit finds no tolerance-pair divergence for the fine
family or the tighter coarse pair. The loose coarse branch is already on a
different platform-produced parent by step 1024; the historical cadence did
not preserve the immediately preceding common platform state. The difference
is only roundoff scale and causes no active-set or topology difference.

Source inspection also shows that `anisotropic_kkt_tolerance` is validated
against the allowed values but is not consumed by the exact simplex projection.
It cannot currently alter support admission, retirement, candidate ordering,
or topology. The checkpoint-v1 schema did not store multiplier arrays, so the
audit records that limitation rather than inventing endpoint multiplier data;
the preserved scalar projection residuals remain at or below `4.44e-16`.

## Production choice

The selected production tolerance is `1e-10`. It is bitwise identical to
`1e-12` at both matched endpoints. The coarse `1e-12` case had materially
higher measured cost, while `1e-10` already gives the converged state. The
next gate is the intentional 128-step local/HPC platform pair from one exact
checksummed checkpoint. Native and anisotropic Qiu work remains blocked until
that platform comparison passes.

The complete immutable audit package is
`results/compact_support_kkt_tolerance_audit_20260917T013929Z`.
