# Native Qiu-SI execution compatibility repair

## Preserved source

The source archive `PF_Codes.zip` remains immutable at SHA-256
`2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90`.
The pristine four-reference driver and function module hashes are
`03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17`
and `3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf`.
The patched source is generated separately as source variant
`numba_compatibility_v1`; the archive is never modified.

## Pristine failure audit

`find_gb` is decorated with `@nb.jit(nopython=False, parallel=True)`. The
archive contains no environment lock or package-version declaration. All
ordinary files have the same 2025-05-26 ZIP timestamp. The only dependency
documentation says that Python Numba is required. The source relies on the
former object-mode-fallback meaning of `nopython=False`, which is evidence for
a Numba release before fallback removal. Two bounded historical candidates are
therefore tested: the unmodified Anaconda 2022.05 and 2021.11 module stacks.

At approximately line 118, pristine code constructs a tuple containing a
Python list of two one-dimensional triple-junction coordinate arrays and a
two-dimensional boundary-coordinate array, then passes it to
`np.concatenate`. NumPy coerces the list to a homogeneous two-row array.
Numba 0.61.0 instead types it as a list plus a 2-D array and rejects the
heterogeneous tuple during nopython compilation. No accepted step occurs.

The exact native initialization supplies `float64` phase and stress arrays,
integer support/count arrays, integer coordinate arrays, a 17-entry
misorientation array, and overlap lists containing ordered integer pairs and
triples. `phi` is `(17,500,500)`, `nf` is `(500,500)`, `mf` is
`(1000,500,500)`, and each stress array is `(500,500)`. A checksummed fixture
captures these arguments immediately before the first `find_gb` call.

`find_gb` does not mutate `phi`, `nf`, `mf`, orientations, or overlap lists. It
does mutate the supplied line-stress arrays and, at the configured stress
cadence, bulk-stress arrays. It returns two lists of variable-length integer
coordinate arrays. Row order is scientific state: later tangent construction,
reference selection, line-density integration, stress evaluation, and stress
extension traverse these rows in order. Duplicate rows are intentionally
retained until the unchanged native boundary sorter applies its own behavior.

## Independent reference and patch

The independent non-JIT oracle explicitly allocates the result and copies the
first triple-junction row, last triple-junction row, and original boundary rows
in that order. It performs no sort or deduplication, retains the boundary
dtype, supports an empty `(0,n)` boundary, and copies periodic-edge coordinates
as ordinary values.

The production repair applies the same two-pass homogeneous allocation inside
`find_gb`. It changes no criterion, periodic rule, ordering, beta/reference
logic, stress calculation, barrier, or phase-field update. The driver and its
only call site are byte-identical to the archive. Synthetic empty, singleton,
unequal-length, duplicate, and periodic-edge fixtures match the independent
reference exactly, including in Numba nopython mode.

## Environment and full-fixture result

The two bounded historical tests both execute the pristine initialization and
ten accepted steps. `anaconda/2022.05` provides Python 3.9.12, NumPy 1.21.5,
Numba 0.55.1, and llvmlite 0.38.0; the retained smoke completes in 41.39 seconds.
`anaconda/2021.11` provides Python 3.9.7, NumPy 1.20.3, Numba 0.54.1, and
llvmlite 0.37.0; the retained smoke completes in 40.71 seconds. Both compile `find_gb`
by deprecated object-mode fallback. They are historically runnable but fail
the preflight's mandatory no-object-mode condition.

The actual 500x500 initialization fixture has SHA-256
`bf39aa567b85ce0128520e4ad590a5af16eb55070a23c941029f68acf3b5b507`.
It contains 64 qualifying triple-junction prepend operations. Every operation
matches the independent oracle exactly. Running pristine and patched
`find_gb` in the same 2022.05 environment produces exact equality for the GB
coordinates and IDs, selected references, beta values, reference-resolved
line density, all three line and bulk stress components, capillary and barrier
terms, and signed elastic pair values.

The final `numba_compatibility_v1` function-module SHA-256 is
`1d656039b8dca20bac8f056ad195fdc201775331e73e0f6f5da1e866622f1f78`.
In the qualified modern environment (Python 3.13.5, NumPy 2.1.3, Numba
0.61.0, llvmlite 0.44.0), the repaired array assembly compiles in nopython
mode. Full `find_gb` compilation then stops in unchanged adjacent sorter code:
first at the undecorated `clockorcounterclock` call, and, when that call is
diagnostically annotated, at the sorter's built-in `abs` call on an array.
Those changes lie outside the authorized array-assembly-only patch.

The array repair is therefore proven equivalent, but neither permitted source
variant has a no-object-mode execution environment. The compatibility result
is `QIU_SI_NATIVE_ENVIRONMENT_INCOMPATIBLE`; no production source variant is
selected.
