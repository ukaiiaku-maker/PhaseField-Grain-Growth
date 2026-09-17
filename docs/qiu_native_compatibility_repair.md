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
