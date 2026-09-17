#!/usr/bin/env python3
"""Generate and validate the minimal native Qiu-SI Numba compatibility patch."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np


DRIVER_MEMBER = "PF_Codes/Idealized_microstructure/[100]_4ref/Bicrystal-4reference-el-pf.py"
FUNCTIONS_MEMBER = "PF_Codes/Idealized_microstructure/[100]_4ref/functions_4ref_new.py"
PRISTINE_ARCHIVE_SHA256 = "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90"
PRISTINE_DRIVER_SHA256 = "03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17"
PRISTINE_FUNCTIONS_SHA256 = "3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf"
FAILING_EXPRESSION = "                        gb_unsort[qq] = np.concatenate(([TJ1, TJ2], gb_unsort[qq]))"
PATCHED_EXPRESSION = """                        # Numba typing compatibility repair: preserve the pristine
                        # [TJ1, TJ2, original rows...] order with one homogeneous array.
                        old_gb = gb_unsort[qq]
                        joined_gb = np.empty((old_gb.shape[0] + 2, old_gb.shape[1]), dtype=old_gb.dtype)
                        joined_gb[0, :] = TJ1
                        joined_gb[1, :] = TJ2
                        joined_gb[2:, :] = old_gb
                        gb_unsort[qq] = joined_gb"""


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def prepend_junction_rows_reference(boundary: np.ndarray, first: np.ndarray, last: np.ndarray) -> np.ndarray:
    """Independent non-JIT oracle for pristine ``[first,last,boundary...]`` order.

    Rows are copied without sorting or deduplication. Empty ``(0, n)`` boundary
    arrays are valid. The boundary dtype and trailing shape control the output.
    """
    boundary = np.asarray(boundary)
    first = np.asarray(first)
    last = np.asarray(last)
    if boundary.ndim != 2:
        raise ValueError("boundary must be two-dimensional")
    if first.shape != boundary.shape[1:] or last.shape != boundary.shape[1:]:
        raise ValueError("junction rows must match the boundary trailing shape")
    result = np.empty((boundary.shape[0] + 2, boundary.shape[1]), dtype=boundary.dtype)
    for column in range(boundary.shape[1]):
        result[0, column] = first[column]
        result[1, column] = last[column]
    for row in range(boundary.shape[0]):
        for column in range(boundary.shape[1]):
            result[row + 2, column] = boundary[row, column]
    return result


def prepend_junction_rows_compatibility(boundary: np.ndarray, first: np.ndarray, last: np.ndarray) -> np.ndarray:
    """Production operation inserted into the patched native function."""
    result = np.empty((boundary.shape[0] + 2, boundary.shape[1]), dtype=boundary.dtype)
    result[0, :] = first
    result[1, :] = last
    result[2:, :] = boundary
    return result


def patch_functions_source(pristine: str) -> str:
    if pristine.count(FAILING_EXPRESSION) != 1:
        raise ValueError("pristine failing expression is absent or ambiguous")
    return pristine.replace(FAILING_EXPRESSION, PATCHED_EXPRESSION)


def generate(archive: Path, output: Path) -> dict:
    archive_bytes = archive.read_bytes()
    if sha256_bytes(archive_bytes) != PRISTINE_ARCHIVE_SHA256:
        raise ValueError("native archive hash mismatch")
    with ZipFile(archive) as stream:
        driver = stream.read(DRIVER_MEMBER)
        functions = stream.read(FUNCTIONS_MEMBER)
    if sha256_bytes(driver) != PRISTINE_DRIVER_SHA256 or sha256_bytes(functions) != PRISTINE_FUNCTIONS_SHA256:
        raise ValueError("pristine member hash mismatch")
    patched = patch_functions_source(functions.decode()).encode()
    output.mkdir(parents=True, exist_ok=False)
    (output / "Bicrystal-4reference-el-pf.py").write_bytes(driver)
    (output / "functions_4ref_new.py").write_bytes(patched)
    diff = "".join(difflib.unified_diff(
        functions.decode().splitlines(keepends=True), patched.decode().splitlines(keepends=True),
        fromfile="pristine/functions_4ref_new.py", tofile="numba_compatibility_v1/functions_4ref_new.py",
    ))
    record = {
        "schema": "qiu-native-compatibility-patch-v1",
        "source_variant": "numba_compatibility_v1",
        "reason": "Numba 0.61 cannot type heterogeneous list/2-D-array np.concatenate input",
        "affected_function": "find_gb",
        "pristine": {"archive_sha256": PRISTINE_ARCHIVE_SHA256, "driver_sha256": PRISTINE_DRIVER_SHA256, "functions_sha256": PRISTINE_FUNCTIONS_SHA256},
        "patched": {"driver_sha256": sha256_bytes(driver), "functions_sha256": sha256_bytes(patched)},
        "unchanged_call_sites": ["Bicrystal-4reference-el-pf.py:158"],
        "behavior": {"row_order": "TJ first row, TJ last row, then existing boundary rows", "duplicates": "preserved", "sorting": "none", "dtype": "existing boundary dtype", "empty_boundary": "returns the two junction rows"},
        "exact_unified_diff": diff,
        "validation_fixtures": [],
        "environment_identities": [],
    }
    (output / "native_compatibility_patch.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(generate(args.archive, args.output), indent=2))


if __name__ == "__main__":
    main()
