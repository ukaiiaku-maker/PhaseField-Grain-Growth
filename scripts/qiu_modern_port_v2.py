#!/usr/bin/env python3
"""Generate the minimal modern-Numba Qiu compatibility port from V1."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import shutil
from pathlib import Path


V1_FUNCTIONS_SHA256 = "1d656039b8dca20bac8f056ad195fdc201775331e73e0f6f5da1e866622f1f78"
DRIVER_SHA256 = "03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17"
V2_FUNCTIONS_SHA256 = "00c9c9c51df14ed14b312197594bc855d7f5e2612d62775f35cf6f3cb71ae0cf"


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def replace_once(source: str, old: str, new: str) -> str:
    count = source.count(old)
    if count != 1:
        raise ValueError(f"modern-port anchor count is {count}, expected one: {old[:80]!r}")
    return source.replace(old, new)


def patch_functions_source(source: str) -> str:
    source = replace_once(source, "        gb_x_all.append(np.array(gb_ij_x))\n        gb_y_all.append(np.array(gb_ij_y))", "        gb_x_all.append(gb_ij_x.copy())\n        gb_y_all.append(gb_ij_y.copy())")
    source = replace_once(source, "@nb.jit(parallel=True)\ndef sort_gb", "@nb.njit(fastmath=False)\ndef sort_gb")
    source = source.replace("abs(disxxx)", "np.abs(disxxx)").replace("abs(disyyy)", "np.abs(disyyy)")
    source = replace_once(source, "        DIS = list(dis)\n        min_index = DIS.index(min(DIS))\n\n        if min(DIS) < 1e-4:", "        # Numba compatibility: np.argmin preserves Python list.index(min(...)) first-tie behavior.\n        min_index = np.argmin(dis)\n        minimum_distance = dis[min_index]\n\n        if minimum_distance < 1e-4:")
    source = replace_once(source, "        elif min(DIS) >= 1e-4:", "        elif minimum_distance >= 1e-4:")
    source = replace_once(source, "        loc = np.zeros(20, dtype=int)", "        loc = np.zeros(20, dtype=np.int64)")
    source = replace_once(source, "#@nb.njit(parallel=True)\ndef clockorcounterclock", "@nb.njit(fastmath=False)\ndef clockorcounterclock")
    source = replace_once(source, "        tan_vec = np.array([gb_x_dd, gb_y_dd, 0.])\n\n", "        # Numba compatibility: retain the exact z-component cross-product algebra as scalars.\n")
    source = replace_once(source, "        norm_vec = -np.array([norm_xx, norm_yy, 0.])\n        \n        \n        vvv += np.sign(np.dot(np.cross(tan_vec, norm_vec), np.array([0., 0., 1.])))", "        cross_z = -gb_x_dd * norm_yy + gb_y_dd * norm_xx\n        vvv += np.sign(cross_z)")
    source = replace_once(source, "    if abs(x) < NNN / 2:\n        return x", "    # The three finite scalar branches are exhaustive for native integer offsets.\n    return x")
    return source


def generate(v1_dir: Path, output: Path, comparison: Path | None = None) -> dict:
    driver = (v1_dir / "Bicrystal-4reference-el-pf.py").read_bytes()
    functions = (v1_dir / "functions_4ref_new.py").read_bytes()
    if sha256_bytes(driver) != DRIVER_SHA256 or sha256_bytes(functions) != V1_FUNCTIONS_SHA256:
        raise ValueError("V1 source identity mismatch")
    patched = patch_functions_source(functions.decode()).encode()
    if sha256_bytes(patched) != V2_FUNCTIONS_SHA256:
        raise ValueError("generated V2 source identity mismatch")
    output.mkdir(parents=True, exist_ok=False)
    (output / "Bicrystal-4reference-el-pf.py").write_bytes(driver)
    (output / "functions_4ref_new.py").write_bytes(patched)
    diff = "".join(difflib.unified_diff(
        functions.decode().splitlines(keepends=True), patched.decode().splitlines(keepends=True),
        fromfile="numba_compatibility_v1/functions_4ref_new.py",
        tofile="numba_compatibility_v2/functions_4ref_new.py",
    ))
    validation = json.loads(comparison.read_text()) if comparison else None
    record = {
        "schema": "qiu-modern-numba-compatibility-patch-v2",
        "source_variant": "numba_compatibility_v2",
        "parent_variant": "numba_compatibility_v1",
        "reason": "make the archived boundary ordering path compile in Numba 0.61 nopython mode without changing finite-input semantics",
        "execution_policy": {"authoritative_numba_threads": 1, "reason": "stress_field_extend performs in-place neighbor updates and is nondeterministic under prange scheduling"},
        "parent": {"driver_sha256": DRIVER_SHA256, "functions_sha256": V1_FUNCTIONS_SHA256},
        "patched": {"driver_sha256": sha256_bytes(driver), "functions_sha256": sha256_bytes(patched)},
        "semantic_changes": ["array copies replace unsupported array constructors", "first-tie argmin replaces list.index(min)", "scalar cross-product algebra replaces temporary vectors", "PBC_scalar finite branches have a non-Optional return type"],
        "exact_unified_diff": diff,
        "actual_fixture_validation": validation,
    }
    (output / "native_compatibility_patch_v2.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("v1_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--comparison", type=Path)
    args = parser.parse_args()
    print(json.dumps(generate(args.v1_dir, args.output, args.comparison), indent=2))


if __name__ == "__main__":
    main()
