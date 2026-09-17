import importlib.util
import json
from pathlib import Path
import sys
import types

import numpy as np

from scripts.qiu_modern_port_v2 import V2_FUNCTIONS_SHA256, generate, sha256_bytes


ROOT = Path(__file__).parents[2]
V1 = ROOT / "native_compatibility" / "numba_compatibility_v1"
V2 = ROOT / "native_compatibility" / "numba_compatibility_v2"
COMPARISON = V2 / "modern_port_v2_actual_comparison.json"


def load_v2():
    sys.modules.setdefault("sparse", types.ModuleType("sparse"))
    spec = importlib.util.spec_from_file_location("qiu_modern_v2_functions", V2 / "functions_4ref_new.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v2_is_reproducibly_generated(tmp_path):
    output = tmp_path / "v2"
    record = generate(V1, output, COMPARISON)
    assert sha256_bytes((output / "functions_4ref_new.py").read_bytes()) == V2_FUNCTIONS_SHA256
    assert (output / "functions_4ref_new.py").read_bytes() == (V2 / "functions_4ref_new.py").read_bytes()
    assert record["actual_fixture_validation"]["passed"] is True
    assert record["execution_policy"]["authoritative_numba_threads"] == 1


def test_pbc_and_cross_product_helpers_compile_nopython_with_reference_semantics():
    functions = load_v2()
    values = np.array([-11, -5, -4, 0, 4, 5, 11], dtype=np.int64)
    expected = np.array([-1, 5, -4, 0, 4, -5, 1], dtype=np.int64)
    actual = np.array([functions.PBC_scalar(x, 10) for x in values])
    assert np.array_equal(actual, expected)
    assert functions.PBC_scalar.nopython_signatures

    phi = np.zeros((1, 8, 8), dtype=np.float64)
    phi[0, 2, 2] = 2.0
    x = np.array([1, 2, 3], dtype=np.int64)
    y = np.array([2, 2, 2], dtype=np.int64)
    loc = np.array([1], dtype=np.int64)
    expected_turn = functions.clockorcounterclock.py_func(x, y, phi, 8, 8, loc, 0, 1, 1, 0)
    assert functions.clockorcounterclock(x, y, phi, 8, 8, loc, 0, 1, 1, 0) == expected_turn
    assert functions.clockorcounterclock.nopython_signatures


def test_actual_initialization_fixture_comparison_is_exact_and_single_threaded():
    comparison = json.loads(COMPARISON.read_text())
    assert comparison["passed"] is True
    assert comparison["execution_threads"] == 1
    assert comparison["repeat_exact"] is True
    assert comparison["fields"]
    assert all(field["exact"] for field in comparison["fields"])
    assert comparison["patched_find_gb_nopython_signatures"]
