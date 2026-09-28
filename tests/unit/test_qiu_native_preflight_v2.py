import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).parents[2] / "native_compatibility" / "numba_compatibility_v1"))

from scripts.qiu_native_preflight_v2_instrumentation import CHECKPOINT_STEPS
from scripts.build_qiu_native_preflight_v2 import instrument


def test_exact_native_preflight_checkpoint_sequence():
    assert CHECKPOINT_STEPS == (0, 1, 10, 100, 1000)


def test_instrumentation_requires_unique_anchors():
    try:
        instrument("not a native driver")
    except ValueError as error:
        assert "anchor count" in str(error)
    else:
        raise AssertionError("invalid source was accepted")
