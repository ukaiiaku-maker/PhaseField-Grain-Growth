from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "native_compatibility" / "numba_compatibility_v1"))

from scripts.qiu_native_preflight_v3_instrumentation import CHECKPOINT_STEPS
from scripts.build_qiu_native_preflight_v3 import instrument


def test_historical_preflight_checkpoint_sequence_and_restart():
    assert CHECKPOINT_STEPS == (0, 1, 10, 100, 500, 1000)
    source = Path("native_compatibility/numba_compatibility_v1/Bicrystal-4reference-el-pf.py").read_text()
    generated = instrument(source)
    assert "if accepted_step in (100, 500):" in generated
    assert "diagnostic_sigma11_R1=sigma11_R1.copy()" in generated
    assert "QIU_NATIVE_HISTORICAL_OBJECT_MODE" in generated
    assert "native dispatcher lacks a nopython signature" not in generated


def test_v3_instrumentation_requires_unique_anchors():
    try:
        instrument("not a native driver")
    except ValueError as error:
        assert "anchor count" in str(error)
    else:
        raise AssertionError("invalid source was accepted")
