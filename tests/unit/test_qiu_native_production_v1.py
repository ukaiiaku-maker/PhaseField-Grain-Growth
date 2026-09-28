from pathlib import Path

import hashlib

from scripts.build_qiu_native_production_v1 import build, instrument


ROOT = Path(__file__).parents[2]
SOURCE = ROOT / "native_compatibility" / "numba_compatibility_v1"
INSTRUMENTATION = ROOT / "scripts" / "qiu_native_preflight_v3_instrumentation.py"


def test_production_instrumentation_uses_absolute_segment_targets_and_exact_end_restart():
    generated = instrument((SOURCE / "Bicrystal-4reference-el-pf.py").read_text())
    assert "QIU_SEGMENT_TARGET_STEP" in generated
    assert "QIU_PRODUCTION_RESTART" in generated
    assert "capture_checkpoint = nstep == nsteps or nstep % checkpoint_cadence == 0" in generated
    assert "if accepted_step == nsteps:" in generated
    assert 'manifest["baseline_eligible"]=True' in generated
    assert 'QIU_NATIVE_PRODUCTION_SEGMENT_COMPLETE' in generated
    assert "nb.set_num_threads(int(os.environ.get('QIU_NUMBA_THREADS', '1')))" in generated


def test_production_instrumented_driver_identity_is_reproducible():
    generated = instrument((SOURCE / "Bicrystal-4reference-el-pf.py").read_text()).encode()
    assert hashlib.sha256(generated).hexdigest() == "aeeee2310b38db054df76dd3fc548a2f0dbf469734e2a2dcd4678e86f1b6c5a4"


def test_production_build_rejects_nonpristine_functions(tmp_path):
    try:
        build(SOURCE, INSTRUMENTATION, tmp_path / "build")
    except ValueError as error:
        assert "source identity mismatch" in str(error)
    else:
        raise AssertionError("patched functions were accepted as pristine")


def test_production_instrumentation_rejects_wrong_source():
    try:
        instrument("not the native source")
    except ValueError as error:
        assert "anchor count" in str(error)
    else:
        raise AssertionError("invalid source was accepted")
