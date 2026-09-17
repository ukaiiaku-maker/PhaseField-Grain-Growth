import hashlib
from pathlib import Path

from scripts.build_qiu_control_v1 import build,instrument


ROOT=Path(__file__).parents[2];SOURCE=ROOT/"native_compatibility/qiu_anisotropy_v1";INSTRUMENTATION=ROOT/"scripts/qiu_native_preflight_v3_instrumentation.py"


def test_control_driver_selects_native_a0_or_frozen_a2_without_changing_other_terms():
    generated=instrument((SOURCE/"Bicrystal-4reference-el-pf.py").read_text())
    assert "update_PF_control(phi,phi_new" in generated
    assert "QIU_CONTROL_CODE" in generated and "QIU_ENERGY_NORMALIZATION" in generated and "QIU_MOBILITY_NORMALIZATION" in generated
    assert "0.65,0.85,16,2.0,1.0910512514090829" in generated
    assert '_manifest["baseline_eligible"]=False' in generated
    assert "QIU_CONTROL_SEGMENT_COMPLETE" in generated


def test_control_build_records_all_four_cases(tmp_path):
    record=build(SOURCE,INSTRUMENTATION,tmp_path/"build")
    assert record["controls"]=={"QIU_SI_4REF_A0_PORT":0,"QIU_SI_4REF_ANISO_E":1,"QIU_SI_4REF_ANISO_M":2,"QIU_SI_4REF_ANISO_EM_INV":3}
    assert record["windows"]==[100,1000,5000]
    assert len(record["instrumented_driver_sha256"])==64
