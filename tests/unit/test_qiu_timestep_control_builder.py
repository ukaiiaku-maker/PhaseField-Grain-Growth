from scripts.build_qiu_timestep_control_v1 import build


def test_half_dt_builder_changes_only_registered_timestep(tmp_path):
    source=tmp_path/"source";source.mkdir()
    # Exercise the actual immutable source through the public builder in a temporary output.
    from pathlib import Path
    repo=Path(__file__).resolve().parents[2]
    result=build(repo/"native_compatibility/qiu_anisotropy_v1",repo/"scripts/qiu_native_preflight_v3_instrumentation.py",tmp_path/"build")
    text=(tmp_path/"build/Bicrystal-4reference-el-pf-control-v1.py").read_text()
    assert "QIU_CONTROL_DT" in text and "dt = 0.1                            # Time Step" not in text
    assert result["accepted_dt"]==0.05 and result["target_steps"]==10000
