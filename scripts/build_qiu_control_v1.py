#!/usr/bin/env python3
"""Build checkpointed A0 and A2 control drivers on the native Qiu network."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

try:
    from scripts.build_qiu_native_production_v1 import instrument as production_instrument
except ModuleNotFoundError:
    from build_qiu_native_production_v1 import instrument as production_instrument


DRIVER_SHA256="03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17"
FUNCTIONS_SHA256="9ce5fb1eec86498cc84cc6bed7c6e25dd37abcb87cecad096064b5a9e1acb78a"


def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(source: str,old: str,new: str) -> str:
    if source.count(old)!=1: raise ValueError(f"control instrumentation anchor count is {source.count(old)}, expected one: {old[:100]!r}")
    return source.replace(old,new)


def instrument(source: str) -> str:
    generated=production_instrument(source)
    generated=replace_once(generated,"nsteps = int(os.environ.get('QIU_SEGMENT_TARGET_STEP', '200000'))","nsteps = int(os.environ.get('QIU_CONTROL_TARGET_STEP', '5000'))")
    generated=replace_once(generated,'preflight_output = os.environ.get("QIU_PRODUCTION_OUTPUT", "production-output")\nrestart_path = os.environ.get("QIU_PRODUCTION_RESTART", "")','preflight_output = os.environ.get("QIU_CONTROL_OUTPUT", "control-output")\nrestart_path = os.environ.get("QIU_CONTROL_RESTART", "")\ncontrol_name = os.environ.get("QIU_CONTROL_NAME", "QIU_SI_4REF_A0_PORT")\ncontrol_code = int(os.environ.get("QIU_CONTROL_CODE", "0"))\nenergy_normalization = float(os.environ.get("QIU_ENERGY_NORMALIZATION", "1"))\nmobility_normalization = float(os.environ.get("QIU_MOBILITY_NORMALIZATION", "1"))\nif control_code not in (0,1,2,3): raise ValueError("invalid QIU_CONTROL_CODE")')
    old="    update_PF(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt)"
    new="    update_PF_control(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt,control_code,0.65,0.85,16,2.0,1.0910512514090829,energy_normalization,mobility_normalization)"
    generated=replace_once(generated,old,new)
    generated=generated.replace('_manifest["baseline_eligible"]=True','_manifest["baseline_eligible"]=False')
    generated=generated.replace('_manifest["instrumentation"]="qiu-native-production-v1"','_manifest["instrumentation"]="qiu-control-v1"')
    generated=generated.replace('_manifest["logical_trajectory"]="QIU_SI_NATIVE_BASELINE_V2"','_manifest["logical_trajectory"]=control_name\n_manifest["control_code"]=control_code\n_manifest["energy_normalization"]=energy_normalization\n_manifest["mobility_normalization"]=mobility_normalization')
    generated=generated.replace('_dispatchers=(Initialize,find_gb,sort_gb,stress_field_line,stress_field_extend,stress_field_bulk,stress_field_bulk_single,update_nfmf,beta,update_PF,E_elastic,grad)','_dispatchers=(Initialize,find_gb,sort_gb,stress_field_line,stress_field_extend,stress_field_bulk,stress_field_bulk_single,update_nfmf,beta,update_PF_control,qiu_anisotropic_pair_drive_cell,qiu_anisotropic_mobility,E_elastic,grad)')
    generated=generated.replace('print("QIU_NATIVE_PRODUCTION_SEGMENT_COMPLETE", start_step, nsteps)','print("QIU_CONTROL_SEGMENT_COMPLETE", control_name, start_step, nsteps)')
    return generated


def build(source_dir: Path,instrumentation: Path,output: Path) -> dict:
    driver=source_dir/"Bicrystal-4reference-el-pf.py";functions=source_dir/"functions_4ref_new.py"
    if sha256(driver)!=DRIVER_SHA256 or sha256(functions)!=FUNCTIONS_SHA256: raise ValueError("anisotropic control source identity mismatch")
    output.mkdir(parents=True,exist_ok=False);generated=output/"Bicrystal-4reference-el-pf-control-v1.py";generated.write_text(instrument(driver.read_text()));shutil.copy2(functions,output/functions.name);shutil.copy2(instrumentation,output/instrumentation.name)
    manifest={"schema":"qiu-control-v1-build","baseline_eligible":False,"source_variant":"qiu_anisotropy_v1","driver_sha256":DRIVER_SHA256,"functions_sha256":FUNCTIONS_SHA256,"instrumented_driver_sha256":sha256(generated),"instrumentation_sha256":sha256(output/instrumentation.name),"controls":{"QIU_SI_4REF_A0_PORT":0,"QIU_SI_4REF_ANISO_E":1,"QIU_SI_4REF_ANISO_M":2,"QIU_SI_4REF_ANISO_EM_INV":3},"windows":[100,1000,5000],"execution_threads":1}
    (output/"control_build_manifest.json").write_text(json.dumps(manifest,indent=2)+'\n');return manifest


def main():
    parser=argparse.ArgumentParser();parser.add_argument("source_dir",type=Path);parser.add_argument("instrumentation",type=Path);parser.add_argument("output",type=Path);args=parser.parse_args();print(json.dumps(build(args.source_dir,args.instrumentation,args.output),indent=2))


if __name__=="__main__": main()
