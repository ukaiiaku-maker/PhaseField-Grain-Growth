#!/usr/bin/env python3
"""Build the bounded, checkpointed Qiu-SI native preflight from qualified sources."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


PRISTINE_DRIVER_SHA256 = "03f8cee8834e5f9669dceb1c831590de55c5ef26f6cc41e0913c0fc89ee2cd17"
FUNCTION_VARIANTS = {
    "pristine_archive": "3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf",
    "numba_compatibility_v1": "1d656039b8dca20bac8f056ad195fdc201775331e73e0f6f5da1e866622f1f78",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"instrumentation anchor count is {source.count(old)}, expected one: {old[:80]!r}")
    return source.replace(old, new)


def instrument(source: str) -> str:
    source = replace_once(source, "import os\nimport numba as nb", "import os\nimport json\nimport numba as nb")
    source = replace_once(
        source,
        "from functions_4ref_new import *\nimport sparse",
        "from functions_4ref_new import *\nimport sparse\n"
        "from qiu_native_preflight_v2_instrumentation import (CHECKPOINT_STEPS, boundary_pairs, morphology, write_checkpoint)",
    )
    source = replace_once(source, "nsteps = 200000              # Max No. Steps", "nsteps = int(os.environ.get('QIU_PREFLIGHT_TARGET_STEP', '1000'))")
    source = replace_once(source, "output_flag = 1                # Enable output", "output_flag = 0                # numeric preflight output only")
    source = replace_once(source, "plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,0)", "if output_flag: plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,0)")
    setup_anchor = 'print("t:",time.time()-t0)\n#########'
    setup = '''print("t:",time.time()-t0)

preflight_output = os.environ.get("QIU_PREFLIGHT_OUTPUT", "preflight-output")
restart_path = os.environ.get("QIU_PREFLIGHT_RESTART", "")
start_step = 1
_, initial_population, _, _, _ = morphology(phi)
if restart_path:
    restart = np.load(restart_path, allow_pickle=False)
    phi[:] = restart["phi"]
    phi_new[:] = restart["phi_new"]
    nf[:] = restart["nf"]
    mf[:] = 0
    mf[:restart["mf_prefix"].shape[0]] = restart["mf_prefix"]
    sigma11_R1[:] = restart["sigma11_R1"]; sigma12_R1[:] = restart["sigma12_R1"]; sigma22_R1[:] = restart["sigma22_R1"]
    sigma11_R2[:] = restart["sigma11_R2"]; sigma12_R2[:] = restart["sigma12_R2"]; sigma22_R2[:] = restart["sigma22_R2"]
    initial_population = restart["initial_population"]
    idov2 = [tuple(x) for x in restart["idov2"]]
    idov3 = [tuple(x) for x in restart["idov3"]]
    start_step = int(restart["accepted_step"]) + 1

def emit_preflight_checkpoint(accepted_step, pre_delta, accepted_delta):
    diagnostic_pairs2, diagnostic_pairs3 = overlaps(phi, number_of_grain)
    nf_diagnostic = np.zeros_like(nf)
    mf_diagnostic = np.zeros((number_of_grain, nx, ny), dtype=mf.dtype)
    update_nfmf(phi, mf_diagnostic, nf_diagnostic, nx, ny, number_of_grain)
    diagnostic_sigma11=np.zeros((nx,ny)); diagnostic_sigma12=np.zeros((nx,ny)); diagnostic_sigma22=np.zeros((nx,ny))
    diagnostic_sigma11_b=np.zeros((nx,ny)); diagnostic_sigma12_b=np.zeros((nx,ny)); diagnostic_sigma22_b=np.zeros((nx,ny))
    diagnostic_gb_x, diagnostic_gb_y = find_gb(phi, nf_diagnostic, mf_diagnostic, nx, ny, dx, dy, number_of_grain, G, diagnostic_sigma11, diagnostic_sigma12, diagnostic_sigma22, diagnostic_sigma11_b, diagnostic_sigma12_b, diagnostic_sigma22_b, sigma12_ext, sigma11_ext, sigma22_ext, ddd, accepted_step, stress_step, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2, misorientation_list, ref_theta_i, diagnostic_pairs2, diagnostic_pairs3)
    pairs = boundary_pairs(phi, diagnostic_pairs2, ddd)
    trajectory_state = None
    if accepted_step == 100:
        trajectory_state = dict(phi=phi, phi_new=phi_new, nf=nf, mf_prefix=mf[:number_of_grain], sigma11_R1=sigma11_R1, sigma12_R1=sigma12_R1, sigma22_R1=sigma22_R1, sigma11_R2=sigma11_R2, sigma12_R2=sigma12_R2, sigma22_R2=sigma22_R2, idov2=np.asarray(idov2,dtype=np.int16), idov3=np.asarray(idov3,dtype=np.int16), accepted_step=np.asarray(accepted_step), initial_population=initial_population)
    write_checkpoint(preflight_output, accepted_step, dt, phi, misorientation_list, eij, pre_delta, accepted_delta, diagnostic_sigma11, diagnostic_sigma12, diagnostic_sigma22, diagnostic_gb_x, diagnostic_gb_y, pairs, ref_theta_i, nf_diagnostic, mf_diagnostic, dx, dy, eta, beta, initial_population=initial_population, native_trial_buffer=phi_new, trajectory_state=trajectory_state)

if start_step == 1:
    zero_delta=np.zeros_like(phi)
    emit_preflight_checkpoint(0, zero_delta, zero_delta)
#########'''
    source = replace_once(source, setup_anchor, setup)
    source = replace_once(source, "for nstep in range(1,nsteps+1):", "for nstep in range(start_step,nsteps+1):")
    update_anchor = "    #update_PF(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt, ETA)\n    update_PF("
    update_replacement = "    capture_checkpoint = nstep in CHECKPOINT_STEPS\n    if capture_checkpoint: phi_before_update = phi.copy()\n    #update_PF(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt, ETA)\n    update_PF("
    source = replace_once(source, update_anchor, update_replacement)
    source = replace_once(
        source,
        "    renorm(phi,phi_new,number_of_grain)\n    #print(\"     Renorm: \",time.time()-t0)",
        "    if capture_checkpoint: pre_renormalization_delta = phi_new - phi_before_update\n"
        "    renorm(phi,phi_new,number_of_grain)\n"
        "    if capture_checkpoint:\n"
        "        accepted_delta = phi - phi_before_update\n"
        "        emit_preflight_checkpoint(nstep, pre_renormalization_delta, accepted_delta)\n"
        "    #print(\"     Renorm: \",time.time()-t0)",
    )
    source += '''
from pathlib import Path as _PreflightPath
_manifest_path=_PreflightPath(preflight_output)/"checkpoint_manifest.json"
_manifest=json.loads(_manifest_path.read_text())
_manifest["baseline_eligible"]=False
_manifest["instrumentation"]="qiu-native-preflight-v2"
_manifest["restart_used"]=bool(restart_path)
_manifest_path.write_text(json.dumps(_manifest,indent=2)+"\\n")
_dispatchers=(Initialize,find_gb,sort_gb,stress_field_line,stress_field_extend,stress_field_bulk,stress_field_bulk_single,update_nfmf,beta,update_PF,E_elastic,grad)
_compilation={f.__name__:{"nopython_signatures":[str(x) for x in f.nopython_signatures],"signatures":[str(x) for x in f.signatures]} for f in _dispatchers}
_object_mode=[name for name,value in _compilation.items() if not value["nopython_signatures"]]
(_PreflightPath(preflight_output)/"compilation_manifest.json").write_text(json.dumps({"dispatchers":_compilation,"object_mode_or_uncompiled":_object_mode},indent=2)+"\\n")
if _object_mode: raise RuntimeError("native dispatcher lacks a nopython signature: "+", ".join(_object_mode))
print("QIU_NATIVE_PREFLIGHT_SEGMENT_COMPLETE", nsteps)
'''
    return source


def build(source_dir: Path, instrumentation: Path, output: Path, source_variant: str) -> dict:
    driver = source_dir / "Bicrystal-4reference-el-pf.py"
    functions = source_dir / "functions_4ref_new.py"
    if sha256(driver) != PRISTINE_DRIVER_SHA256:
        raise ValueError("driver hash mismatch")
    if source_variant not in FUNCTION_VARIANTS:
        raise ValueError("unknown source variant")
    functions_sha256 = FUNCTION_VARIANTS[source_variant]
    if sha256(functions) != functions_sha256:
        raise ValueError("functions hash mismatch")
    output.mkdir(parents=True, exist_ok=False)
    generated = output / "Bicrystal-4reference-el-pf-preflight-v2.py"
    generated.write_text(instrument(driver.read_text()))
    shutil.copy2(functions, output / functions.name)
    shutil.copy2(instrumentation, output / instrumentation.name)
    manifest = {
        "schema": "qiu-native-preflight-v2-build",
        "baseline_eligible": False,
        "source_variant": source_variant,
        "pristine_driver_sha256": PRISTINE_DRIVER_SHA256,
        "functions_sha256": functions_sha256,
        "instrumented_driver_sha256": sha256(generated),
        "instrumentation_sha256": sha256(output / instrumentation.name),
        "checkpoints": [0, 1, 10, 100, 1000],
        "restart_checkpoint": 100,
    }
    (output / "preflight_build_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("instrumentation", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-variant", choices=sorted(FUNCTION_VARIANTS), default="numba_compatibility_v1")
    args = parser.parse_args()
    print(json.dumps(build(args.source_dir, args.instrumentation, args.output, args.source_variant), indent=2))


if __name__ == "__main__":
    main()
