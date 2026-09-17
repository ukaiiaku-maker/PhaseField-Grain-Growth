#!/usr/bin/env python3
"""Prepare an immutable A0/A2 Qiu control window."""
from __future__ import annotations
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np

CONTROLS={"QIU_SI_4REF_A0_PORT":0,"QIU_SI_4REF_ANISO_E":1,"QIU_SI_4REF_ANISO_M":2,"QIU_SI_4REF_ANISO_EM_INV":3}
IDENTITIES={"PF_Codes.zip":"2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90","sparse-0.19.2-py2.py3-none-any.whl":"2fcad7b83dacb5c70123a7de76db3626cbdaf567fcf12bdbf637aca6c4ee1bbe","build/Bicrystal-4reference-el-pf-control-v1.py":"a1d3b971041c4e13aa288c72867690ce5d0e29265b11daab0e52418453d75db1","build/functions_4ref_new.py":"9ce5fb1eec86498cc84cc6bed7c6e25dd37abcb87cecad096064b5a9e1acb78a","build/control_build_manifest.json":"96a816a9423236a25a951cd5a3ef3f5dc097cb34010540846d52946bf20d8dc5","build/qiu_native_preflight_v3_instrumentation.py":"a94dddb70806b86a76dac16a12b76c8b0d00d537940f9d33b4717f857d98b943"}
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def step(path):
 with np.load(path,allow_pickle=False) as z:return int(z["accepted_step"])
def prepare(base,output,control,start,target,energy_norm,mobility_norm,parent=None):
 if control not in CONTROLS or (start,target) not in {(0,100),(100,1000),(1000,5000)}:raise ValueError("unregistered control window")
 if (start==0)!=(parent is None):raise ValueError("parent required exactly for continuation windows")
 for rel,digest in IDENTITIES.items():
  if sha(base/rel)!=digest:raise ValueError(f"base identity mismatch: {rel}")
 parent_sha="";parent_name=""
 if parent:
  if step(parent)!=start:raise ValueError("parent step mismatch")
  parent_sha=sha(parent);parent_name=f"parent/{parent.name}"
 output.mkdir(parents=True,exist_ok=False);shutil.copy2(base/"PF_Codes.zip",output/"PF_Codes.zip");shutil.copy2(base/"sparse-0.19.2-py2.py3-none-any.whl",output/"sparse-0.19.2-py2.py3-none-any.whl");shutil.copytree(base/"build",output/"build")
 if parent:(output/"parent").mkdir();shutil.copy2(parent,output/"parent"/parent.name)
 restart=f'''test "$(sha256sum {parent_name} | awk '{{print $1}}')" = {parent_sha}\nexport QIU_CONTROL_RESTART="$PWD/{parent_name}"''' if parent else "unset QIU_CONTROL_RESTART"
 wrapper=f'''#!/bin/bash
set -Eeuo pipefail
environment=/pub/sdillon1/sw/pfgg-anisotropy/f4aa53f-714dd16b2f3624b8
python_bin="$environment/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 QIU_NUMBA_THREADS=1 PYTHONHASHSEED=0 MPLBACKEND=Agg MPLCONFIGDIR="$PWD/mpl" PYTHONPATH="$PWD/third_party:$PWD/build"
mkdir -p output/checkpoints "$MPLCONFIGDIR" third_party
{restart}
"$python_bin" -m pip install --no-deps --no-index --target third_party sparse-0.19.2-py2.py3-none-any.whl >output/pip.log 2>&1
cp build/control_build_manifest.json output/
export QIU_CONTROL_OUTPUT="$PWD/output/checkpoints" QIU_CONTROL_NAME={control} QIU_CONTROL_CODE={CONTROLS[control]} QIU_CONTROL_TARGET_STEP={target} QIU_CHECKPOINT_CADENCE=250 QIU_ENERGY_NORMALIZATION={energy_norm:.17g} QIU_MOBILITY_NORMALIZATION={mobility_norm:.17g}
/usr/bin/time -v -o output/time.txt "$python_bin" -u build/Bicrystal-4reference-el-pf-control-v1.py >output/stdout.log 2>output/stderr.log
"$python_bin" - <<'PY'
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('output');c=root/'checkpoints';start={start};target={target};control='{control}'
m=json.loads((c/'checkpoint_manifest.json').read_text());assert m['logical_trajectory']==control and m['baseline_eligible'] is False and m['restart_used'] is {bool(parent)};assert all(x['finite'] for x in m['checkpoints'])
comp=json.loads((c/'compilation_manifest.json').read_text());assert comp['object_mode_or_uncompiled']==[]
r=c/f'native-restart-step{{target:06d}}.npz';assert r.is_file()
with np.load(r,allow_pickle=False) as z:assert int(z['accepted_step'])==target
d={{'schema':'qiu-control-segment-v1','classification':'QIU_SI_CONTROL_WINDOW_PASSED','control':control,'control_code':{CONTROLS[control]},'segment_start':start,'segment_target':target,'parent_checkpoint_sha256':'{parent_sha}','end_checkpoint_sha256':hashlib.sha256(r.read_bytes()).hexdigest(),'energy_normalization':{energy_norm!r},'mobility_normalization':{mobility_norm!r},'execution_threads':1,'finite':True}}
(root/'segment_decision.json').write_text(json.dumps(d,indent=2)+'\n')
PY
find output -type f -print0|sort -z|xargs -0 sha256sum >output/payload-files.sha256
'''
 (output/"run_control_segment.sh").write_text(wrapper);(output/"run_control_segment.sh").chmod(0o755)
 inputs='"run_control_segment.sh","PF_Codes.zip","sparse-0.19.2-py2.py3-none-any.whl","build"'+(',"parent"' if parent else '')
 slug=control.lower().replace('qiu_si_4ref_','').replace('_','-')
 (output/"job.toml").write_text(f'''name = "pfgg-qiu-{slug}-{start:04d}-{target:04d}"
entrypoint = ["bash","run_control_segment.sh"]
working_directory = "."
modules = []
environment = {{ PYTHONHASHSEED = "0" }}
input_paths = [{inputs}]
output_paths = ["output"]
result_mode = "archive"
scratch_mode = "local"
input_mode = "archive"
remote_inputs = {{}}
checkpoint_command = []
contains_sensitive_data = false
[resources]
account = "SDILLON1_LAB"
partition = "standard"
time = "06:00:00"
nodes = 1
ntasks = 1
cpus_per_task = 1
memory = "32G"
temporary_space = "100G"
constraint = ""
gres = ""
qos = ""
[notifications]
enabled = false
email = "sdillon1@uci.edu"
events = ["FAIL","END"]
''')
 record={'schema':'qiu-control-segment-preparation-v1','control':control,'control_code':CONTROLS[control],'start':start,'target':target,'parent_checkpoint_sha256':parent_sha,'energy_normalization':energy_norm,'mobility_normalization':mobility_norm,'wrapper_sha256':sha(output/'run_control_segment.sh')};(output/'segment_preparation.json').write_text(json.dumps(record,indent=2)+'\n');return record
def main():
 p=argparse.ArgumentParser();p.add_argument('base',type=Path);p.add_argument('output',type=Path);p.add_argument('--control',choices=CONTROLS,required=True);p.add_argument('--start',type=int,required=True);p.add_argument('--target',type=int,required=True);p.add_argument('--energy-normalization',type=float,required=True);p.add_argument('--mobility-normalization',type=float,required=True);p.add_argument('--parent',type=Path);a=p.parse_args();print(json.dumps(prepare(a.base,a.output,a.control,a.start,a.target,a.energy_normalization,a.mobility_normalization,a.parent),indent=2))
if __name__=='__main__':main()
