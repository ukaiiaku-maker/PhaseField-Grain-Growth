#!/usr/bin/env python3
"""Prepare an immutable half-dt combined-anisotropy qualification segment."""
from __future__ import annotations
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np

CONTROL="QIU_SI_4REF_ANISO_EM_INV_DT_HALF"
IDENTITIES={"PF_Codes.zip":"2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90","sparse-0.19.2-py2.py3-none-any.whl":"2fcad7b83dacb5c70123a7de76db3626cbdaf567fcf12bdbf637aca6c4ee1bbe","build/Bicrystal-4reference-el-pf-control-v1.py":"34ffa08f24aeeaac4f28409170bcbfd89ca39a4564809ffea13ab2b06ee7b27d","build/functions_4ref_new.py":"9ce5fb1eec86498cc84cc6bed7c6e25dd37abcb87cecad096064b5a9e1acb78a","build/control_build_manifest.json":"d63fb5aacbaca387dbf4dece94f6bc24db767fd1752bedc4262fc4bc55e31e17"}
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def step(path):
 with np.load(path,allow_pickle=False) as z:return int(z["accepted_step"])
def prepare(base,output,start,target,energy,mobility,parent=None):
 if (start,target) not in {(0,5000),(5000,10000)}:raise ValueError("unregistered timestep window")
 if (start==0)!=(parent is None):raise ValueError("parent required exactly for continuation")
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
export QIU_CONTROL_OUTPUT="$PWD/output/checkpoints" QIU_CONTROL_NAME={CONTROL} QIU_CONTROL_CODE=3 QIU_CONTROL_TARGET_STEP={target} QIU_CONTROL_DT=0.05 QIU_CHECKPOINT_CADENCE=250 QIU_ENERGY_NORMALIZATION={energy:.17g} QIU_MOBILITY_NORMALIZATION={mobility:.17g}
/usr/bin/time -v -o output/time.txt "$python_bin" -u build/Bicrystal-4reference-el-pf-control-v1.py >output/stdout.log 2>output/stderr.log
"$python_bin" - <<'PY'
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('output');c=root/'checkpoints';start={start};target={target}
m=json.loads((c/'checkpoint_manifest.json').read_text());assert m['logical_trajectory']=='{CONTROL}' and m['restart_used'] is {bool(parent)} and all(x['finite'] for x in m['checkpoints'])
comp=json.loads((c/'compilation_manifest.json').read_text());assert comp['object_mode_or_uncompiled']==[]
r=c/f'native-restart-step{{target:06d}}.npz';assert r.is_file()
with np.load(r,allow_pickle=False) as z:assert int(z['accepted_step'])==target
d={{'schema':'qiu-timestep-segment-v1','classification':'QIU_SI_TIMESTEP_SEGMENT_PASSED','control':'{CONTROL}','accepted_dt':0.05,'segment_start':start,'segment_target':target,'parent_checkpoint_sha256':'{parent_sha}','end_checkpoint_sha256':hashlib.sha256(r.read_bytes()).hexdigest(),'energy_normalization':{energy!r},'mobility_normalization':{mobility!r},'execution_threads':1,'finite':True}}
(root/'segment_decision.json').write_text(json.dumps(d,indent=2)+'\n')
PY
find output -type f -print0|sort -z|xargs -0 sha256sum >output/payload-files.sha256
'''
 (output/"run_timestep_segment.sh").write_text(wrapper);(output/"run_timestep_segment.sh").chmod(0o755)
 inputs='"run_timestep_segment.sh","PF_Codes.zip","sparse-0.19.2-py2.py3-none-any.whl","build"'+(',"parent"' if parent else '')
 (output/"job.toml").write_text(f'''name = "pfgg-qiu-aniso-em-dt-half-{start:05d}-{target:05d}"
entrypoint = ["bash","run_timestep_segment.sh"]
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
time = "12:00:00"
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
 record={'schema':'qiu-timestep-segment-preparation-v1','start':start,'target':target,'accepted_dt':0.05,'parent_checkpoint_sha256':parent_sha,'energy_normalization':energy,'mobility_normalization':mobility,'wrapper_sha256':sha(output/'run_timestep_segment.sh')};(output/'segment_preparation.json').write_text(json.dumps(record,indent=2)+'\n');return record
def main():
 p=argparse.ArgumentParser();p.add_argument('base',type=Path);p.add_argument('output',type=Path);p.add_argument('--start',type=int,required=True);p.add_argument('--target',type=int,required=True);p.add_argument('--energy-normalization',type=float,required=True);p.add_argument('--mobility-normalization',type=float,required=True);p.add_argument('--parent',type=Path);a=p.parse_args();print(json.dumps(prepare(a.base,a.output,a.start,a.target,a.energy_normalization,a.mobility_normalization,a.parent),indent=2))
if __name__=='__main__':main()
