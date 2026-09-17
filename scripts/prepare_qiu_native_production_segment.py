#!/usr/bin/env python3
"""Prepare one immutable historical-runtime Qiu native production segment."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np


ARCHIVE_SHA256 = "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90"
DRIVER_SHA256 = "aeeee2310b38db054df76dd3fc548a2f0dbf469734e2a2dcd4678e86f1b6c5a4"
FUNCTIONS_SHA256 = "3fb625fcb88be515defb813df198f671b45f528e44f32ac7fabae36e1c926aaf"
INSTRUMENTATION_SHA256 = "a94dddb70806b86a76dac16a12b76c8b0d00d537940f9d33b4717f857d98b943"
BUILD_MANIFEST_SHA256 = "f2aaa1b8099c474259983e1aa798171fb044d89a1fb3a0a3dce7ac220bfc68c5"
SPARSE_WHEEL_SHA256 = "95ed0b649a0663b1488756ad4cf242b0a9bb2c9a25bc752a7c6ca9fbe8258966"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def accepted_step(path: Path) -> int:
    with np.load(path, allow_pickle=False) as state:
        return int(np.asarray(state["accepted_step"]).reshape(()))


def wrapper_text(start: int, target: int, parent_name: str, parent_sha256: str) -> str:
    restart = ""
    if parent_name:
        restart = f'''test "$(sha256sum {parent_name} | awk '{{print $1}}')" = {parent_sha256}
export QIU_PRODUCTION_RESTART="$PWD/{parent_name}"
'''
    else:
        restart = "unset QIU_PRODUCTION_RESTART\n"
    expected_restart = f"native-restart-step{target:06d}.npz"
    expected_state = f"native-state-step{target:06d}.npz"
    return f'''#!/bin/bash
set -Eeuo pipefail
module purge
module load anaconda/2022.05
python_bin="$(command -v python)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 QIU_NUMBA_THREADS=1 PYTHONHASHSEED=0
export MPLBACKEND=Agg MPLCONFIGDIR="$PWD/mpl" PYTHONPATH="$PWD/third_party:$PWD/build"
mkdir -p output/checkpoints "$MPLCONFIGDIR" third_party
test "$(sha256sum PF_Codes.zip | awk '{{print $1}}')" = {ARCHIVE_SHA256}
test "$(sha256sum build/Bicrystal-4reference-el-pf-production-v1.py | awk '{{print $1}}')" = {DRIVER_SHA256}
test "$(sha256sum build/functions_4ref_new.py | awk '{{print $1}}')" = {FUNCTIONS_SHA256}
test "$(sha256sum build/qiu_native_preflight_v3_instrumentation.py | awk '{{print $1}}')" = {INSTRUMENTATION_SHA256}
test "$(sha256sum build/production_build_manifest.json | awk '{{print $1}}')" = {BUILD_MANIFEST_SHA256}
test "$(sha256sum sparse-0.13.0-py2.py3-none-any.whl | awk '{{print $1}}')" = {SPARSE_WHEEL_SHA256}
{restart}"$python_bin" -m pip install --no-deps --no-index --target third_party sparse-0.13.0-py2.py3-none-any.whl >output/pip.log 2>&1
"$python_bin" - <<'PY' >output/environment_manifest.json
import contextlib,importlib.metadata,io,json,platform,sys,numpy
capture=io.StringIO()
with contextlib.redirect_stdout(capture): numpy.__config__.show()
names=('numpy','scipy','numba','llvmlite','matplotlib','sparse')
print(json.dumps({{'runtime_variant':'anaconda_2022_05_historical_object_mode','python':platform.python_version(),'executable':sys.executable,'packages':{{x:importlib.metadata.version(x) for x in names}},'blas':capture.getvalue(),'execution_threads':1}},indent=2))
PY
cp build/production_build_manifest.json output/
export QIU_PRODUCTION_OUTPUT="$PWD/output/checkpoints" QIU_SEGMENT_TARGET_STEP={target} QIU_CHECKPOINT_CADENCE=250
/usr/bin/time -v -o output/time-segment-{start:06d}-{target:06d}.txt "$python_bin" -u build/Bicrystal-4reference-el-pf-production-v1.py >output/stdout-segment-{start:06d}-{target:06d}.log 2>output/stderr-segment-{start:06d}-{target:06d}.log
test -f output/checkpoints/{expected_restart}
test -f output/checkpoints/{expected_state}
"$python_bin" - <<'PY'
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('output'); cdir=root/'checkpoints'; start={start}; target={target}
manifest=json.loads((cdir/'checkpoint_manifest.json').read_text())
expected=(([0] if start==0 else []) + list(range(((start//250)+1)*250,target+1,250)))
assert [x['step'] for x in manifest['checkpoints']]==expected
assert manifest['baseline_eligible'] is True and manifest['restart_used'] is {str(bool(parent_name))}
assert manifest['logical_trajectory']=='QIU_SI_NATIVE_BASELINE_V2'
assert all(x['finite'] for x in manifest['checkpoints'])
for item in manifest['checkpoints']:
 p=cdir/item['path']; assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
restart=cdir/'{expected_restart}'
with np.load(restart,allow_pickle=False) as state: assert int(state['accepted_step'])==target
decision={{'schema':'qiu-native-production-segment-v1','classification':'QIU_SI_NATIVE_BASELINE_SEGMENT_PASSED','logical_trajectory':'QIU_SI_NATIVE_BASELINE_V2','segment_start':start,'segment_target':target,'parent_checkpoint_sha256':'{parent_sha256}','end_checkpoint':restart.name,'end_checkpoint_sha256':hashlib.sha256(restart.read_bytes()).hexdigest(),'execution_threads':1,'finite':True}}
(root/'segment_decision.json').write_text(json.dumps(decision,indent=2)+'\n')
artifacts=[]
for p in sorted(x for x in root.rglob('*') if x.is_file() and x.name not in {{'artifact_manifest.json','payload-files.sha256'}}):
 artifacts.append({{'path':str(p.relative_to(root)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}})
(root/'artifact_manifest.json').write_text(json.dumps({{'schema':'qiu-native-production-segment-artifacts-v1','artifacts':artifacts}},indent=2)+'\n')
PY
find output -type f ! -name payload-files.sha256 -print0 | sort -z | xargs -0 sha256sum >output/payload-files.sha256
'''


def prepare(base: Path, output: Path, start: int, target: int, walltime: str, parent: Path | None) -> dict:
    if start < 0 or target <= start or target > 200_000 or target % 250 or start % 250:
        raise ValueError("segment bounds must be aligned positive progress within step 200000")
    if (start == 0) != (parent is None):
        raise ValueError("only the first segment may omit its parent checkpoint")
    parent_sha = ""
    parent_name = ""
    if parent is not None:
        if accepted_step(parent) != start:
            raise ValueError("parent checkpoint accepted_step does not match segment start")
        parent_sha = sha256(parent)
        parent_name = f"parent/{parent.name}"
    identities = {
        base / "PF_Codes.zip": ARCHIVE_SHA256,
        base / "sparse-0.13.0-py2.py3-none-any.whl": SPARSE_WHEEL_SHA256,
        base / "build/Bicrystal-4reference-el-pf-production-v1.py": DRIVER_SHA256,
        base / "build/functions_4ref_new.py": FUNCTIONS_SHA256,
        base / "build/qiu_native_preflight_v3_instrumentation.py": INSTRUMENTATION_SHA256,
        base / "build/production_build_manifest.json": BUILD_MANIFEST_SHA256,
    }
    for path, expected in identities.items():
        if sha256(path) != expected:
            raise ValueError(f"base identity mismatch: {path}")
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(base / "PF_Codes.zip", output / "PF_Codes.zip")
    shutil.copy2(base / "sparse-0.13.0-py2.py3-none-any.whl", output / "sparse-0.13.0-py2.py3-none-any.whl")
    shutil.copytree(base / "build", output / "build")
    if parent is not None:
        (output / "parent").mkdir(); shutil.copy2(parent, output / "parent" / parent.name)
    wrapper = output / "run_native_production_segment.sh"
    wrapper.write_text(wrapper_text(start, target, parent_name, parent_sha)); wrapper.chmod(0o755)
    inputs = '"run_native_production_segment.sh", "PF_Codes.zip", "sparse-0.13.0-py2.py3-none-any.whl", "build"'
    if parent is not None: inputs += ', "parent"'
    (output / "job.toml").write_text(f'''name = "pfgg-qiu-native-{start:06d}-{target:06d}"
entrypoint = ["bash", "run_native_production_segment.sh"]
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
time = "{walltime}"
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
events = ["FAIL", "END"]
''')
    record = {"schema":"qiu-native-production-segment-preparation-v1","start":start,"target":target,"walltime":walltime,"parent_checkpoint":parent_name,"parent_checkpoint_sha256":parent_sha,"wrapper_sha256":sha256(wrapper)}
    (output / "segment_preparation.json").write_text(json.dumps(record,indent=2)+"\n")
    return record


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("base",type=Path);parser.add_argument("output",type=Path);parser.add_argument("--start",type=int,required=True);parser.add_argument("--target",type=int,required=True);parser.add_argument("--walltime",default="12:00:00");parser.add_argument("--parent",type=Path)
    args=parser.parse_args();print(json.dumps(prepare(args.base,args.output,args.start,args.target,args.walltime,args.parent),indent=2))


if __name__ == "__main__": main()
