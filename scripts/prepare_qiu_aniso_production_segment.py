#!/usr/bin/env python3
"""Prepare one immutable exact-restart segment of combined Qiu anisotropy."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np


CONTROL = "QIU_SI_4REF_ANISO_EM_INV"
TERMINAL_STEP = 200_000
IDENTITIES = {
    "PF_Codes.zip": "2740e0af26bd3acd24bfbbca44b8b3cd56e2e3bacf97ef5648ab69b84989bf90",
    "sparse-0.19.2-py2.py3-none-any.whl": "2fcad7b83dacb5c70123a7de76db3626cbdaf567fcf12bdbf637aca6c4ee1bbe",
    "build/Bicrystal-4reference-el-pf-control-v1.py": "a1d3b971041c4e13aa288c72867690ce5d0e29265b11daab0e52418453d75db1",
    "build/functions_4ref_new.py": "9ce5fb1eec86498cc84cc6bed7c6e25dd37abcb87cecad096064b5a9e1acb78a",
    "build/control_build_manifest.json": "96a816a9423236a25a951cd5a3ef3f5dc097cb34010540846d52946bf20d8dc5",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checkpoint_step(path: Path) -> int:
    with np.load(path, allow_pickle=False) as archive:
        return int(archive["accepted_step"])


def prepare(base: Path, output: Path, start: int, target: int, energy: float,
            mobility: float, parent: Path, walltime: str = "12:00:00") -> dict:
    if start < 5_000 or target <= start or target > TERMINAL_STEP or start % 250 or target % 250:
        raise ValueError("production bounds must be 250-step aligned within 5000..200000")
    if checkpoint_step(parent) != start:
        raise ValueError("parent checkpoint step mismatch")
    for relative, expected in IDENTITIES.items():
        if sha256(base / relative) != expected:
            raise ValueError(f"base identity mismatch: {relative}")

    parent_sha = sha256(parent)
    parent_name = f"parent/{parent.name}"
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(base / "PF_Codes.zip", output / "PF_Codes.zip")
    shutil.copy2(base / "sparse-0.19.2-py2.py3-none-any.whl", output / "sparse-0.19.2-py2.py3-none-any.whl")
    shutil.copytree(base / "build", output / "build")
    (output / "parent").mkdir()
    shutil.copy2(parent, output / parent_name)

    wrapper = f'''#!/bin/bash
set -Eeuo pipefail
environment=/pub/sdillon1/sw/pfgg-anisotropy/f4aa53f-714dd16b2f3624b8
python_bin="$environment/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 QIU_NUMBA_THREADS=1 PYTHONHASHSEED=0 MPLBACKEND=Agg MPLCONFIGDIR="$PWD/mpl" PYTHONPATH="$PWD/third_party:$PWD/build"
mkdir -p output/checkpoints "$MPLCONFIGDIR" third_party
test "$(sha256sum {parent_name} | awk '{{print $1}}')" = {parent_sha}
export QIU_CONTROL_RESTART="$PWD/{parent_name}"
"$python_bin" -m pip install --no-deps --no-index --target third_party sparse-0.19.2-py2.py3-none-any.whl >output/pip.log 2>&1
cp build/control_build_manifest.json output/
export QIU_CONTROL_OUTPUT="$PWD/output/checkpoints" QIU_CONTROL_NAME={CONTROL} QIU_CONTROL_CODE=3 QIU_CONTROL_TARGET_STEP={target} QIU_CHECKPOINT_CADENCE=250 QIU_ENERGY_NORMALIZATION={energy:.17g} QIU_MOBILITY_NORMALIZATION={mobility:.17g}
/usr/bin/time -v -o output/time.txt "$python_bin" -u build/Bicrystal-4reference-el-pf-control-v1.py >output/stdout.log 2>output/stderr.log
"$python_bin" - <<'PY'
import hashlib,json
from pathlib import Path
import numpy as np
root=Path('output'); checkpoints=root/'checkpoints'; start={start}; target={target}
manifest=json.loads((checkpoints/'checkpoint_manifest.json').read_text())
assert manifest['logical_trajectory']=='{CONTROL}' and manifest['baseline_eligible'] is False and manifest['restart_used'] is True
assert all(item['finite'] for item in manifest['checkpoints'])
compilation=json.loads((checkpoints/'compilation_manifest.json').read_text()); assert compilation['object_mode_or_uncompiled']==[]
restart=checkpoints/f'native-restart-step{{target:06d}}.npz'; assert restart.is_file()
with np.load(restart,allow_pickle=False) as archive: assert int(archive['accepted_step'])==target
decision={{'schema':'qiu-aniso-production-segment-v1','classification':'QIU_SI_ANISO_PRODUCTION_SEGMENT_PASSED','logical_trajectory':'QIU_SI_4REF_ANISO_EM_INV_FULL','control':'{CONTROL}','segment_start':start,'segment_target':target,'parent_checkpoint_sha256':'{parent_sha}','end_checkpoint_sha256':hashlib.sha256(restart.read_bytes()).hexdigest(),'energy_normalization':{energy!r},'mobility_normalization':{mobility!r},'execution_threads':1,'finite':True}}
(root/'segment_decision.json').write_text(json.dumps(decision,indent=2)+'\n')
PY
find output -type f -print0|sort -z|xargs -0 sha256sum >output/payload-files.sha256
'''
    (output / "run_aniso_production_segment.sh").write_text(wrapper)
    (output / "run_aniso_production_segment.sh").chmod(0o755)
    (output / "job.toml").write_text(f'''name = "pfgg-qiu-aniso-em-full-{start:06d}-{target:06d}"
entrypoint = ["bash","run_aniso_production_segment.sh"]
working_directory = "."
modules = []
environment = {{ PYTHONHASHSEED = "0" }}
input_paths = ["run_aniso_production_segment.sh","PF_Codes.zip","sparse-0.19.2-py2.py3-none-any.whl","build","parent"]
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
events = ["FAIL","END"]
''')
    record = {
        "schema": "qiu-aniso-production-segment-preparation-v1", "control": CONTROL,
        "start": start, "target": target, "parent_checkpoint_sha256": parent_sha,
        "energy_normalization": energy, "mobility_normalization": mobility,
        "wrapper_sha256": sha256(output / "run_aniso_production_segment.sh"),
    }
    (output / "segment_preparation.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--start", type=int, required=True); parser.add_argument("--target", type=int, required=True)
    parser.add_argument("--energy-normalization", type=float, required=True)
    parser.add_argument("--mobility-normalization", type=float, required=True)
    parser.add_argument("--parent", type=Path, required=True); parser.add_argument("--walltime", default="12:00:00")
    args = parser.parse_args()
    print(json.dumps(prepare(args.base, args.output, args.start, args.target,
                             args.energy_normalization, args.mobility_normalization,
                             args.parent, args.walltime), indent=2))


if __name__ == "__main__":
    main()
