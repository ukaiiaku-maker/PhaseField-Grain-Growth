#!/usr/bin/env python3
"""Run one exact compact-support platform-confirmation branch."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy

from grain_growth_pf.pf.compact_support import exact_candidate_graph
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.qualify_compact_support import config, orientation_harmonics
from scripts.run_anisotropic_reduced_pilot import morphology


def sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b""): digest.update(block)
    return digest.hexdigest()


def array_sha(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("checkpoint",type=Path); parser.add_argument("output",type=Path)
    parser.add_argument("--expected-checkpoint-sha256",required=True); parser.add_argument("--dt",type=float,required=True)
    parser.add_argument("--kkt-tolerance",type=float,required=True); parser.add_argument("--steps",type=int,default=128)
    parser.add_argument("--scientific-source",required=True); parser.add_argument("--branch",required=True)
    args=parser.parse_args()
    if sha256(args.checkpoint)!=args.expected_checkpoint_sha256: raise SystemExit("checkpoint SHA-256 mismatch")
    if args.steps!=128: raise SystemExit("platform confirmation requires exactly 128 steps")
    args.output.mkdir(parents=True,exist_ok=False)
    with np.load(args.checkpoint,allow_pickle=False) as archive:
        metadata=json.loads(str(archive["checkpoint_state_json"])); state={key:archive[key].copy() for key in ("eta","active_phases","orientations","mobility_scale")}
    solver=MultiphaseFieldSolver(state["eta"],config(args.dt,args.kkt_tolerance),orientations=state["orientations"])
    solver.load_state_dict({**state,"time":metadata["time"],"step_number":metadata["step"]})
    initial_step=solver.step_number; initial_time=solver.time; history=[]
    for _ in range(args.steps):
        diag=solver.step(compute_energy=False); support=diag.compact_support
        history.append({"step":solver.step_number,"time":solver.time,"pre_step_energy":float(solver._last_pre_step_energy),
                        "kkt_residual":support.kkt_residual,"entries":support.entries,"retirements":support.retirements,
                        "line_search_scale":support.line_search_scale,"active_mean":support.active_mean,
                        "active_p95":support.active_p95,"active_max":support.active_max})
    final_energy=float(solver._anisotropic_energy()); graph,counts=exact_candidate_graph(solver.eta,True); exact=solver.eta>0
    state=solver.state_dict(); final_metadata={"schema":"compact-support-platform-checkpoint-v2","step":solver.step_number,"time":solver.time,
        "scientific_source":args.scientific_source,"parent_sha256":args.expected_checkpoint_sha256,"branch":args.branch}
    result_path=args.output/"final_state.npz"
    np.savez_compressed(result_path,eta=state["eta"],active_phases=state["active_phases"],orientations=state["orientations"],
                        mobility_scale=state["mobility_scale"],checkpoint_state_json=np.asarray(json.dumps(final_metadata,sort_keys=True)))
    result={"schema":"compact-support-platform-confirmation-v2","branch":args.branch,"scientific_source":args.scientific_source,
        "parent_checkpoint":str(args.checkpoint),"parent_sha256":args.expected_checkpoint_sha256,"dt":args.dt,
        "kkt_tolerance":args.kkt_tolerance,"steps_executed":args.steps,"initial_step":initial_step,"final_step":solver.step_number,
        "initial_time":initial_time,"final_time":solver.time,"eta_sha256":array_sha(solver.eta),
        "active_set_sha256":array_sha(exact),"candidate_graph_sha256":array_sha(graph),"candidate_counts_sha256":array_sha(counts),
        "orientations_sha256":array_sha(solver.orientations),"mobility_scale_sha256":array_sha(solver.mobility_scale),
        "final_energy":final_energy,"maximum_kkt_residual":max(x["kkt_residual"] for x in history),
        "history_sha256":hashlib.sha256(json.dumps(history,sort_keys=True,separators=(",",":")).encode()).hexdigest(),
        "morphology":morphology(solver.eta),"harmonics":orientation_harmonics(solver.eta,solver.orientations),
        "environment":{"python":platform.python_version(),"platform":platform.platform(),"numpy":np.__version__,"scipy":scipy.__version__},
        "final_state_sha256":sha256(result_path),"history":history}
    (args.output/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("branch","final_step","final_time","eta_sha256","final_state_sha256")},indent=2))


if __name__=="__main__": main()
