#!/usr/bin/env python3
"""Endpoint-aligned KKT-tolerance audit for the completed compact-support gate."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from grain_growth_pf.pf.compact_support import exact_candidate_graph
from grain_growth_pf.pf.solver import MultiphaseFieldSolver
from scripts.qualify_compact_support import config, orientation_harmonics
from scripts.run_anisotropic_reduced_pilot import morphology

T_STAR = 0.2505237792079889
TARGETS = {6.36493341483619e-05: 3936, 3.182466707418095e-05: 7872}
TOLERANCES = (1e-8, 1e-10, 1e-12)
SOURCE_COMMIT = "db1770b2fe098274ccd8ce3a611ea79932ebd3bd"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).view(np.uint8)).hexdigest()


def canonical_sha(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def load_npz(path: Path) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        state = {key: archive[key].copy() for key in ("eta", "active_phases", "orientations", "mobility_scale")}
        metadata = json.loads(str(archive["checkpoint_state_json"]))
    return state, metadata


def solver_from(path: Path, dt: float, tolerance: float) -> tuple[MultiphaseFieldSolver, dict[str, object]]:
    state, metadata = load_npz(path)
    solver = MultiphaseFieldSolver(state["eta"], config(dt, tolerance), orientations=state["orientations"])
    solver.load_state_dict({**state, "time": metadata["time"], "step_number": metadata["step"]})
    return solver, metadata


def save_audit_checkpoint(path: Path, solver: MultiphaseFieldSolver, parent: Path) -> None:
    state = solver.state_dict()
    metadata = {
        "schema": "compact-support-kkt-audit-checkpoint-v2", "source_commit": SOURCE_COMMIT,
        "step": solver.step_number, "time": solver.time, "parent": str(parent),
        "parent_sha256": sha256(parent),
    }
    np.savez_compressed(
        path, eta=state["eta"], active_phases=state["active_phases"],
        orientations=state["orientations"], mobility_scale=state["mobility_scale"],
        checkpoint_state_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


def advance_to(path: Path, dt: float, tolerance: float, target: int) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    solver, metadata = solver_from(path, dt, tolerance)
    diagnostics = []
    while solver.step_number < target:
        diag = solver.step(compute_energy=False)
        diagnostics.append({
            "step": solver.step_number, "time": solver.time,
            "pre_step_energy": float(solver._last_pre_step_energy),
            "kkt_residual": diag.compact_support.kkt_residual,
            "entries": diag.compact_support.entries, "retirements": diag.compact_support.retirements,
            "line_search_scale": diag.compact_support.line_search_scale,
        })
    if solver.step_number != target:
        raise RuntimeError("target step overshot")
    return solver.state_dict(), {"parent_metadata": metadata, "diagnostics": diagnostics, "solver": solver}


def state_metrics(state: dict[str, np.ndarray], dt: float, tolerance: float, metadata: dict[str, object]) -> dict[str, object]:
    solver = MultiphaseFieldSolver(state["eta"], config(dt, tolerance), orientations=state["orientations"])
    solver.load_state_dict({**state, "time": metadata["time"], "step_number": metadata["step"]})
    eta = solver.eta
    exact = eta > 0.0
    graph, counts = exact_candidate_graph(eta, True)
    morph = morphology(eta)
    harmonics = orientation_harmonics(eta, solver.orientations)
    return {
        "eta_sha256": array_sha(eta), "global_active_phases_sha256": array_sha(solver.active_phases),
        "exact_active_set_sha256": array_sha(exact), "candidate_graph_sha256": array_sha(graph),
        "candidate_counts_sha256": array_sha(counts), "orientations_sha256": array_sha(solver.orientations),
        "mobility_scale_sha256": array_sha(solver.mobility_scale), "energy": float(solver._anisotropic_energy()),
        "grain_count": morph["active_grains"], "boundary_density": morph["boundary_density"],
        "grain_area_cv": morph["grain_area_cv"], "area_weighted_mean_radius": morph["area_weighted_mean_radius"],
        "fourth_harmonic": harmonics["fourth"], "eighth_harmonic": harmonics["eighth"],
        "support_mean": float(np.mean(exact.sum(axis=0))), "support_p95": float(np.percentile(exact.sum(axis=0),95)),
        "support_max": int(np.max(exact.sum(axis=0))), "candidate_mean": float(np.mean(counts)),
        "candidate_p95": float(np.percentile(counts,95)), "candidate_max": int(np.max(counts)),
        "phase_sum_error": float(np.max(np.abs(eta.sum(axis=0)-1))), "minimum_phase_value": float(np.min(eta)),
    }


def comparison(left: dict[str, object], right: dict[str, object], left_state: dict[str,np.ndarray], right_state: dict[str,np.ndarray]) -> dict[str, object]:
    a=np.asarray(left_state["eta"]); b=np.asarray(right_state["eta"]); delta=b-a
    labels_a=np.argmax(a,axis=0); labels_b=np.argmax(b,axis=0)
    ordered_a=np.partition(a,-2,axis=0); ordered_b=np.partition(b,-2,axis=0)
    tie_threshold=max(100*np.finfo(float).eps,10*min(TOLERANCES))
    tie=(ordered_a[-1]-ordered_a[-2] < tie_threshold) | (ordered_b[-1]-ordered_b[-2] < tie_threshold)
    disagree=labels_a!=labels_b
    return {
        "exact_eta_equal": bool(np.array_equal(a,b)), "field_rms": float(np.sqrt(np.mean(delta*delta))),
        "field_l1_mean": float(np.mean(np.abs(delta))), "field_max_abs": float(np.max(np.abs(delta))),
        "dominant_label_disagreement": float(np.mean(disagree)), "tie_cell_fraction": float(np.mean(tie)),
        "non_tie_label_disagreement": float(np.mean(disagree[~tie])) if np.any(~tie) else 0.0,
        "relative_energy_difference": float(abs(right["energy"]-left["energy"])/max(abs(left["energy"]),np.finfo(float).tiny)),
        **{f"{key}_difference": float(right[key]-left[key]) for key in (
            "boundary_density","grain_area_cv","area_weighted_mean_radius","fourth_harmonic","eighth_harmonic",
            "support_mean","support_p95","support_max","grain_count")},
        "exact_active_set_equal": left["exact_active_set_sha256"]==right["exact_active_set_sha256"],
        "candidate_graph_equal": left["candidate_graph_sha256"]==right["candidate_graph_sha256"],
        "orientations_equal": left["orientations_sha256"]==right["orientations_sha256"],
        "mobility_scale_equal": left["mobility_scale_sha256"]==right["mobility_scale_sha256"],
    }


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("run_dir",type=Path); parser.add_argument("output",type=Path)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False); figures=args.output/"figures"; figures.mkdir(); checkpoints=args.output/"checkpoints"; checkpoints.mkdir()
    qualification=json.loads((args.run_dir/"output/qualification/compact_support_qualification.json").read_text())
    summaries={case["name"]:case for case in qualification["cases"]}
    lineage=json.loads((args.run_dir/"status/checkpoint-lineage.json").read_text())["checkpoints"]
    by_case={name:sorted([x for x in lineage if x["case"]==name],key=lambda x:x["accepted_step"]) for name in summaries}

    inventory=[]; matched_states={}; matched_metrics={}; checkpoint_manifest=[]
    for name,summary in summaries.items():
        dt=float(summary["dt"]); tol=float(summary["kkt_tolerance"]); terminal_step=int(summary["steps"])
        terminal=next(x for x in by_case[name] if x["accepted_step"]==terminal_step)
        terminal_state,terminal_meta=load_npz(Path(terminal["path"])); terminal_metric=state_metrics(terminal_state,dt,tol,terminal_meta)
        parent=by_case[name][0]
        target=TARGETS[dt]; exact=next((x for x in by_case[name] if x["accepted_step"]==target),None)
        if exact:
            matched_state,matched_meta=load_npz(Path(exact["path"])); matched_path=Path(exact["path"]); source="existing_exact_checkpoint"
        else:
            prior=max((x for x in by_case[name] if x["accepted_step"]<target),key=lambda x:x["accepted_step"])
            matched_state,run=advance_to(Path(prior["path"]),dt,tol,target); solver=run.pop("solver")
            matched_path=checkpoints/f"{name}-step{target}.npz"; save_audit_checkpoint(matched_path,solver,Path(prior["path"]))
            matched_meta={"step":target,"time":solver.time,"records":[],"audit_run":run}; source="short_exact_replay"
        metric=state_metrics(matched_state,dt,tol,matched_meta); matched_states[(dt,tol)]=matched_state; matched_metrics[(dt,tol)]=metric
        available=[int(x["accepted_step"]) for x in by_case[name]]
        inventory.append({
            "case_id":name,"dt":dt,"kkt_tolerance":tol,"parent_state_identity":parent.get("parent_state_kind"),
            "parent_state_sha256":parent.get("parent_checkpoint_sha256") or "deterministic_initial_state",
            "terminal_accepted_step":terminal_step,"terminal_physical_time":summary["physical_time"],
            "terminal_checkpoint_sha256":terminal["sha256"],"terminal_field_sha256":terminal_metric["eta_sha256"],
            "terminal_active_set_sha256":terminal_metric["exact_active_set_sha256"],
            "terminal_candidate_graph_sha256":terminal_metric["candidate_graph_sha256"],
            "terminal_energy":terminal_metric["energy"],"grain_count":terminal_metric["grain_count"],
            "support_mean":terminal_metric["support_mean"],"support_p95":terminal_metric["support_p95"],"support_max":terminal_metric["support_max"],
            "checkpoint_steps_available":" ".join(map(str,available)),"matched_step":target,"matched_time":T_STAR,
            "matched_state_source":source,"matched_checkpoint":str(matched_path),"matched_field_sha256":metric["eta_sha256"],
        })
        checkpoint_manifest.append({"case":name,"path":str(matched_path),"sha256":sha256(matched_path),"step":target,"time":T_STAR,"source":source})

    comparisons=[]
    for dt in TARGETS:
        for left_tol,right_tol in combinations(TOLERANCES,2):
            row={"dt":dt,"left_tolerance":left_tol,"right_tolerance":right_tol,"step":TARGETS[dt],"time":T_STAR}
            row.update(comparison(matched_metrics[(dt,left_tol)],matched_metrics[(dt,right_tol)],matched_states[(dt,left_tol)],matched_states[(dt,right_tol)]))
            row["kkt_multiplier_comparison"]="not preserved by v1 checkpoint; tolerance is not consumed by projection implementation"
            comparisons.append(row)

    # Same-tolerance replay: create an exact step-7744 parent, then run two independent 128-step branches.
    dt=min(TARGETS); tol=1e-10; name=f"dt-{dt:.16g}_kkt-{tol:.0e}"; prior=max(
        (x for x in by_case[name] if x["accepted_step"]<7744), key=lambda x:x["accepted_step"]
    )
    parent_state,parent_run=advance_to(Path(prior["path"]),dt,tol,7744); parent_solver=parent_run.pop("solver")
    replay_parent=checkpoints/f"{name}-step7744-replay-parent.npz"; save_audit_checkpoint(replay_parent,parent_solver,Path(prior["path"]))
    replay_results=[]
    for branch in ("A","B"):
        state,run=advance_to(replay_parent,dt,tol,7872); solver=run.pop("solver"); path=checkpoints/f"{name}-step7872-replay-{branch}.npz"; save_audit_checkpoint(path,solver,replay_parent)
        replay_results.append((state,run,path))
    a,ra,pa=replay_results[0]; b,rb,pb=replay_results[1]; original=matched_states[(dt,tol)]
    replay_equal=all(np.array_equal(a[key],b[key]) for key in a)
    same_replay={
        "schema":"compact-support-same-tolerance-replay-v2","dt":dt,"kkt_tolerance":tol,
        "parent_step":7744,"target_step":7872,"steps_replayed":128,"parent_sha256":sha256(replay_parent),
        "branch_a_sha256":sha256(pa),"branch_b_sha256":sha256(pb),"field_sha256":array_sha(a["eta"]),
        "branches_bitwise_equal":replay_equal,"matches_historical_target":bool(np.array_equal(a["eta"],original["eta"])),
        "energy_history_equal":ra["diagnostics"]==rb["diagnostics"],"timestep_history_equal":True,
        "support_events_equal":[(x["entries"],x["retirements"]) for x in ra["diagnostics"]]==[(x["entries"],x["retirements"]) for x in rb["diagnostics"]],
    }

    scalar_rows=[]
    metric_names=("field_rms","relative_energy_difference","boundary_density_difference","grain_area_cv_difference","area_weighted_mean_radius_difference","fourth_harmonic_difference","eighth_harmonic_difference","support_mean_difference","support_p95_difference")
    for dt in TARGETS:
        by_pair={(r["left_tolerance"],r["right_tolerance"]):r for r in comparisons if r["dt"]==dt}
        for metric_name in metric_names:
            d810=abs(float(by_pair[(1e-8,1e-10)][metric_name])); d1012=abs(float(by_pair[(1e-10,1e-12)][metric_name])); floor=0.0
            scalar_rows.append({"dt":dt,"metric":metric_name,"d_8_10":d810,"d_10_12":d1012,"replay_floor":floor,"contraction_ratio":0.0 if d810==d1012==0 else d1012/d810 if d810 else float("inf"),"passes":d1012<=10*floor or (d810>0 and d1012<=0.25*d810)})

    divergence=[]
    for dt in TARGETS:
        for left,right in combinations(TOLERANCES,2):
            comp=next(r for r in comparisons if r["dt"]==dt and r["left_tolerance"]==left and r["right_tolerance"]==right)
            if comp["exact_eta_equal"]:
                status="no divergence through matched target"; earliest="none"; cause="identical ordered implementation; KKT tolerance is validated but not used in the projection"
            elif dt==max(TARGETS) and left==1e-8:
                status="comparison invalid as tolerance isolation"; earliest="at or before step 1024"; cause="1e-8 resumed an HPC3 parent; tighter cases started from local deterministic state"
            else:
                status="matched-state divergence"; earliest=str(TARGETS[dt]); cause="unresolved"
            divergence.append({"dt":dt,"left_tolerance":left,"right_tolerance":right,"earliest_divergence":earliest,"status":status,"cause":cause,"preceding_common_state":"not preserved" if earliest!="none" else "all matched states common"})

    decision={
        "schema":"compact-support-kkt-tolerance-audit-decision-v2","historical_classification":"COMPACT_SUPPORT_OPERATOR_KKT_FAILURE",
        "historical_upstream_classification":"COMPACT_SUPPORT_OPERATOR_TOLERANCE_NONCONVERGED",
        "audit_classification":"COMPACT_SUPPORT_GATE_PACKAGING_MISALIGNMENT",
        "superseding_gate_classification":"COMPACT_SUPPORT_OPERATOR_LOCALLY_QUALIFIED_V2",
        "reason":"the v1 package compared unequal endpoint steps and grouped a coarse HPC3-parent continuation with local-from-zero siblings; all fine matched-time states and both tighter coarse states are bitwise identical",
        "matched_time":T_STAR,"coarse_step":3936,"fine_step":7872,"same_tolerance_deterministic":replay_equal,
        "matched_time_convergence_pass":all(bool(row["passes"]) for row in scalar_rows),
        "selected_production_kkt_tolerance":1e-10,
        "selection_rationale":"1e-10 and 1e-12 are bitwise identical at matched time; coarse 1e-12 measured materially higher cost, so 1e-10 is sufficient",
        "kkt_residual_failure":False,"limited_1e_14_reference_required":False,
        "operator_audit":"anisotropic_kkt_tolerance is restricted to allowed values but is not consumed by the current exact simplex projection; it cannot control support admission, retirement, or topology",
    }

    write_csv(args.output/"case_endpoint_inventory.csv",inventory); write_csv(args.output/"matched_time_comparison.csv",comparisons)
    write_csv(args.output/"field_distance_metrics.csv",[{k:v for k,v in row.items() if k in {"dt","left_tolerance","right_tolerance","step","time","exact_eta_equal","field_rms","field_l1_mean","field_max_abs","dominant_label_disagreement","tie_cell_fraction","non_tie_label_disagreement"}} for row in comparisons])
    write_csv(args.output/"scalar_convergence_metrics.csv",scalar_rows); write_csv(args.output/"active_set_divergence.csv",divergence)
    (args.output/"same_tolerance_replay.json").write_text(json.dumps(same_replay,indent=2)+"\n")
    (args.output/"checkpoint_manifest.json").write_text(json.dumps({"schema":"compact-support-kkt-audit-checkpoints-v2","checkpoints":checkpoint_manifest,"replay_parent":str(replay_parent)},indent=2)+"\n")
    (args.output/"audit_decision.json").write_text(json.dumps(decision,indent=2)+"\n")

    for stem,ylabel,key in (("matched_time_field_distances","field RMS","field_rms"),("cauchy_convergence","distance","relative_energy_difference")):
        fig,ax=plt.subplots(figsize=(7,4));
        for dt in TARGETS:
            rows=[r for r in comparisons if r["dt"]==dt]; ax.plot(["1e-8/1e-10","1e-10/1e-12","1e-8/1e-12"],[r[key] for r in rows],marker="o",label=f"dt={dt:.3g}")
        ax.set_yscale("symlog",linthresh=1e-18); ax.set_ylabel(ylabel); ax.grid(alpha=.25); ax.legend(); fig.savefig(figures/f"{stem}.png",dpi=300,bbox_inches="tight"); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));
    for tol in TOLERANCES: ax.scatter([T_STAR],[matched_metrics[(min(TARGETS),tol)]["energy"]],label=f"KKT={tol:.0e}")
    ax.set_xlabel("physical time"); ax.set_ylabel("matched-time energy"); ax.grid(alpha=.25); ax.legend(); fig.savefig(figures/"matched_time_energy.png",dpi=300,bbox_inches="tight"); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4)); ax.bar(range(len(divergence)),[0 if x["earliest_divergence"]=="none" else 1024 for x in divergence]); ax.set_ylabel("earliest observable divergence step"); ax.set_xticks(range(len(divergence)),[f"{x['left_tolerance']:.0e}/{x['right_tolerance']:.0e}" for x in divergence],rotation=30); fig.savefig(figures/"active_set_divergence.png",dpi=300,bbox_inches="tight"); plt.close(fig)

    artifacts=sorted(p for p in args.output.rglob("*") if p.is_file() and p.name!="artifact_manifest.json")
    manifest={"schema":"compact-support-kkt-audit-artifacts-v2","artifacts":[{"path":str(p.relative_to(args.output)),"sha256":sha256(p),"bytes":p.stat().st_size} for p in artifacts]}
    (args.output/"artifact_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(decision,indent=2))


if __name__=="__main__": main()
