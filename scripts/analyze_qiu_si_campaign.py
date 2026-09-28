#!/usr/bin/env python3
"""Postprocess verified native and anisotropic Qiu-SI checkpoint archives."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb
from PIL import Image


ORIENTATION_BINS = 36
ETA = 0.8
DT = 0.1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows_csv(path: Path, fields: Iterable[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def state_members(archive: Path) -> list[tuple[int, str]]:
    with tarfile.open(archive, "r:gz") as stream:
        found = []
        for name in stream.getnames():
            match = __import__("re").search(r"native-state-step(\d+)\.npz$", name)
            if match:
                found.append((int(match.group(1)), name))
    return sorted(found)


def load_npz(stream: tarfile.TarFile, member: str) -> dict[str, np.ndarray]:
    payload = stream.extractfile(member).read()
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        return {key: np.asarray(archive[key]) for key in archive.files}


def interface_statistics(labels: np.ndarray, orientations: np.ndarray) -> tuple[np.ndarray, dict[int, float], float]:
    angles=[]; weights=[]
    harmonic_numerator={2:0j,4:0j,6:0j,8:0j}; total=0.0
    energy=0.0
    for axis, normal in ((0, 0.0), (1, math.pi/2)):
        other=np.roll(labels,-1,axis=axis); mask=labels!=other
        first=labels[mask]; second=other[mask]
        if first.size == 0: continue
        mean=0.5*(orientations[first]+orientations[second])
        relative=np.mod(normal-mean,np.pi)
        delta=np.abs((orientations[first]-orientations[second]+np.pi/4)%(np.pi/2)-np.pi/4)
        misorientation=0.65+0.35*np.sin(2*delta)**2
        # Sharp-interface counterpart of the frozen A2 support law.
        a=np.abs(np.cos(relative))**16+np.abs(np.sin(relative))**16
        support=a**(1/16)
        gamma=1.0910512514090829*misorientation*(0.15+0.85*support)
        energy += float(np.sum(gamma)); angles.append(relative); weights.append(gamma)
        total += float(np.sum(gamma))
        for order in harmonic_numerator:
            harmonic_numerator[order] += np.sum(gamma*np.exp(1j*order*relative))
    if not angles:
        return np.zeros(ORIENTATION_BINS), {order:0.0 for order in harmonic_numerator}, 0.0
    angle=np.concatenate(angles); weight=np.concatenate(weights); edges=np.linspace(0,np.pi,ORIENTATION_BINS+1)
    histogram,_=np.histogram(angle,bins=edges,weights=weight)
    density=histogram/(max(float(histogram.sum()),np.finfo(float).tiny)*np.diff(edges))
    harmonics={order:float(abs(value)/total) for order,value in harmonic_numerator.items()}
    return density,harmonics,energy


def work_terms(state: dict[str, np.ndarray]) -> tuple[float,float,float]:
    x=state["pair_cell_x"].astype(int); y=state["pair_cell_y"].astype(int)
    i=state["pair_phase_i"].astype(int); j=state["pair_phase_j"].astype(int)
    delta=state["accepted_delta"]
    difference=delta[i,x,y]-delta[j,x,y]
    elastic=float(np.sum(state["signed_elastic_pair_force"]*difference))
    barrier=float(np.sum(state["barrier_pair_term"]*difference))
    capillary=float(np.sum(state["capillary_pair_term"]*difference))
    return capillary,elastic,barrier


def diffuse_energy_proxy(phi: np.ndarray, spacing: np.ndarray) -> tuple[float,float]:
    dx=float(spacing[0]); dy=float(spacing[1])
    gx=(np.roll(phi,-1,axis=1)-phi)/dx; gy=(np.roll(phi,-1,axis=2)-phi)/dy
    gradient=0.25*float(np.sum(gx*gx+gy*gy))*dx*dy
    obstacle=(np.pi**2/(4*ETA**2))*float(np.sum(phi*(1-phi)))*dx*dy
    return gradient,obstacle


def render(state: dict[str,np.ndarray]) -> Image.Image:
    labels=state["labels"].astype(int); orientations=state["orientations"]
    hue=np.mod(orientations[labels],np.pi)/np.pi
    hsv=np.stack((hue,np.full_like(hue,.72),np.full_like(hue,.92)),axis=-1)
    rgb=(255*hsv_to_rgb(hsv)).astype(np.uint8)
    boundary=(labels!=np.roll(labels,1,0))|(labels!=np.roll(labels,1,1));rgb[boundary]=18
    return Image.fromarray(rgb)


def configure() -> None:
    mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,
                         "axes.spines.right":False,"figure.dpi":120,"savefig.dpi":300})


def save(fig, output: Path, stem: str) -> list[str]:
    paths=[]
    for suffix in ("png","pdf"):
        path=output/f"{stem}.{suffix}";fig.savefig(path,bbox_inches="tight",facecolor="white");paths.append(path.name)
    plt.close(fig);return paths


def analyze(manifest_path: Path, output: Path) -> dict:
    manifest=json.loads(manifest_path.read_text());output.mkdir(parents=True,exist_ok=True)
    (output/"figures").mkdir(exist_ok=True);(output/"movies").mkdir(exist_ok=True)
    configure(); metrics=[]; normals=[]; harmonics=[]; final_images={}; movie_frames=defaultdict(list)
    archive_rows=[]; seen=defaultdict(set)
    for case in manifest["cases"]:
        name=case["case"]
        for archive_name in case["archives"]:
            archive=Path(archive_name); archive_rows.append({"case":name,"archive":str(archive),"sha256":sha256(archive),"bytes":archive.stat().st_size})
            members=state_members(archive)
            with tarfile.open(archive,"r:gz") as stream:
                for step,member in members:
                    if step in seen[name]: continue
                    seen[name].add(step);state=load_npz(stream,member);summary=json.loads(str(state["summary_json"]))
                    phi=state["phi"]; spacing=state["output_scaling"]; gradient,obstacle=diffuse_energy_proxy(phi,spacing)
                    density,amps,sharp_energy=interface_statistics(state["labels"],state["orientations"])
                    cap_work,elastic_work,barrier_work=work_terms(state)
                    accepted=state["accepted_delta"]
                    row={"case":name,"step":step,"physical_time":float(state["physical_time"]),
                         "grain_count":summary["grain_count"],"boundary_density":summary["boundary_density"],
                         "compactness_mean":summary["compactness_mean"],"active_support_mean":summary["active_support_mean"],
                         "active_support_max":summary["active_support_max"],"phase_sum_max_abs_error":summary["phase_sum_max_abs_error"],
                         "stress11_rms":float(np.sqrt(np.mean(state["sigma11"]**2))),"stress12_rms":float(np.sqrt(np.mean(state["sigma12"]**2))),
                         "stress22_rms":float(np.sqrt(np.mean(state["sigma22"]**2))),"stress_max_abs":float(max(np.max(np.abs(state[x])) for x in ("sigma11","sigma12","sigma22"))),
                         "line_density_l1":float(np.sum(np.abs(state["gb_records"]["line_density_x"])+np.abs(state["gb_records"]["line_density_y"]))),
                         "beta_rms":float(np.sqrt(np.mean(state["pair_beta_first"]**2))) if state["pair_beta_first"].size else 0.0,
                         "reference_count":int(np.unique(np.concatenate((state["pair_reference_first"],state["pair_reference_second"]))).size),
                         "diffuse_gradient_energy_proxy":gradient,"diffuse_obstacle_energy_proxy":obstacle,
                         "sharp_anisotropic_interface_energy_proxy":sharp_energy,
                         "capillary_increment_work":cap_work,"elastic_increment_work":elastic_work,"barrier_increment_work":barrier_work,
                         "accepted_increment_l2":float(np.sqrt(np.sum(accepted*accepted))),
                         "dissipation_proxy":float(np.sum(accepted*accepted)/DT)}
                    row["drive_work_sum"]=cap_work+elastic_work+barrier_work;metrics.append(row)
                    centers=np.degrees((np.arange(ORIENTATION_BINS)+.5)*np.pi/ORIENTATION_BINS)
                    normals.extend({"case":name,"step":step,"angle_degrees":float(angle),"density_per_radian":float(value)} for angle,value in zip(centers,density))
                    harmonics.extend({"case":name,"step":step,"order":order,"magnitude":value} for order,value in amps.items())
                    image=render(state);final_images[name]=(step,image.copy())
                    if case.get("movie",False): movie_frames[name].append((step,image.copy()))
    metrics.sort(key=lambda row:(row["case"],row["step"]));normals.sort(key=lambda row:(row["case"],row["step"],row["angle_degrees"]));harmonics.sort(key=lambda row:(row["case"],row["step"],row["order"]))
    metric_fields=list(metrics[0])
    rows_csv(output/"checkpoint_metrics.csv",metric_fields,metrics);rows_csv(output/"interface_normal_distributions.csv",normals[0].keys(),normals);rows_csv(output/"anisotropic_harmonics.csv",harmonics[0].keys(),harmonics);rows_csv(output/"archive_checksums.csv",archive_rows[0].keys(),archive_rows)

    by_case=defaultdict(list)
    for row in metrics: by_case[row["case"]].append(row)
    native={r["step"]:r for r in by_case["QIU_SI_NATIVE"]};aniso={r["step"]:r for r in by_case["QIU_SI_ANISO_EM_INV"]}
    compare=[]
    for step in sorted(native.keys()&aniso.keys()):
        a=native[step];b=aniso[step]
        compare.append({"step":step,"physical_time":a["physical_time"],"grain_count_delta":b["grain_count"]-a["grain_count"],
                        "boundary_density_delta":b["boundary_density"]-a["boundary_density"],"compactness_delta":b["compactness_mean"]-a["compactness_mean"],
                        "stress_rms_ratio":b["stress12_rms"]/max(a["stress12_rms"],np.finfo(float).tiny),
                        "line_density_ratio":b["line_density_l1"]/max(a["line_density_l1"],np.finfo(float).tiny),
                        "interface_energy_proxy_ratio":b["sharp_anisotropic_interface_energy_proxy"]/max(a["sharp_anisotropic_interface_energy_proxy"],np.finfo(float).tiny)})
    if compare: rows_csv(output/"matched_time_native_anisotropic.csv",compare[0].keys(),compare)

    colors=plt.cm.tab10(np.linspace(0,1,max(1,len(by_case))))
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
    for color,(name,values) in zip(colors,sorted(by_case.items())):
        t=[r["physical_time"] for r in values];axes[0,0].plot(t,[r["grain_count"] for r in values],label=name,color=color);axes[0,1].plot(t,[r["boundary_density"] for r in values],color=color);axes[1,0].plot(t,[r["compactness_mean"] for r in values],color=color);axes[1,1].plot(t,[r["active_support_mean"] for r in values],color=color)
    for ax,label in zip(axes.flat,("Grain count","Boundary density","Mean compactness","Mean active support")):ax.set_ylabel(label);ax.grid(alpha=.2)
    axes[1,0].set_xlabel("Physical time");axes[1,1].set_xlabel("Physical time");axes[0,0].legend(fontsize=6)
    figures=save(fig,output/"figures","kinetics_morphology")
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for color,(name,values) in zip(colors,sorted(by_case.items())):
        t=[r["physical_time"] for r in values];axes[0].plot(t,[r["stress12_rms"] for r in values],label=name,color=color);axes[1].plot(t,[r["line_density_l1"] for r in values],label=name,color=color)
    axes[0].set(xlabel="Physical time",ylabel="RMS shear stress");axes[1].set(xlabel="Physical time",ylabel="Line-density L1");axes[0].legend(fontsize=6)
    figures+=save(fig,output/"figures","stress_line_density")
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for color,(name,values) in zip(colors,sorted(by_case.items())):
        t=[r["physical_time"] for r in values];axes[0].plot(t,[r["diffuse_gradient_energy_proxy"]+r["diffuse_obstacle_energy_proxy"] for r in values],label=name,color=color);axes[1].plot(t,np.cumsum([r["drive_work_sum"] for r in values]),label=name,color=color)
    axes[0].set(xlabel="Physical time",ylabel="Diffuse energy proxy");axes[1].set(xlabel="Physical time",ylabel="Cumulative resolved drive work");axes[0].legend(fontsize=6)
    figures+=save(fig,output/"figures","energy_work_balance")
    final_normal={}
    for name,values in by_case.items():
        step=values[-1]["step"];final_normal[name]=[row for row in normals if row["case"]==name and row["step"]==step]
    fig,ax=plt.subplots(figsize=(7,4))
    for name,values in sorted(final_normal.items()):ax.plot([r["angle_degrees"] for r in values],[r["density_per_radian"] for r in values],label=name)
    ax.set(xlabel="Crystal-frame interface normal (degrees)",ylabel="Weighted density");ax.legend(fontsize=6)
    figures+=save(fig,output/"figures","interface_normal_distribution")
    names=sorted(final_images);fig,axes=plt.subplots(1,len(names),figsize=(3*len(names),3),squeeze=False)
    for ax,name in zip(axes[0],names):step,image=final_images[name];ax.imshow(image);ax.set_title(f"{name}\nstep {step}");ax.axis("off")
    figures+=save(fig,output/"figures","terminal_morphology_contact_sheet")
    movies=[]
    for name,frames in movie_frames.items():
        frames.sort();images=[image.resize((500,500),Image.Resampling.NEAREST) for _,image in frames]
        path=output/"movies"/f"{name.lower()}.gif";images[0].save(path,save_all=True,append_images=images[1:],duration=80,loop=0,optimize=False);movies.append(path.name)
    energy_balance={}
    for name,values in by_case.items():
        total=np.asarray([r["diffuse_gradient_energy_proxy"]+r["diffuse_obstacle_energy_proxy"] for r in values]);work=np.asarray([r["drive_work_sum"] for r in values])
        residual=np.diff(total)-work[1:]
        energy_balance[name]={"initial_energy_proxy":float(total[0]),"final_energy_proxy":float(total[-1]),"cumulative_drive_work":float(work.sum()),"balance_residual_l2":float(np.sqrt(np.sum(residual*residual))),"states":len(values)}
    summary={"schema":"qiu-si-campaign-analysis-v1","classification":"QIU_SI_FINAL_ANALYSIS_COMPLETE","manifest_sha256":sha256(manifest_path),"cases":{name:{"first_step":values[0]["step"],"last_step":values[-1]["step"],"states":len(values)} for name,values in by_case.items()},"energy_work_accounting":energy_balance,"matched_native_anisotropic_states":len(compare),"figures":figures,"movies":movies,"limitations":["energy columns are explicitly reported discrete proxies because the archived Qiu driver does not expose a scalar thermodynamic potential","single deterministic native network; no ensemble uncertainty interval"]}
    (output/"analysis_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    readme=f"""# Qiu-SI native and anisotropic campaign analysis\n\nClassification: `{summary['classification']}`\n\nThe tables contain checkpoint-resolved kinetics, morphology, support, stress, line density, beta/reference, interface-normal harmonics, and energy/work diagnostics. Energy quantities are labeled proxies in both column names and the machine-readable summary. The matched-time table contains the native-versus-combined comparison at every shared checkpoint.\n"""
    (output/"README.md").write_text(readme)
    generated=[path for path in output.rglob("*") if path.is_file() and path.name!="checksums.sha256"]
    (output/"checksums.sha256").write_text("".join(f"{sha256(path)}  {path.relative_to(output)}\n" for path in sorted(generated)))
    return summary


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("manifest",type=Path);parser.add_argument("output",type=Path);args=parser.parse_args()
    print(json.dumps(analyze(args.manifest,args.output),indent=2))


if __name__=="__main__":main()
