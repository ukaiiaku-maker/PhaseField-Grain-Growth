"""Numeric checkpoint instrumentation for the bounded native Qiu-SI preflight."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import resource
import time

import numba as nb
import numpy as np

try:
    from functions_4ref_new import E_elastic
except ModuleNotFoundError:
    # Static/unit-test imports do not stage the native runtime dependencies.
    E_elastic = None


CHECKPOINT_STEPS = (0, 1, 10, 100, 1000)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@nb.njit(cache=False)
def pair_diagnostics(phi, nf, mf, nx, ny, dx, dy, eta, ref_theta_i, sigma11, sigma12, sigma22, eij, misorientation_list):
    count = 0
    for l in range(nx):
        for m in range(ny):
            count += nf[l, m] * nf[l, m]
    cell_x = np.empty(count, np.int32); cell_y = np.empty(count, np.int32)
    phase_i = np.empty(count, np.int16); phase_j = np.empty(count, np.int16)
    elastic = np.empty(count, np.float64); capillary = np.empty(count, np.float64); barrier = np.empty(count, np.float64)
    beta_first = np.empty(count, np.float64); beta_second = np.empty(count, np.float64)
    reference_first = np.empty(count, np.float64); reference_second = np.empty(count, np.float64)
    cursor = 0
    for l in range(nx):
        l_p = (l + 1) % nx; l_m = (l - 1) % nx
        for m in range(ny):
            m_p = (m + 1) % ny; m_m = (m - 1) % ny
            for n1 in range(nf[l, m]):
                i = mf[n1, l, m]
                stencil_i = (phi[i,l_m,m_m] + 4*phi[i,l_m,m] + 4*phi[i,l_p,m] + phi[i,l_p,m_p] + phi[i,l_m,m_p] + 4*phi[i,l,m_m] + 4*phi[i,l,m_p] + phi[i,l_p,m_m] - 20*phi[i,l,m])/(6*dx*dx)
                for n2 in range(nf[l, m]):
                    j = mf[n2, l, m]
                    stencil_j = (phi[j,l_m,m_m] + 4*phi[j,l_m,m] + 4*phi[j,l_p,m] + phi[j,l_p,m_p] + phi[j,l_m,m_p] + 4*phi[j,l,m_m] + 4*phi[j,l,m_p] + phi[j,l_p,m_m] - 20*phi[j,l,m])/(6*dx*dx)
                    ppp = (phi[j,l,m]*stencil_i - phi[i,l,m]*stencil_j) + np.pi**2/(2*eta**2)*(phi[i,l,m]-phi[j,l,m])
                    # update_PF defines E_el=-E_elastic and then subtracts E_el;
                    # the signed force actually added to the pair drive is E_elastic.
                    signed_elastic = E_elastic(phi,i,j,l,m,l_p,l_m,m_p,m_m,ref_theta_i,misorientation_list[i],misorientation_list[j],sigma11[l,m],sigma12[l,m],sigma22[l,m],dx,dy)
                    gradient_i_x=(phi[i,l_p,m]-phi[i,l_m,m])/(2*dx); gradient_i_y=(phi[i,l,m_p]-phi[i,l,m_m])/(2*dy)
                    gradient_j_x=(phi[j,l_p,m]-phi[j,l_m,m])/(2*dx); gradient_j_y=(phi[j,l,m_p]-phi[j,l,m_m])/(2*dy)
                    normal_x=phi[i,l,m]*gradient_j_x-phi[j,l,m]*gradient_i_x; normal_y=phi[i,l,m]*gradient_j_y-phi[j,l,m]*gradient_i_y
                    mean_orientation=(misorientation_list[i]+misorientation_list[j])/2
                    rotated_x=np.cos(mean_orientation)*normal_x+np.sin(mean_orientation)*normal_y
                    rotated_y=-np.sin(mean_orientation)*normal_x+np.cos(mean_orientation)*normal_y
                    ratio=0.0
                    if rotated_y != 0.0: ratio=-rotated_x/rotated_y
                    if rotated_y != 0.0 and ratio >= 0.0 and ratio < 1.0: lower=0
                    elif rotated_y != 0.0 and ratio >= 1.0: lower=1
                    elif rotated_y == 0.0 or ratio < -1.0: lower=2
                    else: lower=3
                    upper=(lower+1)%4
                    theta=misorientation_list[i]-misorientation_list[j]
                    if theta >= np.pi/2 and theta < np.pi: theta-=np.pi/2
                    if theta >= np.pi and theta < 3*np.pi/2: theta-=np.pi
                    if theta >= 3*np.pi/2 and theta < 2*np.pi: theta-=3*np.pi/2
                    if theta <= -np.pi/2 and theta > -np.pi: theta+=np.pi/2
                    if theta <= -np.pi and theta > -3*np.pi/2: theta+=np.pi
                    if theta <= -3*np.pi/2 and theta > -2*np.pi: theta+=3*np.pi/2
                    b1=0.0
                    if abs(theta)<=36.7*np.pi/180: b1=2*np.tan(theta/2)
                    b2=-b1
                    cell_x[cursor]=l; cell_y[cursor]=m; phase_i[cursor]=i; phase_j[cursor]=j
                    elastic[cursor]=signed_elastic; capillary[cursor]=ppp; barrier[cursor]=np.pi/eta*np.sqrt(phi[i,l,m]*phi[j,l,m])*eij[i,j]
                    beta_first[cursor]=b1; beta_second[cursor]=b2
                    reference_first[cursor]=lower*ref_theta_i+mean_orientation; reference_second[cursor]=upper*ref_theta_i+mean_orientation
                    cursor += 1
    return cell_x,cell_y,phase_i,phase_j,elastic,capillary,barrier,beta_first,beta_second,reference_first,reference_second


def boundary_pairs(phi, overlap_pairs, criterion):
    """Return native ``index_gb`` order without changing its pair discovery."""
    return np.asarray(
        [(i, j) for i, j in overlap_pairs if np.any(np.abs(phi[i] * phi[j] - 0.25) < criterion)],
        dtype=np.int16,
    ).reshape((-1, 2))


def pbc_scalar(value: float, size: int) -> float:
    if value >= size / 2: return value - size
    if value <= -size / 2: return value + size
    return value


def line_records(gb_x, gb_y, boundary_pairs, orientations, ref_theta_i, nx, ny, dx, dy, beta_function):
    records=[]
    for boundary_id,(xs,ys) in enumerate(zip(gb_x,gb_y)):
        pair=boundary_pairs[boundary_id]; or1=orientations[pair[0]]; or2=orientations[pair[1]]; mor=(or1+or2)/2
        beta1,beta2=beta_function(or1,or2)
        for point in range(len(xs)):
            if len(xs)==1: ax,ay=0.0,0.0
            elif point==0: ax,ay=xs[1]-xs[0],ys[1]-ys[0]
            elif point==len(xs)-1: ax,ay=xs[point]-xs[point-1],ys[point]-ys[point-1]
            else: ax,ay=(xs[point+1]-xs[point-1])/2,(ys[point+1]-ys[point-1])/2
            dx_ds=-pbc_scalar(ax,nx)*dx; dy_ds=-pbc_scalar(ay,ny)*dy
            ddy=-(np.cos(mor)*dx_ds+np.sin(mor)*dy_ds); ddx=-np.sin(mor)*dx_ds+np.cos(mor)*dy_ds
            ratio=np.inf if ddx==0 else ddy/ddx
            if ddx!=0 and 0<=ratio<1: lower=0
            elif ddx!=0 and ratio>=1: lower=1
            elif ddx==0 or ratio < -1: lower=2
            else: lower=3
            upper=(lower+1)%4
            phi_k=lower*ref_theta_i+mor; phi_k1=(upper if lower<3 else 0)*ref_theta_i+mor
            d_x=(np.sin(phi_k1)*dx_ds-np.cos(phi_k1)*dy_ds)/np.sin(phi_k1-phi_k)
            d_y=(-np.sin(phi_k)*dx_ds+np.cos(phi_k)*dy_ds)/np.sin(phi_k1-phi_k)
            records.append((boundary_id,point,pair[0],pair[1],int(xs[point]),int(ys[point]),lower,upper,beta1,beta2,d_x,d_y))
    dtype=[('boundary_id','i4'),('point','i4'),('phase_i','i2'),('phase_j','i2'),('x','i4'),('y','i4'),('lower_reference','i1'),('upper_reference','i1'),('beta1','f8'),('beta2','f8'),('line_density_x','f8'),('line_density_y','f8')]
    return np.asarray(records,dtype=dtype)


def morphology(phi):
    labels=np.argmax(phi,axis=0); counts=np.bincount(labels.ravel(),minlength=phi.shape[0]); grain_count=int(np.count_nonzero(counts))
    edges=(labels!=np.roll(labels,1,0))|(labels!=np.roll(labels,1,1)); boundary_density=float(np.mean(edges))
    perimeters=np.zeros(phi.shape[0],dtype=np.int64)
    for shifted in (np.roll(labels,1,0),np.roll(labels,-1,0),np.roll(labels,1,1),np.roll(labels,-1,1)):
        diff=labels!=shifted; perimeters += np.bincount(labels[diff],minlength=phi.shape[0])
    compactness=np.zeros(phi.shape[0],dtype=float); valid=(counts>0)&(perimeters>0); compactness[valid]=4*np.pi*counts[valid]/perimeters[valid]**2
    return labels,counts,compactness,grain_count,boundary_density


def atomic_npz(path: Path, **arrays) -> None:
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('wb') as stream: np.savez_compressed(stream,**arrays)
    os.replace(temporary,path)


def atomic_json(path: Path, value) -> None:
    temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2)+'\n');os.replace(temporary,path)


def write_checkpoint(output, step, dt, phi, orientations, eij, pre_delta, post_delta, sigma11, sigma12, sigma22, gb_x, gb_y, boundary_pairs, ref_theta_i, nf_diag, mf_diag, dx, dy, eta, beta_function, initial_population=None, native_trial_buffer=None, trajectory_state=None):
    started=time.perf_counter();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    labels,counts,compactness,grain_count,boundary_density=morphology(phi); active=phi>0
    pairs=pair_diagnostics(phi,nf_diag,mf_diag,phi.shape[1],phi.shape[2],dx,dy,eta,ref_theta_i,sigma11,sigma12,sigma22,eij,orientations)
    lines=line_records(gb_x,gb_y,boundary_pairs,orientations,ref_theta_i,phi.shape[1],phi.shape[2],dx,dy,beta_function)
    if initial_population is None:
        maximum_population_loss = 0.0
    else:
        initial_population=np.asarray(initial_population)
        present=initial_population>0
        maximum_population_loss=float(np.max((initial_population[present]-counts[present])/initial_population[present]))
    summary={'step':step,'physical_time':step*dt,'grid_shape':[int(phi.shape[1]),int(phi.shape[2])],'cell_spacing':[dx,dy],'grain_count':grain_count,'boundary_density':boundary_density,'compactness_mean':float(np.mean(compactness[counts>0])),'maximum_population_loss':maximum_population_loss,'finite':bool(all(np.isfinite(x).all() for x in (phi,sigma11,sigma12,sigma22,pre_delta,post_delta))),'phase_sum_max_abs_error':float(np.max(np.abs(phi.sum(axis=0)-1))),'minimum_phase':float(phi.min()),'maximum_phase':float(phi.max()),'active_support_mean':float(active.sum(axis=0).mean()),'active_support_max':int(active.sum(axis=0).max()),'pair_record_count':int(len(pairs[0])),'line_record_count':int(len(lines)),'capillary_pair_term_sum':float(np.sum(pairs[5])),'signed_elastic_pair_force_sum':float(np.sum(pairs[4])),'barrier_pair_term_sum':float(np.sum(pairs[6])),'runtime_seconds':time.perf_counter()-started,'maximum_resident_kbytes':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)}
    checkpoint=output/f'native-state-step{step:06d}.npz'
    checkpoint_arrays=dict(accepted_step=np.asarray(step),physical_time=np.asarray(step*dt),output_scaling=np.asarray([dx,dy]),phi=phi,active_support=active,native_nf=nf_diag,native_mf_prefix=mf_diag,labels=labels,orientations=orientations,barrier=eij,pre_renormalization_delta=pre_delta,accepted_delta=post_delta,sigma11=sigma11,sigma12=sigma12,sigma22=sigma22,gb_records=lines,pair_cell_x=pairs[0],pair_cell_y=pairs[1],pair_phase_i=pairs[2],pair_phase_j=pairs[3],signed_elastic_pair_force=pairs[4],capillary_pair_term=pairs[5],barrier_pair_term=pairs[6],pair_beta_first=pairs[7],pair_beta_second=pairs[8],pair_reference_first=pairs[9],pair_reference_second=pairs[10],grain_population=counts,compactness=compactness,summary_json=np.asarray(json.dumps(summary,sort_keys=True)))
    if native_trial_buffer is not None: checkpoint_arrays['native_trial_buffer']=native_trial_buffer
    atomic_npz(checkpoint,**checkpoint_arrays)
    if trajectory_state is not None: atomic_npz(output/f'native-restart-step{step:06d}.npz',**trajectory_state)
    summary['path']=checkpoint.name;summary['sha256']=sha256(checkpoint);summary['bytes']=checkpoint.stat().st_size
    manifest_path=output/'checkpoint_manifest.json';manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {'schema':'qiu-native-preflight-v2-checkpoints','checkpoints':[]}
    manifest['checkpoints']=[x for x in manifest['checkpoints'] if x['step']!=step]+[summary];manifest['checkpoints'].sort(key=lambda x:x['step']);atomic_json(manifest_path,manifest)
    return summary
