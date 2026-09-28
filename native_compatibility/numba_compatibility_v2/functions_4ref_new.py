import numpy as np
from numpy.random import *
import matplotlib
import matplotlib.pyplot as plt
import math
import random
import os
import numba as nb
from scipy import signal
from numba.typed import List
import time
import random
import sparse



@nb.njit(fastmath=False, parallel=True)
def Initialize(phi,gb,nf,mf,nx,ny,dx,dy,eta,number_of_grain,x1_grain,x2_grain,r_nuclei):
    for i in range(1,number_of_grain):
        x_nuclei = x1_grain[i-1]
        y_nuclei = x2_grain[i-1]

        for m in range(ny):
            for l in range(nx):
                if l > nx-1:
                    l = l - nx
                if l < 0:
                    l = l + nx
                if m > ny-1:
                    m = m - ny
                if m < 0:
                    m = m + ny
                    
                ddddx = PBC_scalar(l-x_nuclei, nx)
                ddddy = PBC_scalar(m-y_nuclei, ny)
                rr_nuclei = r_nuclei[i-1]    
                r = np.sqrt( (ddddx*dx)**2 +(ddddy*dy)**2 ) - rr_nuclei
                
                tmp = r*np.pi/eta
                phi_tmp = 0.5*(1.-np.sin(tmp))
                if tmp >= np.pi/2.:
                    phi_tmp=0.
                if tmp <= -np.pi/2.:
                    phi_tmp=1.
                    nf[l,m] = nf[l,m]-1
                if phi_tmp > 0:
                    nf_tmp = nf[l,m]+1
                    nf[l,m] = nf_tmp
                    mf[nf_tmp,l,m] = i
                    phi[i,l,m] = phi_tmp
                    phi[0,l,m] = phi[0,l,m]-phi[i,l,m]

def UpdateGB(phi,gb,nx,ny):
    for m in range(0,ny):
        for l in range(0,nx):
            gb[l,m] = np.sum(phi[:,l,m]*phi[:,l,m])

def check_grains(phin,phio,number_of_grain):
    check=False
    for i in range(number_of_grain):
        if(np.abs((phin[i,:,:].max()-phio[i,:,:].max()))>1E-2):
            check=True
    return check
    
#@nb.jit(nopython = False,parallel=True)
def overlaps(phi,number_of_grain):
    nonzero_product_indices2 = []
    nonzero_product_indices3 = []
    for i in range(number_of_grain):
        for j in range(i + 1, number_of_grain):
            product = phi[i, :, :] * phi[j, :, :]
            if np.any(product != 0):
                nonzero_product_indices2.append((i, j))
                for k in range(j + 1, number_of_grain):
                    product3 = phi[i, :, :] * phi[j, :, :] * phi[k, :, :]
                    if np.any(product3 != 0):
                        nonzero_product_indices3.append((i, j, k))
    return nonzero_product_indices2,nonzero_product_indices3

@nb.jit(nopython = False,parallel=True)
def find_gb(phi,nf,mf,nx,ny,dx,dy,number_of_grain, G, sigma11, sigma12, sigma22, sigma11_b, sigma12_b, sigma22_b, sigma12_ext, sigma11_ext, sigma22_ext, criterion, nstep, stress_step,sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,misorientation_list,ref_theta_i,idov2,idov3):
    gb_x_all = []
    gb_y_all = []
    gb_unsort = []
    normal_index = []
    TJ_unsort = []
    index_gb = []
    index_gb_array = []
    index_TJ = []
    
    for i,j in idov2:
        phiphi = phi[i,:,:] * phi[j,:,:]
        gb_ij = np.argwhere(np.abs(phi[i] * phi[j] - 0.25) < criterion)
        if len(gb_ij) > 0.:
            index_gb.append([i, j])
            index_gb_array.append(np.array([i, j]))
            gb_unsort.append(gb_ij)
            normal_index.append(i)

    for i,j,k in idov3:
        phiphiphi = phi[i,:,:] * phi[j,:,:] * phi[k,:,:]
        TJ_ijk = np.argwhere(np.abs(phiphiphi - 0.027) < criterion)
        if len(TJ_ijk) > 1:
            index_TJ.append([i, j, k])
            TJ_unsort.append(TJ_ijk)
                    
    #print("Part I: ",time.time()-t0)
    #t0=time.time()
    flag = False
    if len(TJ_unsort) > 0 and len(gb_unsort) >= 3:
        flag = True
        for pp in range(len(index_TJ)):
            for qq in range(len(index_gb)):
                if set(index_gb[qq]) < set(index_TJ[pp]):
                    if len(TJ_unsort[pp]) > 1:
                        TJ1 = TJ_unsort[pp][0]
                        TJ2 = TJ_unsort[pp][-1]
                        # Numba typing compatibility repair: preserve the pristine
                        # [TJ1, TJ2, original rows...] order with one homogeneous array.
                        old_gb = gb_unsort[qq]
                        joined_gb = np.empty((old_gb.shape[0] + 2, old_gb.shape[1]), dtype=old_gb.dtype)
                        joined_gb[0, :] = TJ1
                        joined_gb[1, :] = TJ2
                        joined_gb[2:, :] = old_gb
                        gb_unsort[qq] = joined_gb
                    elif len(TJ_unsort[pp]) == 1:
                        TJ1 = TJ_unsort[pp][0]
                        # This branch is unreachable for the archived TJ discovery rule
                        # (which retains only arrays with >1 row), but Numba still types it.
                        old_gb = gb_unsort[qq]
                        joined_gb = np.empty((old_gb.shape[0] + 1, old_gb.shape[1]), dtype=old_gb.dtype)
                        joined_gb[0, :] = TJ1
                        joined_gb[1:, :] = old_gb
                        gb_unsort[qq] = joined_gb

    for dd in range(len(gb_unsort)):
        gb_ij = gb_unsort[dd]

        normal_ii = normal_index[dd]
        gb_ij_x, gb_ij_y = sort_gb(gb_ij, flag, normal_ii,nx,ny,dx,dy,phi)

        beta1, beta2 = beta(misorientation_list[index_gb[dd][0]], misorientation_list[index_gb[dd][1]])

        if nstep % stress_step == 0:
            stress_field_bulk(phi,nx,ny,gb_ij_x, gb_ij_y, sigma11_b, sigma12_b, sigma22_b, beta1, beta2,
                              misorientation_list[index_gb[dd][0]], misorientation_list[index_gb[dd][1]],
                              ref_theta_i, flag, nstep, G, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2,
                              sigma12_R2, sigma22_R2,dx,dy)

        gb_x_all.append(gb_ij_x.copy())
        gb_y_all.append(gb_ij_y.copy())

    stress_field_line(phi,nx,ny,gb_x_all, gb_y_all, sigma11, sigma12, sigma22, misorientation_list, index_gb_array,
                      ref_theta_i, flag, nstep, G, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2, sigma12_ext, sigma11_ext, sigma22_ext,dx,dy)

    #print("Part II: ",time.time()-t0)
    #t0=time.time()
    extend_loop = 0
    for extend_loop in range(5):
        stress_field_extend(phi, nx, ny, dx, dy, gb_x_all, gb_y_all, nf, mf, sigma11, sigma12, sigma22)
    #print("Part III: ",time.time()-t0)
    return gb_x_all, gb_y_all
    
    

@nb.njit(fastmath=False)
def sort_gb(gb_ij, flag, normal_ii, nx, ny, dx, dy, phi):
    gb_ij_x = gb_ij[:, 0]
    gb_ij_y = gb_ij[:, 1]

    gb_ij_x_new = np.empty(1, dtype=gb_ij_x.dtype)
    gb_ij_y_new = np.empty(1, dtype=gb_ij_y.dtype)

    gb_ij_x_new[0] = gb_ij_x[0]
    gb_ij_y_new[0] = gb_ij_y[0]

    gb_ij_x = np.delete(gb_ij_x, 0)
    gb_ij_y = np.delete(gb_ij_y, 0)

    lll = 0
    while gb_ij_x.size > 0:
        disxxx = gb_ij_x - gb_ij_x_new[lll]
        Disxxx = np.where(np.abs(disxxx) > nx / 2, disxxx - np.sign(disxxx) * nx, disxxx)
        disyyy = gb_ij_y - gb_ij_y_new[lll]
        Disyyy = np.where(np.abs(disyyy) > ny / 2, disyyy - np.sign(disyyy) * ny, disyyy)
        dis = Disxxx ** 2 + Disyyy ** 2

        # Numba compatibility: np.argmin preserves Python list.index(min(...)) first-tie behavior.
        min_index = np.argmin(dis)
        minimum_distance = dis[min_index]

        if minimum_distance < 1e-4:
            gb_ij_x = np.delete(gb_ij_x, min_index)
            gb_ij_y = np.delete(gb_ij_y, min_index)
        elif minimum_distance >= 1e-4:
            lll += 1
            gb_ij_x_new = np.append(gb_ij_x_new, gb_ij_x[min_index])
            gb_ij_y_new = np.append(gb_ij_y_new, gb_ij_y[min_index])
            gb_ij_x = np.delete(gb_ij_x, min_index)
            gb_ij_y = np.delete(gb_ij_y, min_index)

    #JJU = ((gb_ij_x_new[0]-gb_ij_x_new[len(gb_ij_x_new)//5])**2 + (gb_ij_y_new[0]-gb_ij_y_new[len(gb_ij_x_new)//5])**2) > ((gb_ij_x_new[0]-gb_ij_x_new[-1])**2 + (gb_ij_y_new[0]-gb_ij_y_new[-1])**2)
    gb_ij_x_new = np.delete(gb_ij_x_new, 0)
    #gb_ij_x_new = np.delete(gb_ij_x_new, 0)
    #gb_ij_x_new = np.delete(gb_ij_x_new, -1)
    gb_ij_x_new = np.delete(gb_ij_x_new, -1)
    
    gb_ij_y_new = np.delete(gb_ij_y_new, 0)
    #gb_ij_y_new = np.delete(gb_ij_y_new, 0)
    #gb_ij_y_new = np.delete(gb_ij_y_new, -1)
    gb_ij_y_new = np.delete(gb_ij_y_new, -1)    

    if len(gb_ij_x_new) > 1:
        
        signn = 0
        
        loc = np.zeros(20, dtype=np.int64)
        for qq in range(len(loc)):
            loc[qq] = int(int(qq*(len(gb_ij_x_new)-1)//20)+1)
        #np.random.randint(len(gb_ij_x_new)-1, size=10)+1
        #loc = [len(gb_ij_x_new)//4 + 1]

        signn = clockorcounterclock(gb_ij_x_new, gb_ij_y_new, phi, nx, ny, loc, normal_ii, dx, dy, signn)

        if signn < 0.:
            #print('Flip')
            gb_ij_x_new = np.flip(gb_ij_x_new)
            gb_ij_y_new = np.flip(gb_ij_y_new)
        #else: print('not Flip')
    return gb_ij_x_new, gb_ij_y_new



@nb.njit(fastmath=False)
def clockorcounterclock(gb_ij_x_new, gb_ij_y_new, phi, nx, ny, loc, normal_ii, dx, dy, vvv):
    
    for i_loc in range(len(loc)):
        indecc = loc[i_loc]
        gb_x_dd = PBC_scalar(gb_ij_x_new[indecc] - gb_ij_x_new[indecc - 1], nx)
        gb_y_dd = PBC_scalar(gb_ij_y_new[indecc] - gb_ij_y_new[indecc - 1], ny)
        # Numba compatibility: retain the exact z-component cross-product algebra as scalars.

        gb_ij_l = int(gb_ij_x_new[indecc])
        gb_ij_m = int(gb_ij_y_new[indecc])
        gb_ij_l_m = (gb_ij_l - 1) % nx
        gb_ij_m_m = (gb_ij_m - 1) % ny
    
        norm_xx = (phi[normal_ii, gb_ij_l, gb_ij_m] - phi[normal_ii, gb_ij_l_m, gb_ij_m]) / dx
        norm_yy = (phi[normal_ii, gb_ij_l, gb_ij_m] - phi[normal_ii, gb_ij_l, gb_ij_m_m]) / dy
        LL_N = np.sqrt(norm_xx ** 2 + norm_yy ** 2)
        cross_z = -gb_x_dd * norm_yy + gb_y_dd * norm_xx
        vvv += np.sign(cross_z)
    
    return vvv




@nb.njit(nopython=True,parallel=True)
def stress_field_line(phi, nx,ny, gb_x, gb_y, sigma11, sigma12, sigma22, misorientation_list, index_gb, ref_theta_i, flag, nstep, G, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,sigma12_ext,sigma11_ext,sigma22_ext,dx,dy):
    coeff = np.sqrt(2)/2
    for gb_id in range(len(gb_x)): 
#        print(index_gb[gb_id][0], index_gb[gb_id][1])

        for oo in range(0, len(gb_x[gb_id])):
            mm = int(gb_x[gb_id][oo])        
            nn = int(gb_y[gb_id][oo])    

            II_12 = 0
            II_11 = 0
            II_22 = 0

    
            for num_k in range(0, len(gb_x)): 
                
                
                OR1 = misorientation_list[index_gb[num_k][0]]
                OR2 = misorientation_list[index_gb[num_k][1]]
                MOR = (OR1 + OR2)/2
                beta1, beta2 = beta(misorientation_list[index_gb[num_k][0]], misorientation_list[index_gb[num_k][1]])
                I1_12 = 0
                I2_12 = 0
                I3_12 = 0
                I4_12 = 0        
        
                I1_11 = 0
                I2_11 = 0
                I3_11 = 0
                I4_11 = 0 
                
                I1_22 = 0
                I2_22 = 0
                I3_22 = 0
                I4_22 = 0 
                            
                for kk in range(0, len(gb_x[num_k])):
                    
                    Ddis_x0 = mm - gb_x[num_k][kk]
                    Ddis_y0 = nn - gb_y[num_k][kk]
                    if Ddis_x0 < 0:      Ddis_x0 = Ddis_x0 + nx           
                    if Ddis_y0 < 0:      Ddis_y0 = Ddis_y0 + ny          
                    
                    tao12_1ref = sigma12_R1[Ddis_x0, Ddis_y0]
                    tao11_1ref = sigma11_R1[Ddis_x0, Ddis_y0]
                    tao22_1ref = sigma22_R1[Ddis_x0, Ddis_y0]
                    tao12_3ref = sigma12_R2[Ddis_x0, Ddis_y0]
                    tao11_3ref = sigma11_R2[Ddis_x0, Ddis_y0]
                    tao22_3ref = sigma22_R2[Ddis_x0, Ddis_y0]

                    
                    tao12_1 = tao12_1ref*np.cos(MOR) + tao12_3ref*np.sin(MOR)
                    tao11_1 = tao11_1ref*np.cos(MOR) + tao11_3ref*np.sin(MOR)
                    tao22_1 = tao22_1ref*np.cos(MOR) + tao22_3ref*np.sin(MOR)
        
                    tao12_3 = tao12_1ref*np.cos(MOR+np.pi/2) + tao12_3ref*np.sin(MOR+np.pi/2)
                    tao11_3 = tao11_1ref*np.cos(MOR+np.pi/2) + tao11_3ref*np.sin(MOR+np.pi/2)
                    tao22_3 = tao22_1ref*np.cos(MOR+np.pi/2) + tao22_3ref*np.sin(MOR+np.pi/2)
                   
                    tao11_2 = coeff*(tao11_1+tao11_3)
                    tao11_4 = coeff*(-tao11_1+tao11_3)
                    tao12_2 = coeff*(tao12_1+tao12_3)
                    tao12_4 = coeff*(-tao12_1+tao12_3)
                    tao22_2 = coeff*(tao22_1+tao22_3)
                    tao22_4 = coeff*(-tao22_1+tao22_3)
        
                    if kk == 0:
                        dx_dss = gb_x[num_k][1] - gb_x[num_k][0]
                        dy_dss = gb_y[num_k][1] - gb_y[num_k][0]
                        dx_ds = -PBC_scalar(dx_dss, nx) *dx
                        dy_ds = -PBC_scalar(dy_dss, ny) *dy
        
                    if kk == len(gb_x[num_k])-1:
                        dx_dss = gb_x[num_k][kk] - gb_x[num_k][kk-1]
                        dy_dss = gb_y[num_k][kk] - gb_y[num_k][kk-1]
                        dx_ds = -PBC_scalar(dx_dss, nx) *dx
                        dy_ds = -PBC_scalar(dy_dss, ny) *dy
        
                    if kk != len(gb_x[num_k])-1 and kk != 0:
                        dx_dss = gb_x[num_k][kk+1] - gb_x[num_k][kk-1]
                        dy_dss = gb_y[num_k][kk+1] - gb_y[num_k][kk-1]
                        dx_ds = -PBC_scalar(dx_dss, nx)/2 *dx
                        dy_ds = -PBC_scalar(dy_dss, ny)/2 *dy
                        
                        
                    #ddx_dds = np.cos(MOR) * dx_ds + np.sin(MOR) * dy_ds
                    #ddy_dds =-np.sin(MOR) * dx_ds + np.cos(MOR) * dy_ds

                    ddy_dds = -(np.cos(MOR) * dx_ds + np.sin(MOR) * dy_ds)
                    ddx_dds =-np.sin(MOR) * dx_ds + np.cos(MOR) * dy_ds

                    
                    if ddx_dds != 0. and ddy_dds/ddx_dds >= 0. and ddy_dds/ddx_dds < 1.: 
    
                        nnni = 0.
                        phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                        phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                        
                        dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                        dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                        
                        I1_12 += dY_dS  * tao12_1
                        I2_12 += -dX_dS * tao12_2
                        I1_11 += dY_dS  * tao11_1
                        I2_11 += -dX_dS * tao11_2
                        I1_22 += dY_dS  * tao22_1
                        I2_22 += -dX_dS * tao22_2
    #                    LLL = np.array([1/(np.sin(phi_k1-phi_k)) * beta1, 1/(np.sin(phi_k1-phi_k)) * beta2])
    
                        
                    if ddx_dds != 0. and ddy_dds/ddx_dds >= 1.:
                        nnni = 1.0
                        phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                        phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                        
                        dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                        dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                        
                        I2_12 += dY_dS  * tao12_2
                        I3_12 += -dX_dS * tao12_3
                        I2_11 += dY_dS  * tao11_2
                        I3_11 += -dX_dS * tao11_3
                        I2_22 += dY_dS  * tao22_2
                        I3_22 += -dX_dS * tao22_3
    #                    LLL = np.array([1/(np.sin(phi_k1-phi_k)) * beta2, 1/(np.sin(phi_k1-phi_k)) * beta3])
    
                        
                    if ddx_dds == 0. or (ddx_dds != 0. and ddy_dds/ddx_dds < -1.):              
                        nnni = 2.0
                        phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                        phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                        
                        dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                        dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                        
                        I3_12 += dY_dS  * tao12_3
                        I4_12 += -dX_dS * tao12_4
                        I3_11 += dY_dS  * tao11_3
                        I4_11 += -dX_dS * tao11_4
                        I3_22 += dY_dS  * tao22_3
                        I4_22 += -dX_dS * tao22_4
    #                    LLL = np.array([1/(np.sin(phi_k1-phi_k)) * beta2, 1/(np.sin(phi_k1-phi_k)) * beta3])
    
    
                    if ddx_dds != 0. and ddy_dds/ddx_dds >= -1. and ddy_dds/ddx_dds < 0.:                
                        nnni = 3.0
                        phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                        phi_k1 = 0.0*ref_theta_i              + (OR1 + OR2)/2
                        
                        dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                        dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                        
                        I4_12 += -dY_dS  * tao12_4
                        I1_12 += dX_dS * tao12_1
                        I4_11 += -dY_dS  * tao11_4
                        I1_11 += dX_dS * tao11_1
                        I4_22 += -dY_dS  * tao22_4
                        I1_22 += dX_dS * tao22_1
    #                    LLL = np.array([1/(np.sin(phi_k1-phi_k)) * beta2, 1/(np.sin(phi_k1-phi_k)) * beta3])
     
                beta3 = beta1
                beta4 = beta2      
                
                
                II_12 += beta1 * I1_12 + beta2 * I2_12 + beta3 * I3_12 + beta4 * I4_12
                II_11 += beta1 * I1_11 + beta2 * I2_11 + beta3 * I3_11 + beta4 * I4_11
                II_22 += beta1 * I1_22 + beta2 * I2_22 + beta3 * I3_22 + beta4 * I4_22
    
                
            sigma11_v = G/coeff * II_11 + sigma11_ext
            sigma12_v = G/coeff * II_12 + sigma12_ext
            sigma22_v = G/coeff * II_22 + sigma22_ext
                
            sigma11[mm, nn] += sigma11_v
            sigma12[mm, nn] += sigma12_v
            sigma22[mm, nn] += sigma22_v
        

@nb.njit(parallel=True)
def judge(x, y, gb_x, gb_y, Flag):
    Flag = True
    for pp in range(0, len(gb_x)):
        for qq in range(0, len(gb_x[pp])):
            if x == int(gb_x[pp][qq]) and y == int(gb_y[pp][qq]): Flag = False
    return Flag


@nb.njit(nopython=True, parallel=True)
def stress_field_extend(phi, nx, ny, dx, dy, gb_ij_x, gb_ij_y, nf, mf, sigma11, sigma12, sigma22):
    for m in nb.prange(ny):
    #for m in range(ny):
        for l in range(nx):
            l_p = (l + 1) % nx
            l_m = (l - 1 + nx) % nx
            m_p = (m + 1) % ny
            m_m = (m - 1 + ny) % ny
            Flaagg = -10
            Flaagg = judge(l, m, gb_ij_x, gb_ij_y, Flaagg)
            
            if Flaagg:
                for n1 in range(nf[l, m]):
                    i = mf[n1, l, m]
                    q, qq, normalx, normaly = grad(phi, i, l, m, l_p, l_m, m_p, m_m, dx, dy)
                    I = -int(np.sign(np.sign(normalx) * np.sign(phi[i, l, m] - 0.5)))
                    J = -int(np.sign(np.sign(normaly) * np.sign(phi[i, l, m] - 0.5)))
                    
                    llnew = (l + I) % nx
                    mmnew = (m + J) % ny
                    sigma11[l, m] += 1e-0 * (np.abs(normalx) * (sigma11[llnew, m] - sigma11[l, m]) +
                                              np.abs(normaly) * (sigma11[l, mmnew] - sigma11[l, m]))
                    sigma12[l, m] += 1e-0 * (np.abs(normalx) * (sigma12[llnew, m] - sigma12[l, m]) +
                                              np.abs(normaly) * (sigma12[l, mmnew] - sigma12[l, m]))
                    sigma22[l, m] += 1e-0 * (np.abs(normalx) * (sigma22[llnew, m] - sigma22[l, m]) +
                                              np.abs(normaly) * (sigma22[l, mmnew] - sigma22[l, m]))
    return sigma11, sigma12, sigma22
    #print("Stress Extend: ",t0-time.time())


@nb.njit(parallel=True)
def stress_field_bulk(phi, nx,ny, gb_x, gb_y, sigma11_b, sigma12_b, sigma22_b, beta1, beta2, OR1, OR2, ref_theta_i, flag, nstep, G, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,dx,dy):
    coeff = np.sqrt(2)/2
    for m in range(ny):
        for l in range(nx):
            I1_12 = 0
            I2_12 = 0
            I3_12 = 0
            I4_12 = 0        
    
            I1_11 = 0
            I2_11 = 0
            I3_11 = 0
            I4_11 = 0 
            
            I1_22 = 0
            I2_22 = 0
            I3_22 = 0
            I4_22 = 0 

            for kk in range(0, len(gb_x)):
                
                MOR = (OR1 + OR2)/2
                
                
                Ddis_x0 = l - gb_x[kk]
                Ddis_y0 = m - gb_y[kk]

                if Ddis_x0 < 0:      Ddis_x0 = Ddis_x0 + nx
                if Ddis_y0 < 0:      Ddis_y0 = Ddis_y0 + ny
                
                
                tao12_1ref = sigma12_R1[Ddis_x0, Ddis_y0]
                tao11_1ref = sigma11_R1[Ddis_x0, Ddis_y0]
                tao22_1ref = sigma22_R1[Ddis_x0, Ddis_y0]
                tao12_3ref = sigma12_R2[Ddis_x0, Ddis_y0]
                tao11_3ref = sigma11_R2[Ddis_x0, Ddis_y0]
                tao22_3ref = sigma22_R2[Ddis_x0, Ddis_y0]
                
                tao12_1 = tao12_1ref*np.cos(MOR) + tao12_3ref*np.sin(MOR)
                tao11_1 = tao11_1ref*np.cos(MOR) + tao11_3ref*np.sin(MOR)
                tao22_1 = tao22_1ref*np.cos(MOR) + tao22_3ref*np.sin(MOR)
    
                tao12_3 = tao12_1ref*np.cos(MOR+np.pi/2) + tao12_3ref*np.sin(MOR+np.pi/2)
                tao11_3 = tao11_1ref*np.cos(MOR+np.pi/2) + tao11_3ref*np.sin(MOR+np.pi/2)
                tao22_3 = tao22_1ref*np.cos(MOR+np.pi/2) + tao22_3ref*np.sin(MOR+np.pi/2)

                tao11_2 = coeff*(tao11_1+tao11_3)
                tao11_4 = coeff*(-tao11_1+tao11_3)
                tao12_2 = coeff*(tao12_1+tao12_3)
                tao12_4 = coeff*(-tao12_1+tao12_3)
                tao22_2 = coeff*(tao22_1+tao22_3)
                tao22_4 = coeff*(-tao22_1+tao22_3)

                if kk == 0:
                    dx_dss = gb_x[1] - gb_x[0]
                    dy_dss = gb_y[1] - gb_y[0]
                    dx_ds = -PBC_scalar(dx_dss, nx) *dx
                    dy_ds = -PBC_scalar(dy_dss, ny) *dy
    
                if kk == len(gb_x)-1:
                    dx_dss = gb_x[kk] - gb_x[kk-1]
                    dy_dss = gb_y[kk] - gb_y[kk-1]
                    dx_ds = -PBC_scalar(dx_dss, nx) *dx
                    dy_ds = -PBC_scalar(dy_dss, ny) *dy
    
                if kk != len(gb_x)-1 and kk != 0:
                    dx_dss = gb_x[kk+1] - gb_x[kk-1]
                    dy_dss = gb_y[kk+1] - gb_y[kk-1]
                    dx_ds = -PBC_scalar(dx_dss, nx)/2 *dx
                    dy_ds = -PBC_scalar(dy_dss, ny)/2 *dy
                   
        
                ddy_dds = -(np.cos(MOR) * dx_ds + np.sin(MOR) * dy_ds)
                ddx_dds =-np.sin(MOR) * dx_ds + np.cos(MOR) * dy_ds

                if ddx_dds != 0. and ddy_dds/ddx_dds >= 0. and ddy_dds/ddx_dds < 1.: 

                    nnni = 0.
                    phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                    phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                    
                    dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                    dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                    
                    I1_12 += dY_dS  * tao12_1
                    I2_12 += -dX_dS * tao12_2
                    I1_11 += dY_dS  * tao11_1
                    I2_11 += -dX_dS * tao11_2
                    I1_22 += dY_dS  * tao22_1
                    I2_22 += -dX_dS * tao22_2
                    beta10 = beta1
                    beta20 = beta2

                    
                if ddx_dds != 0. and ddy_dds/ddx_dds >= 1.:
                    nnni = 1.0
                    phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                    phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                    
                    dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                    dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                    
                    I2_12 += dY_dS  * tao12_2
                    I3_12 += -dX_dS * tao12_3
                    I2_11 += dY_dS  * tao11_2
                    I3_11 += -dX_dS * tao11_3
                    I2_22 += dY_dS  * tao22_2
                    I3_22 += -dX_dS * tao22_3
                    beta10 = beta1
                    beta20 = beta2

                    
                if ddx_dds == 0. or (ddx_dds != 0. and ddy_dds/ddx_dds < -1.):              
                    nnni = 2.0
                    phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                    phi_k1 = (nnni+1.)*ref_theta_i        + (OR1 + OR2)/2
                    
                    dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                    dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                   
                    I3_12 += dY_dS  * tao12_3
                    I4_12 += -dX_dS * tao12_4
                    I3_11 += dY_dS  * tao11_3
                    I4_11 += -dX_dS * tao11_4
                    I3_22 += dY_dS  * tao22_3
                    I4_22 += -dX_dS * tao22_4
                    beta10 = beta1
                    beta20 = beta2

                if ddx_dds != 0. and ddy_dds/ddx_dds >= -1. and ddy_dds/ddx_dds < 0.:                
                    nnni = 3.0
                    phi_k = nnni*ref_theta_i              + (OR1 + OR2)/2
                    phi_k1 = 0.0*ref_theta_i              + (OR1 + OR2)/2
                    
                    dX_dS = 1/(np.sin(phi_k1 - phi_k)) * (np.sin(phi_k1) * dx_ds - np.cos(phi_k1) * dy_ds)
                    dY_dS = 1/(np.sin(phi_k1 - phi_k)) * (-np.sin(phi_k) * dx_ds + np.cos(phi_k) * dy_ds)
                    
                    I4_12 += -dY_dS  * tao12_4
                    I1_12 += dX_dS * tao12_1
                    I4_11 += -dY_dS  * tao11_4
                    I1_11 += dX_dS * tao11_1
                    I4_22 += -dY_dS  * tao22_4
                    I1_22 += dX_dS * tao22_1
                    beta10 = beta1
                    beta20 = beta2

            beta30 = beta10
            beta40 = beta20

            sigma11_v_b = G/coeff * (beta10 * I1_11 + beta20 * I2_11 + beta30 * I3_11 + beta40 * I4_11)
            sigma12_v_b = G/coeff * (beta10 * I1_12 + beta20 * I2_12 + beta30 * I3_12 + beta40 * I4_12) #+ tao_ext
            sigma22_v_b = G/coeff * (beta10 * I1_22 + beta20 * I2_22 + beta30 * I3_22 + beta40 * I4_22) #+ tao_ext
                
            sigma11_b[l, m] += sigma11_v_b
            sigma12_b[l, m] += sigma12_v_b
            sigma22_b[l, m] += sigma22_v_b




@nb.njit(parallel=True)
def stress_field_bulk_single(phi, nx, ny, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,dx,dy):
    for m in nb.prange(ny):
        for l in range(nx):

            Ddis_x0 = l
            Ddis_y0 = m
            DDis_x0 = PBC_scalar(Ddis_x0, nx)
            DDis_y0 = PBC_scalar(Ddis_y0, ny)
            
            tao12_1ref = 0
            tao11_1ref = 0
            tao22_1ref = 0
            tao12_2ref = 0
            tao11_2ref = 0
            tao22_2ref = 0
            ppxy = np.array([[0,0], [nx,0], [-nx,0], [0,nx], [0,-nx]])
            
            for pppp in range(0, len(ppxy)):
                DDis_x = (DDis_x0 + ppxy[pppp][0])*dx
                DDis_y = (DDis_y0 + ppxy[pppp][1])*dy
                rho2 = (DDis_x**2 + DDis_y**2 + 0.01**2)
    
                tao12_1ref += (DDis_x/rho2)*(1-2*DDis_y**2/rho2)
                tao11_1ref += -(DDis_y/rho2)*(1+2*DDis_x**2/rho2)
                tao22_1ref += (DDis_y/rho2)*(1-2*DDis_y**2/rho2)
    
                tao12_2ref += (DDis_y/rho2)*(1-2*DDis_y**2/rho2)
                tao11_2ref += (DDis_x/rho2)*(1-2*DDis_y**2/rho2)
                tao22_2ref += (DDis_x/rho2)*(1+2*DDis_y**2/rho2)
                                        
                
            sigma11_R1[l, m] = tao11_1ref
            sigma12_R1[l, m] = tao12_1ref
            sigma22_R1[l, m] = tao22_1ref
            sigma11_R2[l, m] = tao11_2ref
            sigma12_R2[l, m] = tao12_2ref
            sigma22_R2[l, m] = tao22_2ref

    return sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2


#@nb.jit(parallel=True)
def renorm(phi,phi_new, number_of_grain):
    phi_new[phi_new<0.0]=0.0
    phi_new[phi_new>1.0]=1.0
    a=phi_new[0,:,:]
    for j in nb.prange(0,number_of_grain):
        a=phi_new[0,:,:]
        for i in range(1,number_of_grain):
            a=a+phi_new[i,:,:]
        phi[j]=phi_new[j,:,:]/a



@nb.njit(parallel=True)
def PBC_scalar(x, NN):
    NNN = NN
    if x >= NNN / 2:
        return x - NN
    if x <= -NNN / 2:
        return x + NN
    # The three finite scalar branches are exhaustive for native integer offsets.
    return x
    
    
    
@nb.njit(parallel=True)
def PBC_list(x, NN):
    x = np.array(x)
    return list(np.where(abs(x)>NN/2, x - np.sign(x)*NN, x))




@nb.njit(nopython = True, fastmath=False, parallel=True)
def update_nfmf(phi,mf,nf,nx,ny,number_of_grain):
    #for m in range(ny):
    for m in nb.prange(ny):
        m_p = (m + 1) % ny
        m_m = (m - 1) % ny
        for l in range(nx):
            l_p = (l + 1) % nx
            l_m = (l - 1) % nx
            n = 0
            for i in range(number_of_grain):
                if phi[i,l,m] > 0.0 or (phi[i,l,m] == 0.0 and (phi[i,l_p,m] > 0.0 or phi[i,l_m,m] > 0.0 or phi[i,l,m_p] > 0.0 or phi[i,l,m_m] > 0.0)):
                    n += 1
                    mf[n-1,l,m] = i
            nf[l,m] = n
            
            
@nb.njit(parallel=True)    
def beta(ORR1, ORR2):
    TTHETA = ORR1 - ORR2
    if TTHETA >= np.pi/2 and TTHETA < np.pi:
        TTHETA = TTHETA - np.pi/2
    if TTHETA >= np.pi and TTHETA < 3*np.pi/2:
        TTHETA = TTHETA - 2*np.pi/2
    if TTHETA >= 3*np.pi/2 and TTHETA < 2*np.pi:
        TTHETA = TTHETA - 3*np.pi/2
        
    if TTHETA <= -np.pi/2 and TTHETA > -np.pi:
        TTHETA = TTHETA + np.pi/2    
    if TTHETA <= -np.pi and TTHETA > -3*np.pi/2:
        TTHETA = TTHETA + 2*np.pi/2    
    if TTHETA <= -3*np.pi/2  and TTHETA > -2*np.pi:
        TTHETA = TTHETA + 3*np.pi/2    
        
        
    if abs(TTHETA) <= 36.7*np.pi/180: 
        beta1 = 2.*np.tan(TTHETA/2.)    
        beta2 = - beta1
        
    if abs(TTHETA) > 36.7*np.pi/180 and abs(TTHETA) <= 53.7*np.pi/180:
        beta1 = 2.*np.tan(TTHETA/2. - np.sign(TTHETA)/2*53.7*np.pi/180)   *0.0
        beta2 = -beta1
    
    if abs(TTHETA) > 53.7*np.pi/180: 
        beta1 = 2.*np.tan(TTHETA/2. - np.sign(TTHETA)*np.pi/4)   *0.0   
        beta2 = - beta1
    return beta1, beta2



@nb.njit(nopython=False,parallel=True)
def update_PF(phi, phi_new, nx, ny, dx, dy, g1, g2, eta, ref_theta_i, sigma11, sigma12, sigma22, pmobi, A, mf, nf, eij, misorientation_list, dt):
    for m in range(ny):
        for l in range(nx):
            l_p = (l + 1) % nx
            l_m = (l - 1) % nx
            m_p = (m + 1) % ny
            m_m = (m - 1) % ny
            for n1 in range(nf[l, m]):
                i = mf[n1, l, m]
                dpi = 0.0
                # Stencil calculation for phi[i]
                stencil_i = (phi[i, l_m, m_m] + 4 * phi[i, l_m, m] + 4 * phi[i, l_p, m] +
                             phi[i, l_p, m_p] + phi[i, l_m, m_p] + 4 * phi[i, l, m_m] +
                             4 * phi[i, l, m_p] + phi[i, l_p, m_m] - 20 * phi[i, l, m]) / (6 * dx * dx)

                for n2 in range(nf[l, m]):
                    j = mf[n2, l, m]

                    # Stencil calculation for phi[j]
                    stencil_j = (phi[j, l_m, m_m] + 4 * phi[j, l_m, m] + 4 * phi[j, l_p, m] +
                                 phi[j, l_p, m_p] + phi[j, l_m, m_p] + 4 * phi[j, l, m_m] +
                                 4 * phi[j, l, m_p] + phi[j, l_p, m_m] - 20 * phi[j, l, m]) / (6 * dx * dx)
                    
                    # Capillarity force
                    ppp = ((phi[j, l, m] * stencil_i - phi[i, l, m] * stencil_j) + np.pi**2 / (2 * eta**2) * (phi[i, l, m] - phi[j, l, m]))

                    # Elastic energy
                    E_el = -E_elastic(phi, i, j, l, m, l_p, l_m, m_p, m_m, ref_theta_i, misorientation_list[i], misorientation_list[j], sigma11[l, m], sigma12[l, m], sigma22[l, m], dx, dy)

                    # Mobility calculation
                    if i != j: mij = pmobi
                    else: mij = 0
                    phii_phij = phi[i, l, m] * phi[j, l, m]

                    dpi += mij * (ppp - E_el + np.pi / eta * np.sqrt(phii_phij) * eij[i, j])
                
                phi_new[i, l, m] = phi[i, l, m] + dpi * dt
                
                #if phi_new[i, l, m] >= 0.5:
                #    A[i] += 1
                    
    return phi_new#, #A
    
    
@nb.njit(parallel=True)
def E_elastic(phi, i, j, l, m, l_p, l_m, m_p, m_m, ref_theta_i, OR1, OR2, sigma11_v, sigma12_v, sigma22_v,dx,dy):
    nabla_i_x, nabla_i_y, nabla_i_xnn, nabla_i_ynn = grad(phi, i, l, m, l_p, l_m, m_p, m_m,dx,dy)
    nabla_j_x, nabla_j_y, nabla_j_xnn, nabla_j_ynn = grad(phi, j, l, m, l_p, l_m, m_p, m_m,dx,dy)
    
    normal_j_x = phi[i, l, m]*nabla_j_x - phi[j, l, m]*nabla_i_x
    normal_j_y = phi[i, l, m]*nabla_j_y - phi[j, l, m]*nabla_i_y

    MOR = (OR1 + OR2)/2

    dnormal_j_x = np.cos(MOR) * normal_j_x + np.sin(MOR) * normal_j_y
    dnormal_j_y =-np.sin(MOR) * normal_j_x + np.cos(MOR) * normal_j_y

    beta1, beta2 = beta(OR1, OR2)
    beta3 = beta1
    beta4 = beta2  
    
    
    if dnormal_j_y != 0. and -dnormal_j_x/dnormal_j_y >= 0. and -dnormal_j_x/dnormal_j_y < 1.:  
        L1 = beta1
        L2 = beta2
        phi_k = 0.0*ref_theta_i              + (OR1 + OR2)/2
        phi_k1 = 1.0*ref_theta_i             + (OR1 + OR2)/2
        
    if dnormal_j_y != 0. and -dnormal_j_x/dnormal_j_y >= 1.:  
        L1 = beta2
        L2 = beta3
        phi_k = 1.0*ref_theta_i              + (OR1 + OR2)/2
        phi_k1 = 2.0*ref_theta_i             + (OR1 + OR2)/2
        
    if dnormal_j_y == 0. or (-dnormal_j_y != 0. and -dnormal_j_x/dnormal_j_y < -1.):  
        L1 = beta3
        L2 = beta4
        phi_k = 2.0*ref_theta_i              + (OR1 + OR2)/2
        phi_k1 = 3.0*ref_theta_i             + (OR1 + OR2)/2
        
        
    if dnormal_j_y != 0. and -dnormal_j_x/dnormal_j_y >= -1. and -dnormal_j_x/dnormal_j_y < 0.:  
        L1 = beta4
        L2 = beta1
        phi_k = 3.0*ref_theta_i              + (OR1 + OR2)/2
        phi_k1 = 0.0*ref_theta_i             + (OR1 + OR2)/2
    
    norm_pro = nabla_i_xnn * nabla_j_xnn +  nabla_i_ynn * nabla_j_ynn
    Lambda1 =  L1
    Lambda2 =  L2


    tao_k = (sigma22_v-sigma11_v)/2. * np.sin(2.*phi_k) + sigma12_v * np.cos(2.*phi_k)
    tao_k1 = (sigma22_v-sigma11_v)/2. * np.sin(2.*phi_k1) + sigma12_v * np.cos(2.*phi_k1)    

    
    E_el  =  (- tao_k * Lambda1 - tao_k1 * Lambda2) * np.hypot(nabla_j_x, nabla_j_y) * norm_pro

    return E_el

            

def plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,nstep):
    UpdateGB(phi,gb,nx,ny)
    fig, ax = plt.subplots(figsize=(3.0, 3), dpi = 300)
    plt.xticks([])
    plt.yticks([])
    plt.imshow(gb, cmap='gray')
    plt.title('t = ' + str(nstep))
    namefile = savename+'/t-'+str(nstep)+'.jpg'
    if output_flag:
        plt.savefig(namefile, bbox_inches = 'tight', dpi = 600, transparent = True, edgecolor='w')
        plt.close()
        #np.save(savename+ '/OP_t'+str(nstep)+'.npy', phi)
        sparse_phi = sparse.COO(phi)
        sparse.save_npz(savename+ '/OP_t'+str(nstep)+'.npz', sparse_phi)

    if not output_flag: plt.show()

def plotstress(gb_x_t,gb_y_t,sigma12_b, sigma11_b, sigma22_b,savename,output_flag,nstep):
    fig, ax = plt.subplots(figsize=(3.0, 3), dpi = 300)
    plt.xticks([])
    plt.yticks([])
    for nuummm in range(len(gb_x_t)):
        plt.scatter(gb_y_t[len(gb_x_t)-1-nuummm], gb_x_t[len(gb_x_t)-1-nuummm], s=0.5, marker = 's', c='black')
    plt.imshow(sigma12_b, cmap='coolwarm', vmin=-0.3e-1, vmax=0.3e-1)
    plt.title(r'$\sigma_{12}$, t = ' + str(nstep))
    namefile = savename+'/stress-t-'+str(nstep)+'.jpg'
    if output_flag:
        plt.savefig(namefile, bbox_inches = 'tight', dpi = 600, transparent = True, edgecolor='w')
        plt.close()
        np.savez(savename+ '/Stress_t'+str(nstep)+'.npz', sigma11_b, sigma22_b, sigma12_b)
    if not output_flag:         plt.show()


@nb.njit(nopython=False, fastmath=False, parallel=True)
def grad(phi, c, l, m, l_p, l_m, m_p, m_m,dx,dy):
    nabla_c_x_l = (phi[c, l_p, m] - phi[c, l_m, m])/(2*dx)
    nabla_c_y_l = (phi[c, l, m_p] - phi[c, l, m_m])/(2*dy)
    nabla_c_x = nabla_c_x_l / np.sqrt(nabla_c_x_l**2 + nabla_c_y_l**2 + 1E-10)
    nabla_c_y = nabla_c_y_l / np.sqrt(nabla_c_x_l**2 + nabla_c_y_l**2 + 1E-10)
    return nabla_c_x_l , nabla_c_y_l, nabla_c_x, nabla_c_y
            
#Only Idealized_4ref
@nb.njit(nopython = True, fastmath=False, parallel=True)
def cal(sin_phi_i, cos_phi_i, ETA):
    sin_phi_i_over_ETA = 5. * sin_phi_i / ETA
    cos_phi_i_over_ETA = 5. * cos_phi_i / ETA
    sinhsin = np.sinh(sin_phi_i_over_ETA)
    sinhcos = np.sinh(cos_phi_i_over_ETA)
    coshsin = np.cosh(sin_phi_i_over_ETA)
    coshcos = np.cosh(cos_phi_i_over_ETA)
    Rn_2_sin = -sin_phi_i * sinhsin / (1. + coshsin) + \
               (cos_phi_i**2 / ETA) * (coshsin / (coshsin + 1.) - sinhsin**2 / (coshsin + 1.)**2)
    Rn_2_cos = -cos_phi_i * sinhcos / (coshcos + 1.) + \
               (sin_phi_i**2 / ETA) * (coshcos / (coshcos + 1.) - sinhcos**2 / (coshcos + 1.)**2)
    return cos_phi_i, sin_phi_i, Rn_2_cos, Rn_2_sin


#Only Idealized_4ref
@nb.njit(nopython=True, fastmath=False, parallel=True)
def cal_inc(a, b, misorientation_list, l, m, phi_i, ref_theta_i):
    nn_i = phi_i // ref_theta_i
    phi_i_k = nn_i*ref_theta_i              + (misorientation_list[a] + misorientation_list[b])/2
    phi_i_k1 = (nn_i+1)*ref_theta_i         + (misorientation_list[a] + misorientation_list[b])/2
    cos_phi_i = np.sin(phi_i_k1 - phi_i)/np.sin(phi_i_k1 - phi_i_k)
    sin_phi_i = np.sin(phi_i - phi_i_k)/np.sin(phi_i_k1 - phi_i_k)
    return cos_phi_i, sin_phi_i


#Only Idealized_4ref
@nb.njit(nopython=True, fastmath=False, parallel=True)
def inclination(a, b, misorientation_list, l, m, l_p, l_m, m_p, m_m,dx,dy):
    nabla_a_x, nabla_a_y, xx, yy = grad(phi, a, l, m, l_p, l_m, m_p, m_m,dx,dy)
    nabla_b_x, nabla_b_y, xxx, yyy = grad(phi, b, l, m, l_p, l_m, m_p, m_m,dx,dy)
    normal_a_x = phi[b, l, m]*nabla_a_x - phi[a, l, m]*nabla_b_x
    normal_a_y = phi[b, l, m]*nabla_a_y - phi[a, l, m]*nabla_b_y
    phi0_i = np.arctan(normal_a_y / (normal_a_x + 1E-8)) #- (misorientation_list[a] + misorientation_list[b])/2
    
    if phi0_i <0 : phi_i = phi0_i + np.pi
    else: phi_i = phi0_i
    
    cos, sin = cal_inc(a, b, misorientation_list, l, m, phi_i, ref_theta_i)

    return cal(sin, cos, ETA)


@nb.njit(nopython=True, parallel=True)
def update_area(A, phi):
    for i in range(len(phi)):
        AAr = np.sum(phi[i])
        if AAr > 0.0:
            A[i] = np.round(AAr,2)
        else:  A[i] = 0    
    return A