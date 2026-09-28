
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
from functions_4ref_new import *
import sparse


nb.set_num_threads(8)

#########
#########  Parameters and Definitions  ##########
#########
Dim = 16                           #Available 20, 50, 100, 400
nx = 500                            # nx/100 is the size of the sistem for nx>100
ny = nx
number_of_grain = Dim + 1#**2 + 1   # Number of grains
dx, dy = 1, 1                       # Mesh discretization
dt = 0.1                            # Time Step
nsteps = 200000              # Max No. Steps
eta = 5.0 * dx                      # Interface Thickenss
r_nuclei = dx*2.5 #*100/150         # radius of the initial grains
ddd = 0.020                          # criterion for GB detection

# par. PF
eee = 20.0#5.0e+7                   #
pmobi = 2.0*(np.pi**2)/(8.*eta)     #

output_flag = 1                # Enable output
n_grain_area = 50                   # n step grain_area
gb_step = 250                        # n step between gb plots / update
stress_step = 250                   # n step between stress plots


G = 0.02                             # scaling stress coeff.
#beta1 = 0.5
#beta2 = -0.5

sigma12_ext = 0.                         # external shear
sigma11_ext = 0.                        # external shear
sigma22_ext = 0.                         # external shear

#LLL = beta2 - beta1

# Par. Microstructure
# Par. Microstructure

OORR1 = 0.0
OORR2 = -22.6
OORR3 = 28.1


misorientation_list = np.array([0., OORR1, OORR1, OORR1, OORR1, OORR1, OORR1, OORR1, OORR1, 
                                OORR2, OORR2, OORR2, OORR2, OORR3, OORR3, OORR3, OORR3])     *np.pi/180
misorientation_list[0] = 0.

r_nuclei = np.array([12*dx, 12*dx, 12*dx, 12*dx, 12*dx, 12*dx, 12*dx, 12*dx, 
            25*dx, 25*dx, 25*dx, 25*dx, 25*dx, 25*dx, 25*dx, 25*dx])*5./1.5 

x1_grain = [1*nx/2, 1*nx/2, 0.,     0., 1*nx/4, 1*nx/4, 3*nx/4, 3*nx/4,        1*nx/2, 1*nx/2, 0.,     0.,          1*nx/4, 3*nx/4, 1*nx/4, 3*nx/4]
x2_grain = [1*nx/2, 0.,     1*nx/2, 0., 1*nx/4, 3*nx/4, 1*nx/4, 3*nx/4,        1*nx/4, 3*nx/4, 1*nx/4, 3*nx/4,      1*nx/2, 1*nx/2, 0.,     0.]





#Par Disconnection
ETA = 0.8                           # Reg. parameter
ref = 4.                            # number of references
ref_theta_i = np.pi/ref             # angle ref
g1 = 1.                             # int-energy ref1
g2 = 1.                             # int-energy ref1

#Create output folder
savename = str(ref)+'_ref-' + str(Dim) + '_grains-' + str(int(misorientation_list[1]*180/np.pi)) + '_angle-nointer'
if output_flag: 
    os.system("mkdir "+savename)
    matplotlib.use('agg')
    
#########
#########  Variable Decl. ##########
#########
phi = np.zeros((number_of_grain,nx,ny))
phi_new = np.zeros((number_of_grain,nx,ny))
mf = np.zeros((1000,nx,ny),dtype = int)
nf = np.zeros((nx,ny),dtype = int)
wij = np.zeros((number_of_grain,number_of_grain)); aij = np.zeros((number_of_grain,number_of_grain)); mij = np.zeros((number_of_grain,number_of_grain))
eij = np.zeros((number_of_grain,number_of_grain))
gb = np.zeros((nx,ny))
sigma11_R1 = np.zeros((nx,ny),dtype = float);sigma12_R1 = np.zeros((nx,ny),dtype = float);sigma22_R1 = np.zeros((nx,ny),dtype = float)
sigma11_R2 = np.zeros((nx,ny),dtype = float);sigma12_R2 = np.zeros((nx,ny),dtype = float);sigma22_R2 = np.zeros((nx,ny),dtype = float)

for i in range(0,number_of_grain):
    for j in range(0,number_of_grain):
        eij[i,j] = 0.0
        mij[i,j] = pmobi
        if i == j:
            mij[i,j] = 0.0
        if i == 0 or j == 0:
            eij[i,j] = eee
        if i < j:
            eij[i,j] = -eij[i,j]

phi[0,:,:] = 1.0
nf[:,:] = 1

##########
##########  Initial Condition   ##########
##########
print("...Initialization...")
#t0=time.time()
Initialize(phi,gb,nf,mf,nx,ny,dx,dy,eta,number_of_grain,x1_grain,x2_grain,r_nuclei)
plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,0)
idov2,idov3=overlaps(phi,number_of_grain)
phi_old=phi
#print("t:",time.time()-t0)

##########
##########   Compute Stress Field of a Dislocation   ###########
##########
print("...Stress Field single dislocation...")
t0=time.time()
stress_field_bulk_single(phi, nx, ny, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,dx,dy)
print("t:",time.time()-t0)
#########
#########   Time Itetation  ###########
#########
print("...Time Iteration...")
#
t0=time.time()
for nstep in range(1,nsteps+1):
    if nstep% 10 == 0:
        print ("step: ", nstep, "  sim. time: ", time.time()-t0)
        t0=time.time()
    
    sigma11 = np.zeros((nx,ny),dtype = float)
    sigma12 = np.zeros((nx,ny),dtype = float)
    sigma22 = np.zeros((nx,ny),dtype = float)
    sigma11_b = np.zeros((nx,ny),dtype = float)
    sigma12_b = np.zeros((nx,ny),dtype = float)
    sigma22_b = np.zeros((nx,ny),dtype = float)
    
    #t0=time.time()
    if( nstep%50 == 0 ):
    #    if(check_grains(phi,phi_old,number_of_grain)):
        idov2,idov3=overlaps(phi,number_of_grain)
    #        phi_old=phi
        
    gb_x_t, gb_y_t = find_gb(phi, nf, mf, nx, ny, dx, dy, number_of_grain, G, sigma11, sigma12, sigma22, sigma11_b, sigma12_b, sigma22_b, sigma12_ext, sigma11_ext, sigma22_ext, ddd, nstep, stress_step, sigma11_R1, sigma12_R1, sigma22_R1, sigma11_R2, sigma12_R2, sigma22_R2,misorientation_list,ref_theta_i,idov2,idov3)
    #print("     Find gb: ",time.time()-t0,number_of_grain)
    A = np.zeros(number_of_grain)
    #t0=time.time()
    update_nfmf(phi,mf,nf,nx,ny,number_of_grain)
    #print("     Update nfmf: ",time.time()-t0)
    #t0=time.time()
    #update_PF(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt, ETA)
    update_PF(phi,phi_new,nx,ny,dx,dy,g1,g2,eta,ref_theta_i,sigma11,sigma12,sigma22,pmobi,A,mf,nf,eij,misorientation_list,dt)

    #print("     Update PF: ",time.time()-t0)
    #t0=time.time()
    renorm(phi,phi_new,number_of_grain)
    #print("     Renorm: ",time.time()-t0)
    #print("#")

    if output_flag:
        if nstep % n_grain_area == 0:   
            A = update_area(A, phi)
            f = open(savename+'/' + str(Dim)+'-polycrystal_grain-area.txt', 'a')
            for kkk in range(0, number_of_grain):
                    f.write(str(A[kkk]) + '\t')
            f.write('\n')
            f.close()
    
    if nstep % stress_step == 0:
        plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,nstep)
        plotstress(gb_x_t,gb_y_t,sigma12_b, sigma11_b, sigma22_b,savename,output_flag,nstep)
    elif nstep % gb_step == 0:
        plotgb(phi,gb,nx,ny,number_of_grain,savename,output_flag,nstep)
    
