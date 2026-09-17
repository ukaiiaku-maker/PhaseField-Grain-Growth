import importlib.util
from pathlib import Path
import sys
import types

import numba as nb
import numpy as np

from grain_growth_pf.pf.anisotropic import _pair_law
from grain_growth_pf.pf.qiu_si import qiu_pair_anisotropy_correction


ROOT=Path(__file__).parents[2]
FUNCTIONS=ROOT/"native_compatibility/qiu_anisotropy_v1/functions_4ref_new.py"


def load_functions():
    sys.modules.setdefault("sparse",types.ModuleType("sparse"))
    spec=importlib.util.spec_from_file_location("qiu_anisotropy_native_functions",FUNCTIONS)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def active_lists(phi):
    phases,nx,ny=phi.shape;nf=np.full((nx,ny),phases,dtype=np.int64);mf=np.zeros((1000,nx,ny),dtype=np.int64)
    for phase in range(phases): mf[phase]=phase
    return nf,mf


def test_a0_control_kernel_is_bitwise_identical_to_modern_native_update():
    nb.set_num_threads(1);f=load_functions();rng=np.random.default_rng(91)
    phi=rng.uniform(0.1,1.0,(3,4,5));phi/=phi.sum(axis=0);nf,mf=active_lists(phi)
    arguments=(4,5,1,1,1.0,1.0,5.0,np.pi/4,np.zeros((4,5)),np.zeros((4,5)),np.zeros((4,5)),2*np.pi**2/(8*5),np.zeros(3),mf,nf,np.asarray([[0.,-20.,-20.],[20.,0.,0.],[20.,0.,0.]]),np.asarray([0.,-0.2,0.4]),0.1)
    native=np.zeros_like(phi);control=np.zeros_like(phi)
    f.update_PF(phi.copy(),native,*arguments)
    f.update_PF_control(phi.copy(),control,*arguments,0,0.65,0.85,16,2.0,1.0910512514090829,1.0,1.0)
    np.testing.assert_array_equal(control,native)
    assert f.update_PF_control.nopython_signatures


def test_local_anisotropic_drive_matches_global_variational_oracle():
    f=load_functions();rng=np.random.default_rng(92);phi=rng.uniform(0.1,0.8,(2,4,5));arguments=(0.1,0.7,1.0,0.65,0.85,16,1.0910512514090829,0.93)
    _,oracle=qiu_pair_anisotropy_correction(phi[0],phi[1],*arguments)
    measured=np.empty((4,5))
    for l in range(4):
        for m in range(5):
            measured[l,m]=f.qiu_anisotropic_pair_drive_cell(phi,0,1,l,m,(l+1)%4,(l-1)%4,(m+1)%5,(m-1)%5,1.0,0.1,0.7,0.65,0.85,16,1.0910512514090829,0.93)
    np.testing.assert_allclose(measured,oracle,rtol=2e-14,atol=2e-14)


def test_native_mobility_kernel_matches_pair_law():
    f=load_functions();rng=np.random.default_rng(93);phi=rng.uniform(0.1,0.8,(2,4,5));l,m=2,3;lp,lm,mp,mm=3,1,4,2
    gix=(phi[0,lp,m]-phi[0,lm,m])/2;giy=(phi[0,l,mp]-phi[0,l,mm])/2;gjx=(phi[1,lp,m]-phi[1,lm,m])/2;gjy=(phi[1,l,mp]-phi[1,l,mm])/2
    theta=np.arctan2(phi[0,l,m]*gjy-phi[1,l,m]*giy,phi[0,l,m]*gjx-phi[1,l,m]*gix)
    pmobi=2*np.pi**2/(8*5)
    expected=_pair_law(theta,0.1,0.7,1.0,pmobi,0.65,0.85,16,2.0,1.0910512514090829,0.93,1.07)[2]
    actual=f.qiu_anisotropic_mobility(phi,0,1,l,m,lp,lm,mp,mm,1.0,0.1,0.7,pmobi,0.65,0.85,16,2.0,1.0910512514090829,0.93,1.07)
    np.testing.assert_allclose(actual,expected,rtol=2e-15,atol=2e-15)
