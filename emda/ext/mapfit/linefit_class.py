"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from timeit import default_timer as timer
from emda.plotter import plot_nlines
from emda.quaternions import get_RM, rotationMatrixToEulerAngles
from emda.mapfit.utils import get_FRS,create_xyz_grid,get_xyz_sum
import fcodes_fast
from emda.restools import cut_resolution_for_linefit
from emda.mapfit.derivatives_newmethod2 import derivatives
from emda.config import *

np.set_printoptions(suppress=True) # Suppress insignificant values for clarity

class LineFit:
    def __init__(self):
        self.e0_lf          = None
        self.step           = None
        self.e1_lf          = None
        self.cbin_idx_lf    = None
        self.cbin_lf        = None
        self.w_grid         = None


    def get_linefit_static_data(self,e0,cbin_idx,res_arr,smax):
        self.e0_lf,self.cbin_idx_lf,self.cbin_lf = cut_resolution_for_linefit(e0,
                                                    cbin_idx,
                                                    res_arr,
                                                    smax)

    def f(self,k):
        from emda.mapfit.utils import get_fsc_wght
        # w = 1.0 for line search
        nx,ny,nz = self.e0_lf.shape
        t = self.step[:3]*k[0]
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        q_init = np.array([1.0, 0.0, 0.0, 0.0])
        tmp = np.insert(self.step[3:]*k[1], 0, 0.0)
        tmp = tmp + q_init
        q = tmp/np.sqrt(np.dot(tmp, tmp))
        rotmat = get_RM(q)
        ers = get_FRS(rotmat,self.e1_lf * st, interp='linear')
        w_grid = get_fsc_wght(self.e0_lf, ers[:,:,:,0],self.cbin_idx_lf,self.cbin_lf)
        fval = np.sum(w_grid * self.e0_lf * np.conjugate(ers[:,:,:,0]))
        return -fval.real

    def calc_fval_for_different_kvalues_at_this_step(self,e1):
        from scipy import optimize
        start = timer()
        nx,ny,nz = e1.shape
        cx = self.e0_lf.shape[0] // 2
        cy = cz = cx
        dx = int((nx - 2*cx)/2)
        dy = int((ny - 2*cy)/2)
        dz = int((nz - 2*cz)/2)
        self.e1_lf = e1[dx:dx+2*cx,dy:dy+2*cy,dz:dz+2*cz]
        assert self.e1_lf.shape == self.e0_lf.shape
        init_guess = [1.0,1.0]
        minimum = optimize.minimize(self.f, init_guess, method='Powell')
        end = timer()
        #print(' time for line search: ', end-start)
        return minimum.x