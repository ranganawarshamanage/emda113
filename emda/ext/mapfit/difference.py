# difference map 
from __future__ import absolute_import, division, print_function, unicode_literals
from emda import iotools
from emda.mapfit import mapaverage
import numpy as np
import fcodes_fast
from emda.config import *

def difference_map(f1, f2, cell, origin=None, res_arr=None, bin_idx=None, nbin=None, s_grid=None):
    from emda import fsc, restools
    nx, ny, nz = f1.shape
    if s_grid is None:
        maxbin = np.amax(np.array([nx//2,ny//2,nz//2]))
        _, s_grid, _ = fcodes_fast.resolution_grid_full(cell,0.0,1,maxbin,nx,ny,nz)
    if bin_idx is None:
        nbin,res_arr,bin_idx = restools.get_resolution_array(cell,f1)
    # estimating covariance between current map vs. static map
    f1f2_fsc,_ = fsc.anytwomaps_fsc_covariance(f1,f2,bin_idx,nbin)
    '''_,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                            f1,f2,bin_idx,nbin,0,nx,ny,nz)'''
    # find the relative B-factor between static and frt_full
    resol = mapaverage.get_resol(f1f2_fsc,res_arr)
    scale, bfac_relative = mapaverage.estimate_scale_and_relativebfac([f1,f2],cell,resol)
    #fobj.write('Relative B-factor between static map and moving(full) map: '+ str(bfac_relative)+'\n')
    # apply B-factor and scale to 2nd map
    f2_scaled = scale * f2 * np.exp(bfac_relative * (s_grid**2)/4.0)
    # calculating maps
    dm1_dm2 = np.real(np.fft.ifftshift(
                         np.fft.ifftn(
                         np.fft.ifftshift(f1 - f2_scaled))))
    dm2_dm1 = np.real(np.fft.ifftshift(
                         np.fft.ifftn(
                         np.fft.ifftshift(f2_scaled - f1))))
    iotools.write_mrc(dm1_dm2,
                      'diffmap_m1-m2.mrc',
                      cell, origin)
    iotools.write_mrc(dm2_dm1,
                      'diffmap_m2-m1.mrc',
                      cell, origin)

def diffmap_scalebypower(f1, f2, cell, smax=0.0, origin=None, res_arr=None, \
    bin_idx=None, nbin=None, s_grid=None):
    from emda.scale_maps import scale_twomaps_by_power, transfer_power
    from emda import restools

    if bin_idx is None:
        nbin,res_arr,bin_idx = restools.get_resolution_array(cell,f1)
    
    scale = scale_twomaps_by_power(f1=f1, f2=f2, bin_idx=bin_idx, \
        res_arr=res_arr)
    if smax >= 0.0:
        mask = res_arr > smax
        for i in range(len(res_arr)):
            res = res_arr[i] * mask[i]
            scl = scale[i] * mask[i]
            print(res, scl)
    else:
        print('resolution cannot be negative. Stopping now...')
        exit()

    scale_grid = transfer_power(bin_idx,res_arr,scale * mask)
    f2_scaled = scale_grid * f2
    mask_grid = transfer_power(bin_idx,res_arr,mask)
    f1 = f1 * mask_grid
    # calculating maps
    dm1_dm2 = np.real(np.fft.ifftshift(
                         np.fft.ifftn(
                         np.fft.ifftshift(f1 - f2_scaled))))
    dm2_dm1 = np.real(np.fft.ifftshift(
                         np.fft.ifftn(
                         np.fft.ifftshift(f2_scaled - f1))))
    iotools.write_mrc(dm1_dm2,
                      'diffmap_m1-m2_pwr.mrc',
                      cell, origin)
    iotools.write_mrc(dm2_dm1,
                      'diffmap_m2-m1_pwr.mrc',
                      cell, origin)
    return dm1_dm2, dm2_dm1

def diffmap_scalebyampli(f1, f2, cell, smax=0.0, origin=None, res_arr=None, \
    bin_idx=None, nbin=None):
    from emda import restools
    import fcodes_fast as fc

    if bin_idx is None:
        nbin,res_arr,bin_idx = restools.get_resolution_array(cell,f1)
    nx, ny, nz = f1.shape
    diffmap = fc.differencemap(fo=f1,fc=f2,bin_idx=bin_idx,res_arr=res_arr, \
        smax=smax,mode=1,nbin=len(res_arr),nx=nx,ny=ny,nz=nz)
    return diffmap

def diffmap_normalisedsf(f1, f2, cell, smax=0.0, origin=None, res_arr=None, \
    bin_idx=None, nbin=None):
    from emda import restools
    import fcodes_fast as fc

    if bin_idx is None:
        nbin,res_arr,bin_idx = restools.get_resolution_array(cell,f1)

    nx, ny, nz = f1.shape

    diffmap = fc.diffmap_norm(fo=f1,fc=f2,bin_idx=bin_idx,res_arr=res_arr, \
        smax=smax,mode=1,nbin=len(res_arr),nx=nx,ny=ny,nz=nz)

    return diffmap