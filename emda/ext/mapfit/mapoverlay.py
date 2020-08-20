"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from timeit import default_timer as timer
from emda.mapfit import map_Class
from emda.mapfit import emfit_Class
from emda.mapfit import linefit_class
from emda.config import *

def output_rotated_maps(emmap1, r_lst, t_lst, Bf_arr=[0.0]):
    import numpy as np
    import fcodes_fast
    from emda.mapfit import utils
    from emda.plotter import plot_nlines,plot_nlines_log
    from scipy.ndimage.interpolation import shift
    from emda import iotools 
    import numpy.ma as ma
    from emda import fsc

    #com = emmap1.com1  
    arr_lst = emmap1.arr_lst
    fo_lst = emmap1.fo_lst   
    bin_idx = emmap1.bin_idx  
    nbin = emmap1.nbin     
    res_arr = emmap1.res_arr  
    cell = emmap1.map_unit_cell 
    origin = emmap1.map_origin
    #assert len(fo_lst) == len(t_lst) == len(r_lst)
    nx,ny,nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])  
    cov_lst = []
    fsc12_lst = []
    fsc12_lst_unaligned = []
    imap_f = 0
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    #data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    #iotools.write_mrc(arr_lst[0],'static_map_realsp.mrc',cell,origin)
    del data2write
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        f1f2_fsc_unaligned,_ = fsc.anytwomaps_fsc_covariance(fo_lst[0],fo,bin_idx,nbin)
        fsc12_lst_unaligned.append(f1f2_fsc_unaligned)            
        imap_f = imap_f + 1
        '''# rotate original map
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        arr = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo * st))))
        #arr = np.fft.fftshift(arr)
        print('before triliner map')
        nx, ny, nz = arr.shape
        data2write = fcodes_fast.trilinear_map(rotmat,arr,nx,ny,nz)
        print('after triliner map')'''
        #
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        frt_full = utils.get_FRS(rotmat,fo * st, interp='cubic')[:,:,:,0]
        frt_lst.append(frt_full)
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        #data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        iotools.write_mrc(data2write,"{0}_{1}.{2}".format('fitted_map',str(imap_f),'mrc'),cell,origin)      
        # estimating covaraince between current map vs. static map
        f1f2_fsc,f1f2_covar = fsc.anytwomaps_fsc_covariance(fo_lst[0],frt_full,bin_idx,nbin)
        cov_lst.append(f1f2_covar)
        fsc12_lst.append(f1f2_fsc)
        plot_nlines(res_arr,
                    [f1f2_fsc_unaligned,f1f2_fsc],
                    "{0}_{1}.{2}".format('fsc',str(imap_f),'eps'),
                    ["FSC before","FSC after"])

def main(maplist, ncycles, t_init, theta_init, smax, masklist, fobj, \
        interp, halfmaps=False):
    from emda.quaternions import get_quaternion, get_RM
    import numpy as np
    from emda.mapfit.utils import create_xyz_grid
    from emda.mapfit import interp_derivatives
    #from emda import scale_maps
    #from emda.iotools import write_mrc
    from emda.mapfit import run_fit
    if len(maplist) < 2: 
        print(" At least 2 maps required!")
        exit()
    try:
        if halfmaps:
            print('Map overlay using halfmaps')
            emmap1 = map_Class.Overlay(maplist,masklist) #halfmaps to overlay
        else:
            print('Map overlay not using halfmaps')
            emmap1 = map_Class.EmmapOverlay(maplist,masklist) #maps to overlay
        fobj.write('Map overlay\n')
    except NameError:
        if halfmaps: 
            print('Map overlay using halfmaps')
            emmap1 = map_Class.Overlay(maplist)
        else:
            print('Map overlay not using halfmaps')
            emmap1 = map_Class.EmmapOverlay(maplist)
        fobj.write('Map overlay\n')
    emmap1.load_maps(fobj)
    emmap1.calc_fsc_from_maps(fobj)
    # converting theta_init to rotmat for initial iteration
    fobj.write('\n')
    fobj.write('Initial fitting parameters:\n')
    fobj.write('    Translation: '+ str(t_init)+' \n')
    fobj.write('    Rotation: '+ str(theta_init)+' \n')
    q = get_quaternion(theta_init)
    q = q/np.sqrt(np.dot(q,q))
    print('Initial quaternion: ', q)
    rotmat = get_RM(q)
    fobj.write('    Rotation matrix: '+ str(rotmat)+' \n')
    fobj.write('\n')
    fobj.write('    # fitting cycles: '+ str(ncycles)+' \n')
    t = t_init
    rotmat_lst = []
    transl_lst = []
    # resolution estimate for line-fit
    dist = np.sqrt((emmap1.res_arr - smax)**2) # smax can be a reasonable value
    slf = np.argmin(dist) + 1 
    if slf % 2 != 0: slf = slf - 1
    slf = min([len(dist), slf])
    #
    dfs_interp = False
    # preparing parameters for minimization
    if dfs_interp:
        cell = emmap1.map_unit_cell
        xyz = create_xyz_grid(cell, emmap1.map_dim)
        vol = cell[0] * cell[1] * cell[2]
    #
    for ifit in range(1, len(emmap1.eo_lst)):
        fobj.write('Fitting set#: ' + str(ifit)+' \n')
        start_fit = timer()
        if dfs_interp:
            #dfs_full = interp_derivatives.dfs_fullmap(emmap1.arr_lst[ifit],xyz,vol)
            dfs_full = interp_derivatives.dfs_fullmap(
                                     np.fft.ifftshift(
                                     np.real(np.fft.ifftn(
                                     #np.fft.ifftshift(f1)))),
                                     np.fft.ifftshift(emmap1.eo_lst[ifit])))),
                                                      xyz,vol)
        rotmat, t = run_fit.run_fit(emmap1,smax,rotmat,t,slf,ncycles,ifit,fobj,interp)
        rotmat_lst.append(rotmat)
        transl_lst.append(t)
        end_fit = timer()
        fobj.write('Final Translation: '+ str(t)+' \n')
        fobj.write('Final Rotation matrix: '+ str(rotmat)+' \n')
        print('time for fitting: ', end_fit-start_fit)

    output_rotated_maps(emmap1,
                        rotmat_lst,
                        transl_lst)

