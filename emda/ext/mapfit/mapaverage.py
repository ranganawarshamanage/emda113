"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from emda.mapfit import map_Class
from emda.mapfit import emfit_Class
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

    com = emmap1.com1  
    fo_lst = emmap1.fo_lst   
    eo_lst = emmap1.eo_lst
    hffsc_lst = emmap1.hffsc_lst 
    signalvar_lst = emmap1.signalvar_lst
    totalvar_lst = emmap1.totalvar_lst
    bin_idx = emmap1.bin_idx  
    nbin = emmap1.nbin     
    res_arr = emmap1.res_arr  
    cell = emmap1.map_unit_cell 
    origin = emmap1.map_origin
    #assert len(fo_lst) == len(t_lst) == len(r_lst)
    nx,ny,nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])  
    ert_lst = []
    ert_lst.append(eo_lst[0])  
    cov_lst = []
    fsc12_lst = []
    imap_f = 0
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    del data2write
    for fo, eo, t, rotmat in zip(fo_lst[1:], eo_lst[1:], t_lst, r_lst):
        imap_f = imap_f + 1
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        frt_full = utils.get_FRS(cell,rotmat,fo * st)[:,:,:,0]
        frt_lst.append(frt_full)
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        iotools.write_mrc(data2write,
                          "{0}_{1}.{2}".format('fitted_map',str(1),'mrc'),
                          cell,
                          origin)
        # estimating covaraince between current map vs. static map
        f1f2_fsc,f1f2_covar = fsc.anytwomaps_fsc_covariance(fo_lst[0],
                                                            frt_full,
                                                            bin_idx,
                                                            nbin)
        '''f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full,bin_idx,nbin,0,nx,ny,nz)'''
        # mask f1f2_covar so that anything beyond zero get zeros
        f1f2_covar = set_array(f1f2_covar)
        cov_lst.append(f1f2_covar)
        fsc12_lst.append(f1f2_fsc)
        # normalised maps
        ert_lst.append(utils.get_FRS(cell,rotmat,eo * st)[:,:,:,0])

    # Here calculates the new singal_vars using covariance
    static_signal = np.zeros((nbin,len(cov_lst)), dtype='float')
    pred_signalvar_lst = []
    # 1. fit a curve to S-matrix and get B value of that curve
    for i in range(len(cov_lst)):
        s_0i_ma = ma.masked_less_equal(cov_lst[i],0.0)
        # apply s_0i mask on sgnl_var_0i
        sgnl_var_0i_ma = ma.masked_less_equal(signalvar_lst[0] * 
                                              signalvar_lst[i+1], 
                                              0.0)
        sgnl_var_0i_ma = ma.masked_where(ma.getmask(s_0i_ma), sgnl_var_0i_ma)
        s_0i = (s_0i_ma/ma.masked_where(ma.getmask(s_0i_ma), 
                        np.sqrt(sgnl_var_0i_ma.filled(0.0)))).filled(0.0)
        bfac, intercept = get_bfac(res_arr,s_0i)
        s = ((1.0/res_arr)**2)/4.0
        pred_s01 =  np.exp(-1.0*abs(bfac) * s + intercept)
        plot_nlines(res_arr,[s_0i, pred_s01],'s_0i.eps',["s_0i","pred_s_0i"])        
        # 2. estimate scale and relative B value between maps
        resol = get_resol(fsc12_lst[i],res_arr)
        scale, bfac_relative = estimate_scale_and_relativebfac([frt_lst[0],frt_lst[i+1]],cell,resol)
        # 3. now predict new signalvar of ii maps using cov_ij
        signal_00 = cov_lst[0]*np.exp(abs(bfac)*s) * np.sqrt(scale*np.exp(bfac_relative*s))
        static_signal[:,i] = signal_00
        if i == 0:
            pred_signalvar_lst.append(signal_00)
        signal_11 = cov_lst[0]*np.exp(abs(bfac)*s) / np.sqrt(scale*np.exp(bfac_relative*s))
        pred_signalvar_lst.append(signal_11)
        # 4. now plot new and old signals
        '''signal_lst2plot = [ signalvar_lst[0],
                            signalvar_lst[i+1],
                            cov_lst[i],
                            signal_00,
                            signal_11 ]
        plot_nlines_log(res_arr,
                        signal_lst2plot,
                        ["old_s00","old_s11","s01","new_s00","new_s11"],
                        'log_variance_all.eps')'''

    pred_signalvar_lst[0] = np.average(static_signal, axis=1)
    plot_nlines_log(res_arr,
                    pred_signalvar_lst)    

    # output sigma values into a file
    '''for i in range(len(signal_00)):
        s1 = signal_00[i]
        s2 = signal_11[i]
        s12 = cov_lst[0][i]
        print(s1, s2, s12, s12**2 - s1*s2)'''

    #pred_signalvar_lst = [signal_00,signal_11]
    # map averaging with normalised maps
    averagemaps = utils.avg_and_diffmaps(ert_lst,
                                         cell,nbin,
                                         pred_signalvar_lst,
                                         totalvar_lst,
                                         cov_lst,
                                         hffsc_lst,
                                         bin_idx,
                                         emmap1.s_grid,
                                         res_arr,
                                         Bf_arr)
    # calc. fsc between static map and averaged maps
    for i in range(2):
        f1f2_fsc,_ = fsc.anytwomaps_fsc_covariance(eo_lst[0],
                                                   averagemaps[:,:,:,i,0],
                                                   bin_idx,
                                                   nbin)
        '''_,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                            eo_lst[0],averagemaps[:,:,:,i,0],bin_idx,nbin,0,nx,ny,nz)'''
        fsc12_lst.append(f1f2_fsc)
    '''plot_nlines(res_arr,
                fsc12_lst,
                'fsc_compare.eps',
                ["static-aligned2","static-avg1","static-avg2"])     ''' 
    #utils.output_maps(averagemaps,emmap1.com1,cell,origin,Bf_arr)


'''def difference_map(emmap1, r_lst, t_lst, Bf_arr=[0.0]):
    # calculate only the difference map
    from emda.mapfit import diffmap
    #diffmap.calc_diffmap(emmap1, r_lst, t_lst, Bf_arr=[0.0])
    diffmap.calc_diffmap_simple(emmap1, r_lst, t_lst, Bf_arr=[0.0])'''

def get_resol(bin_fsc,res_arr):
    import numpy as np
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    if len(bin_fsc) > 0: dist = np.sqrt((bin_fsc - 0.143)**2)
    return res_arr[np.argmin(dist)]

def get_bfac(res_arr,signal,scale=1):
    from scipy import stats
    import numpy as np
    high_res = 5.5 # Angstrom. This is a hard cutoff. Need to change later
    res_ub = np.sqrt((res_arr - high_res)**2)
    ub_res = np.argmin(res_ub)
    res_arr_trunc, signal_trunc = trunc_array(res_arr,signal)
    if len(signal_trunc) < ub_res:
        s = 1/res_arr[:ub_res]
        slope, intercept,_,_,_ = stats.linregress((s*s)/4.0, np.log(signal[:ub_res]))
    else:
        s = 1/res_arr_trunc
        slope, intercept,_,_,_ = stats.linregress((s*s)/4.0, np.log(signal_trunc))
    bfac = slope
    print('B factor = ', bfac, 'Intercept = ', intercept)    
    return bfac, intercept

def estimate_scale_and_relativebfac(fo_lst,uc,resol):
    from emda.mapfit import curve_fit_3
    params = curve_fit_3.main(fo_lst,uc,resol)
    scale, bfac_relative = params
    return scale, bfac_relative

def set_array(arr, thresh=0.0):
    set2zero = False
    i = -1
    for ival in arr:
        i = i + 1
        if ival <= thresh and not set2zero:
            set2zero = True
        if set2zero: arr[i] = 0.0
    return arr

def trunc_array(res_arr,signal):
    import numpy as np
    val_lst = []
    res_lst = []
    j = -1
    fsc_thresh = 0.7 # this is an arbitrary choice.
    for ival in signal:
        if ival < 1.0:
            if ival > fsc_thresh:
                j = j + 1
                val_lst.append(ival)
                res_lst.append(res_arr[j])
            if ival < fsc_thresh:
                break
    return np.array(res_lst, dtype='float', copy=False), np.array(val_lst, dtype='float', copy=False)

def fsc_between_static_and_transfomed_map(staticmap,movingmap,bin_idx,rm,t,cell,nbin):
    import fcodes_fast
    from emda.mapfit import utils
    from emda import fsc
    nx, ny, nz = staticmap.shape
    st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
    frt_full = utils.get_FRS(cell,rm,movingmap * st)[:,:,:,0]
    f1f2_fsc,_ = fsc.anytwomaps_fsc_covariance(staticmap,
                                               frt_full,
                                               bin_idx,
                                               nbin)
    '''f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                staticmap,frt_full,bin_idx,nbin,0,nx,ny,nz)'''
    return f1f2_fsc

def get_ibin(bin_fsc):
    import numpy as np
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    #dist = np.sqrt((bin_fsc - 0.143)**2)
    dist = np.sqrt((bin_fsc - 0.3)**2)
    ibin = np.argmin(dist) + 1 # adding 1 because fResArr starts with zero
    return ibin

def main(maplist, ncycles, t_init, theta_init, smax, masklist, \
    fobj, interp, fit = False, Bf_arr=[0.0]):
    import numpy as np
    from emda import fsc
    from emda.quaternions import get_quaternion, get_RM
    from emda.mapfit import newsignal, run_fit, utils

    if len(maplist) < 4: 
        print('Map averaging requires at least two sets of maps ' 
             'with half sets')
        exit()
    elif (len(maplist) >= 4 and len(maplist) % 2 == 0):
        print('Map average enabled')
        fobj.write('Map average enabled\n')
    try:
        emmap1 = map_Class.EmmapAverage(hfmap_list=maplist,mask_list=masklist)
    except NameError:
        emmap1 = map_Class.EmmapAverage(hfmap_list=maplist)
    emmap1.load_maps(fobj)
    emmap1.calc_fsc_variance_from_halfdata(fobj)
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
    if fit:
        for ifit in range(1, len(emmap1.eo_lst)):
            fobj.write('Fitting set#: ' + str(ifit)+' \n')
            rotmat, t = run_fit.run_fit(emmap1,smax,rotmat,t,slf,ncycles,ifit,fobj,interp)
            rotmat_lst.append(rotmat)
            transl_lst.append(t)

        '''output_rotated_maps(emmap1,
                            rotmat_lst,
                            transl_lst)'''
        '''newsignal.output_rotated_maps(emmap1,
                            rotmat_lst,
                            transl_lst,
                            fobj=fobj)'''
        '''newsignal.output_rotated_maps_using_fo(emmap1,
                            rotmat_lst,
                            transl_lst,
                            fobj=fobj)'''
        newsignal.output_rotated_maps_using_fo_allmaps(emmap1,
                            rotmat_lst,
                            transl_lst,
                            fobj=fobj)
    if not fit:
        '''# using normalized Fourier coeffi
        cov_lst = []
        # estimating covaraince between current map vs. static map
        for frt in emmap1.fo_lst[1:]:
            f1f2_fsc,f1f2_covar = \
                fsc.anytwomaps_fsc_covariance(emmap1.fo_lst[0], \
                    frt,emmap1.bin_idx, emmap1.nbin)
            # mask f1f2_covar so that anything beyond zero get zeros
            f1f2_covar = set_array(f1f2_covar)
            cov_lst.append(f1f2_covar)        
        averagemaps = utils.avg_and_diffmaps(maps2avg=emmap1.eo_lst,
                                        uc=emmap1.map_unit_cell,
                                        nbin=emmap1.nbin,
                                        sgnl_var=emmap1.signalvar_lst,
                                        totl_var=emmap1.totalvar_lst,
                                        covar=cov_lst,
                                        hffsc=emmap1.hffsc_lst,
                                        bin_idx=emmap1.bin_idx,
                                        s_grid=emmap1.s_grid,
                                        res_arr=emmap1.res_arr,
                                        Bf_arr=Bf_arr)'''
        # using un-normalized Fourier coeffi.
        averagemaps = newsignal.avgmaps_using_fo_nofit(emmap1=emmap1, \
            fobj=fobj)
        # averagemap output
        utils.output_maps(averagemaps=averagemaps,com_lst=[], \
            t_lst=[], r_lst=[], unit_cell=emmap1.map_unit_cell, \
                map_origin=emmap1.map_origin,bf_arr=Bf_arr, center=None,)

def run_fit(emmap1, f1, f2, t_init, theta_init, smax, ncycles, fobj):
    import numpy as np
    from emda.mapfit.frequency_marching import frequency_marching
    from emda.quaternions import get_quaternion, get_RM
    fsc_lst = []
    q = get_quaternion(theta_init)
    q = q/np.sqrt(np.dot(q,q))
    print(q)
    rotmat = get_RM(q)
    fobj.write('    Rotation matrix: '+ str(rotmat)+' \n')
    fobj.write('\n')
    fobj.write('    # fitting cycles: '+ str(ncycles)+' \n')
    t = t_init
    fobj.write('Fitting set#:  \n')
    for i in range(5):
        if i==0:
            #smax = 6.0 # A
            dist = np.sqrt((emmap1.res_arr - smax)**2)
            ibin = np.argmin(dist) + 1 # adding 1 because fResArr starts with zero
            ibin_old = ibin
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                                        f1,
                                        f2,
                                        emmap1.bin_idx,
                                        rotmat,
                                        t,
                                        emmap1.map_unit_cell,
                                        emmap1.nbin)
            fsc_lst.append(f1f2_fsc)
        else:
            # Apply initial rotation and translation to calculate fsc
            f1f2_fsc = fsc_between_static_and_transfomed_map(
                                        f1,
                                        f2,
                                        emmap1.bin_idx,
                                        rotmat,
                                        t,
                                        emmap1.map_unit_cell,
                                        emmap1.nbin)
            ibin = get_ibin(f1f2_fsc)
            if ibin_old == ibin: 
                fsc_lst.append(f1f2_fsc)
                fobj.write('\n')
                fobj.write('FSC between static and moving maps\n')
                fobj.write('\n')
                fobj.write('bin#\n')
                fobj.write('resolution (A)\n')
                fobj.write('start FSC\n')
                fobj.write('end FSC\n')
                fobj.write('\n')
                for j in range(len(emmap1.res_arr)):
                    #print(emmap1.res_arr[j], fsc_lst[0][j], fsc_lst[1][j])
                    fobj.write("{:5d} {:6.2f} {:8.4f} {:8.4f}\n".format(
                        j,
                        emmap1.res_arr[j],
                        fsc_lst[0][j],
                        fsc_lst[1][j]))
                break
            else: 
                ibin_old = ibin
        static_cutmap,cBIdx, cbin = frequency_marching(f1,
                                                       emmap1.bin_idx,
                                                       emmap1.res_arr,
                                                       bmax=ibin,
                                                       fobj=fobj)  
        moving_cutmap,_, _ = frequency_marching(f2,
                                                emmap1.bin_idx,
                                                emmap1.res_arr,
                                                bmax=ibin)
        assert static_cutmap.shape == moving_cutmap.shape
        emmap1.ceo_lst = [static_cutmap,moving_cutmap]
        emmap1.cbin_idx = cBIdx
        emmap1.cdim = moving_cutmap.shape
        emmap1.cbin = cbin
        fit = emfit_Class.EmFit(emmap1)
        fit.minimizer(ncycles, t, rotmat, fobj=fobj)
        ncycles = 3
        t = fit.t_accum
        rotmat = fit.rotmat
    return fit.rotmat, fit.t_accum


if (__name__ == "__main__"):
    maplist = ['/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/nat_hf1.mrc',
            '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/nat_hf2.mrc',
            '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/lig_hf1.mrc',
            '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/lig_hf2.mrc']
    masklist = ['/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/nat_msk.mrc',
                '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/nat_msk.mrc',
                '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/lig_msk.mrc',
                '/Users/ranganaw/MRC/REFMAC/Vinoth/reboxed_maps/lig_msk.mrc']
    main(maplist)


