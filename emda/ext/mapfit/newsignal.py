# Applying B-factor and calculating signal
from __future__ import absolute_import, division, print_function, unicode_literals

from emda.mapfit import mapaverage
from emda.config import *

def output_rotated_maps(emmap1, r_lst, t_lst, Bf_arr=[0.0], fobj=None):
    import numpy as np
    import fcodes_fast
    from emda.mapfit import utils
    from emda.plotter import plot_nlines,plot_nlines_log
    from scipy.ndimage.interpolation import shift
    from emda import iotools
    import numpy.ma as ma
    from emda import fsc

    com_lst = emmap1.com_lst 
    fo_lst = emmap1.fo_lst   
    eo_lst = emmap1.eo_lst
    #hffsc_lst = emmap1.hffsc_lst 
    #signalvar_lst = emmap1.signalvar_lst
    #totalvar_lst = emmap1.totalvar_lst
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    #signalvar_lst = []
    #signalvar_lst.append(emmap1.signalvar_lst[0])
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
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
    #fsc12_lst = []
    i = -1
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    fobj.write('\n static_map.mrc was written. \n')
    del data2write
    signalvar_lst = []
    signalvar_static = emmap1.signalvar_lst[0]
    fobj.write('\n ***Map averaging using normalized structure factors*** \n')
    # smoothening signal
    fobj.write('\nStatic map signal \n')
    signalvar_static = get_extended_signal(res_arr,
                                           mapaverage.set_array(signalvar_static),
                                           hffsc_lst[0],
                                           fobj=fobj)
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    for fo, eo, t, rotmat in zip(fo_lst[1:], eo_lst[1:], t_lst, r_lst):
        i = i + 1
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        #frt_full = utils.get_FRS(cell,rotmat,fo * st)[:,:,:,0]
        #frt_lst.append(frt_full)
        '''data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        iotools.write_mrc(data2write,
                          "{0}_{1}.{2}".format('fitted_map',str(1),'mrc'),
                          cell,
                          origin)'''

        # estimating covaraince between current map vs. static map
        '''f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full,bin_idx,nbin,0,nx,ny,nz)
        # mask f1f2_covar so that anything beyond zero get zeros
        f1f2_covar = mapaverage.set_array(f1f2_covar)
        cov_lst.append(f1f2_covar)
        fsc12_lst.append(f1f2_fsc)
        # normalised maps
        ert_lst.append(utils.get_FRS(cell,rotmat,eo * st)[:,:,:,0])'''

        ######
        # apply transformation on half maps of moving map
        fobj.write('Apply transformation on half maps... \n')
        frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i+2] * st)[:,:,:,0]
        frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i+3] * st)[:,:,:,0]
        #frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+2] * st)[:,:,:,0]
        #frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+3] * st)[:,:,:,0]
        frt_full = (frt_hf1 + frt_hf2)/2.0
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
        iotools.write_mrc(data2write,
                          "{0}_{1}.{2}".format('fitted_map',str(i),'mrc'),
                          cell,
                          origin)
        '''# re-estimate sigvar and totvar
        bin_fsc,noisevar,signalvar,totalvar,frt_full,ert_full = fsc.halfmaps_fsc_variance(
                                            frt_hf1,frt_hf2,bin_idx,nbin)
        #frt_full,ert_full,noisevar,signalvar,totalvar,bin_fsc = fcodes_fast.calc_fsc_using_halfmaps(
        #        frt_hf1,frt_hf2,bin_idx,nbin,debug_mode,nx,ny,nz) 
        ert_lst.append(ert_full)
        #frt_lst.append(frt)
        signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))
        hffsc_lst.append(mapaverage.set_array(bin_fsc))
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        iotools.write_mrc(data2write,
                          "{0}_{1}.{2}".format('fitted_map',str(i),'mrc'),
                          cell,
                          origin)'''
        ######

        # estimating covariance between current map vs. static map
        _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full,bin_idx,nbin,0,nx,ny,nz)
        # find the relative B-factor between static and frt_full
        resol = mapaverage.get_resol(f1f2_fsc,res_arr)
        scale, bfac_relative = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_full],cell,resol)
        fobj.write('Relative B-factor between static map and moving(full) map: '+ str(bfac_relative)+'\n')
        # apply B-factor and scale to 2nd map
        s_grid = (emmap1.s_grid**2)/4.0
        frt_full_scaled = scale * frt_full * np.exp(bfac_relative * s_grid)
        # re-estimate covariance
        f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full_scaled,bin_idx,nbin,0,nx,ny,nz)
        # smoothen signal
        fobj.write('\nCovariance between static and moving maps \n')
        f1f2_covar = get_extended_signal(res_arr,mapaverage.set_array(f1f2_covar),f1f2_fsc,fobj=fobj)
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        #fsc12_lst.append(f1f2_fsc)
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf1],cell,resol)
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf2],cell,resol)
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        bin_fsc,noisevar,signalvar,totalvar,frt_full_scaled,ert_full = fsc.halfmaps_fsc_variance(
                                            frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin)
        #_,ert_full,noisevar,signalvar,totalvar,bin_fsc = fcodes_fast.calc_fsc_using_halfmaps(
        #        frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin,debug_mode,nx,ny,nz)
        ert_lst.append(ert_full)
        # smoothen signal
        bin_fsc_full = 2*bin_fsc/(1+bin_fsc)
        fobj.write('\nMoving map signal \n')
        signalvar = get_extended_signal(res_arr,mapaverage.set_array(signalvar),bin_fsc_full,fobj=fobj)       
        signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))
        hffsc_lst.append(mapaverage.set_array(bin_fsc_full))

    var_cov_lst = signalvar_lst + cov_lst
    plot_nlines_log(res_arr,
                    var_cov_lst)
    fvarcov = open("var_cov.txt", "w")
    for i in range(nbin):
        s11 = var_cov_lst[0][i]
        s22 = var_cov_lst[1][i]
        s12 = var_cov_lst[2][i]
        fvarcov.write("{:.2f} {:.4f} {:.4f} {:.4f}\n".format(res_arr[i],s11, s22, s12))
    # map averaging with normalised maps
    averagemaps = utils.avg_and_diffmaps(ert_lst,
                                         cell,nbin,
                                         signalvar_lst,
                                         totalvar_lst,
                                         cov_lst,
                                         hffsc_lst,
                                         bin_idx,
                                         emmap1.s_grid,
                                         res_arr,
                                         Bf_arr)
    # calc. fsc between static map and averaged maps
    '''for i in range(2):
        _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                            eo_lst[0],averagemaps[:,:,:,i,0],bin_idx,nbin,0,nx,ny,nz)
        fsc12_lst.append(f1f2_fsc)
    plot_nlines(res_arr,
                fsc12_lst,
                'fsc_compare.eps',
                ["static-aligned2","static-avg1","static-avg2"])     ''' 
    #utils.output_maps(averagemaps,emmap1.com1,cell,origin,Bf_arr)
    utils.output_maps(averagemaps,emmap1.com_lst,emmap1.box_centr,cell,origin,Bf_arr)

def output_rotated_maps_using_fo(emmap1, r_lst, t_lst, Bf_arr=[0.0], fobj=None):
    import numpy as np
    import fcodes_fast
    from emda.mapfit import utils
    from emda.plotter import plot_nlines,plot_nlines_log
    from scipy.ndimage.interpolation import shift
    from emda import iotools
    import numpy.ma as ma

    com_lst = emmap1.com_lst
    fo_lst = emmap1.fo_lst   
    signalvar_lst = []
    #signalvar_lst.append(emmap1.signalvar_lst[0])
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx  
    nbin = emmap1.nbin     
    res_arr = emmap1.res_arr  
    cell = emmap1.map_unit_cell 
    origin = emmap1.map_origin
    nx,ny,nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0])   
    cov_lst = []
    fsc12_lst = []
    i = -1
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    del data2write
    fobj.write('\n *** Map averaging using un-normalized structure factors *** \n')
    #new_signalvar_lst = []
    #new_signalvar_lst.append(mapaverage.set_array(signalvar_lst[0]))
    #new_totalvat_lst = []
    #new_totalvat_lst.append(mapaverage.set_array(totalvar_lst[0]))
    maps2average_lst = []
    maps2average_lst.append(fo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    # smoothening signal
    fobj.write('\nStatic map signal \n')
    signalvar_static = get_extended_signal(res_arr,
                                           mapaverage.set_array(signalvar_static),
                                           hffsc_lst[0],
                                           fobj=fobj)
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        i = i + 1
        '''data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo))))
        data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        iotools.write_mrc(data2write,'fitted_map_before.mrc',cell,origin)
        del data2write'''
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        frt_full = utils.get_FRS(cell,rotmat,fo * st)[:,:,:,0]
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
        iotools.write_mrc(data2write,'fitted_map_'+str(i+1)+'_.mrc',cell,origin)
        del data2write
        # estimating covaraince between current map vs. static map
        _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full,bin_idx,nbin,0,nx,ny,nz)
        # find the relative B-factor between static and frt_full
        resol = mapaverage.get_resol(f1f2_fsc,res_arr)
        '''scale, bfac_relative = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_full],cell,resol)
        # apply B-factor and scale to 2nd map
        s_grid = (emmap1.s_grid**2)/4.0
        frt_full_scaled = scale * frt_full * np.exp(bfac_relative * s_grid)
        maps2average_lst.append(frt_full_scaled)
        # re-estimate covariance with scaled data
        f1f2_covar,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full_scaled,bin_idx,nbin,0,nx,ny,nz)
        f1f2_covar = mapaverage.set_array(f1f2_covar)
        cov_lst.append(f1f2_covar)'''
        ######
        # apply transformation on half maps of moving map
        '''frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+2] * st)[:,:,:,0]
        frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+3] * st)[:,:,:,0]'''
        frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i+2] * st)[:,:,:,0]
        frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i+3] * st)[:,:,:,0]
        s_grid = (emmap1.s_grid**2)/4.0
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf1],cell,resol)
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf2],cell,resol)
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        # re-estimate sigvar and totvar
        bin_fsc,noisevar,signalvar,totalvar,frt_full_scaled,_ = fsc.halfmaps_fsc_variance(
                                            frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin)
        #frt_full_scaled,_,noisevar,signalvar,totalvar,bin_fsc = fcodes_fast.calc_fsc_using_halfmaps(
        #        frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin,debug_mode,nx,ny,nz) 
        maps2average_lst.append(frt_full_scaled)
        #data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full_scaled))))
        #data2write = shift(data2write, np.subtract(com,emmap1.box_centr))
        #iotools.write_mrc(data2write,'fitted_map_scaled.mrc',cell,origin)
        #del data2write
        # smoothen signal
        bin_fsc_full = 2*bin_fsc/(1+bin_fsc)
        fobj.write('\nMoving map signal \n')
        signalvar = get_extended_signal(res_arr,mapaverage.set_array(signalvar),bin_fsc_full,fobj=fobj)       
        signalvar_lst.append(mapaverage.set_array(signalvar))

        f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full_scaled,bin_idx,nbin,0,nx,ny,nz)
        #f1f2_covar = mapaverage.set_array(f1f2_covar)
        #cov_lst.append(f1f2_covar)
        # smoothen signal
        fobj.write('\nCovariance between static and moving maps \n')
        f1f2_covar = get_extended_signal(res_arr,mapaverage.set_array(f1f2_covar),f1f2_fsc,fobj=fobj)
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        #signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))
        ######

    var_cov_lst = signalvar_lst + cov_lst
    plot_nlines_log(res_arr,
                    var_cov_lst)
    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple(maps2average_lst,
                                                cell, nbin,
                                                signalvar_lst,
                                                totalvar_lst,
                                                cov_lst,
                                                bin_idx,
                                                emmap1.s_grid,
                                                res_arr,
                                                Bf_arr,
                                                emmap1.com_lst,
                                                emmap1.box_centr,
                                                origin)
    # calc. fsc between static map and averaged maps
    '''for i in range(2):
        _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                            fo_lst[0],averagemaps[:,:,:,i,0],bin_idx,nbin,0,nx,ny,nz)
        fsc12_lst.append(f1f2_fsc)
    plot_nlines(res_arr,
                fsc12_lst,
                'fsc_compare.eps',
                ["static-aligned2","static-avg1","static-avg2"])     ''' 
    #utils.output_maps(averagemaps,emmap1.com_lst,emmap1.box_centr,cell,origin,Bf_arr)
    utils.output_maps(averagemaps,emmap1.com_lst,t_lst,r_lst,emmap1.box_centr,cell,origin,Bf_arr)

def output_rotated_maps_using_fo_allmaps(emmap1, r_lst, t_lst, Bf_arr=[0.0], fobj=None):
    import numpy as np
    import fcodes_fast
    from emda.mapfit import utils
    from emda.plotter import plot_nlines,plot_nlines_log
    from scipy.ndimage.interpolation import shift
    from emda import iotools
    import numpy.ma as ma
    from emda import fsc

    com_lst = emmap1.com_lst
    fo_lst = emmap1.fo_lst 
    nmaps = len(fo_lst)  
    signalvar_lst = []
    totalvar_lst = []
    totalvar_lst.append(emmap1.totalvar_lst[0])
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx  
    nbin = emmap1.nbin     
    res_arr = emmap1.res_arr  
    cell = emmap1.map_unit_cell 
    origin = emmap1.map_origin
    nx,ny,nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0]) 
    ert_lst = []
    ert_lst.append(emmap1.eo_lst[0])  
    cov_lst = []
    #fsc12_lst = []
    i = 0
    var_covar_mat = np.zeros((nmaps, nmaps, nbin), dtype='float')
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    del data2write
    fobj.write('\n *** Map averaging using un-normalized structure factors *** \n')
    maps2average_lst = []
    maps2average_lst.append(fo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    '''# smoothening signal
    fobj.write('\nStatic map signal \n')
    signalvar_static = get_extended_signal(res_arr,
                                           mapaverage.set_array(signalvar_static),
                                           hffsc_lst[0],
                                           fobj=fobj)'''
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    var_covar_mat[0,0,:] = mapaverage.set_array(signalvar_static)
    s_grid = (emmap1.s_grid**2)/4.0

    for fo, t, rotmat in zip(fo_lst[1:], t_lst, r_lst):
        #from emda.mapfit import mapaverage
        i = i + 1
        st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,t)
        frt_full = utils.get_FRS(cell,rotmat,fo * st)[:,:,:,0]
        data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(frt_full))))
        data2write = shift(data2write, np.subtract(com_lst[0],emmap1.box_centr))
        iotools.write_mrc(data2write,'fitted_map_'+str(i)+'_.mrc',cell,origin)
        del data2write
        # estimating covaraince between current map vs. static map
        _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full,bin_idx,nbin,0,nx,ny,nz)
        # get resolution from fsc
        resol = mapaverage.get_resol(f1f2_fsc,res_arr)
        ######
        # apply transformation on half maps of moving map
        '''frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+2] * st)[:,:,:,0]
        frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.unmask_fhf_lst[2*i+3] * st)[:,:,:,0]'''
        frt_hf1 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i] * st)[:,:,:,0]
        frt_hf2 = utils.get_FRS(cell,rotmat,emmap1.fhf_lst[2*i+1] * st)[:,:,:,0]
        frt_lst.append((frt_hf1+frt_hf2)/2.0)
        # find the relative B-factor between static and frt_full
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf1],cell,resol)
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],frt_hf2],cell,resol)
        #s_grid = (emmap1.s_grid**2)/4.0
        frt_hf1_scaled = scale1 * frt_hf1 * np.exp(bfac_relative1 * s_grid)
        frt_hf2_scaled = scale2 * frt_hf2 * np.exp(bfac_relative2 * s_grid)
        # re-estimate sigvar and totvar
        bin_fsc,noisevar,signalvar,totalvar,frt_full_scaled,ert_full_scaled = fsc.halfmaps_fsc_variance(
                                            frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin)
        #frt_full_scaled,ert_full_scaled,noisevar,signalvar,totalvar,bin_fsc = fcodes_fast.calc_fsc_using_halfmaps(
        #        frt_hf1_scaled,frt_hf2_scaled,bin_idx,nbin,debug_mode,nx,ny,nz) 
        maps2average_lst.append(frt_full_scaled)
        ert_lst.append(ert_full_scaled)
        # smoothen signal
        bin_fsc_full = 2*bin_fsc/(1+bin_fsc)
        '''fobj.write('\nMoving map signal \n')
        signalvar = get_extended_signal(res_arr,mapaverage.set_array(signalvar),bin_fsc_full,fobj=fobj) '''      
        signalvar_lst.append(mapaverage.set_array(signalvar))
        var_covar_mat[i,i,:] = mapaverage.set_array(signalvar)

        f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],frt_full_scaled,bin_idx,nbin,0,nx,ny,nz)
        # smoothen signal
        '''fobj.write('\nCovariance between static and moving maps \n')
        f1f2_covar = get_extended_signal(res_arr,mapaverage.set_array(f1f2_covar),f1f2_fsc,fobj=fobj)'''
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        var_covar_mat[0,i,:] = mapaverage.set_array(f1f2_covar)
        var_covar_mat[i,0,:] = mapaverage.set_array(f1f2_covar)
        #signalvar_lst.append(mapaverage.set_array(signalvar))
        totalvar_lst.append(mapaverage.set_array(totalvar))
        ######
    # calculating covariances between fullmaps
    fobj.write('\nCovariance between moving maps after fitting! \n')
    t_init = [0.0, 0.0, 0.0]
    theta_init = [(1,0,0), 0.0]
    ncycles = 10
    smax = 6
    for i in range(len(frt_lst)):
        for j in range(len(frt_lst)):
            if i == j: continue
            if i == 0 or j == 0: continue
            if j > i:
                frt_i = frt_lst[i]
                frt_j = frt_lst[j]
                # Can you do fitting between i and j here before calc. covariance
                rotmat, trans = mapaverage.run_fit(emmap1=emmap1,
                                   f1=ert_lst[i],#frt_i,
                                   f2=ert_lst[j],#frt_j,
                                   t_init=t_init,
                                   theta_init=theta_init,
                                   smax=smax,
                                   ncycles=ncycles,
                                   fobj=fobj)
                st,_,_,_ = fcodes_fast.get_st(nx,ny,nz,trans)
                frt_j_new = utils.get_FRS(cell,rotmat,frt_j * st)[:,:,:,0]
                # scale 1st map to 2nd map
                _,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                        frt_i,frt_j_new,bin_idx,nbin,0,nx,ny,nz)
                # get resolution from fsc
                resol = mapaverage.get_resol(f1f2_fsc,res_arr)
                scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_i,frt_j_new],cell,resol)
                frt_j_scaled = scale1 * frt_j_new * np.exp(bfac_relative1 * s_grid)
                # calc covariance 
                f1f2_covar,f1f2_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                        frt_i,frt_j_scaled,bin_idx,nbin,0,nx,ny,nz)
                # smoothen signal
                fobj.write('\nCovariance between map_i and map_j maps \n')
                f1f2_covar = get_extended_signal(res_arr,mapaverage.set_array(f1f2_covar),f1f2_fsc,fobj=fobj)
                # mask f1f2_covar so that anything beyond zero get zeros
                cov_lst.append(mapaverage.set_array(f1f2_covar))
                var_covar_mat[i,j,:] = mapaverage.set_array(f1f2_covar)
            #if i > j:
            #    var_covar_mat[i,j,:]
    #
    var_cov_lst = signalvar_lst + cov_lst
    '''plot_nlines_log(res_arr,
                    var_cov_lst)'''
    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple_allmaps(maps2average_lst,
                                                  cell, nbin,
                                                  var_covar_mat,
                                                  totalvar_lst,
                                                  bin_idx,
                                                  emmap1.s_grid,
                                                  res_arr,
                                                  Bf_arr)

    # calc. fsc between static map and averaged maps 
    utils.output_maps(averagemaps,emmap1.com_lst,t_lst,r_lst,emmap1.box_centr,cell,origin,Bf_arr)


def avgmaps_using_fo_nofit(emmap1, Bf_arr=[0.0], fobj=None):
    import numpy as np
    import fcodes_fast
    from emda.mapfit import utils
    from emda.plotter import plot_nlines,plot_nlines_log
    from scipy.ndimage.interpolation import shift
    from emda import iotools
    import numpy.ma as ma
    from emda import fsc

    fo_lst = emmap1.fo_lst 
    nmaps = len(fo_lst)  
    signalvar_lst = []
    totalvar_lst = []
    totalvar_lst.append(mapaverage.set_array(emmap1.totalvar_lst[0]))
    hffsc_lst = []
    hffsc_lst.append(mapaverage.set_array(emmap1.hffsc_lst[0]))
    bin_idx = emmap1.bin_idx  
    nbin = emmap1.nbin     
    res_arr = emmap1.res_arr  
    cell = emmap1.map_unit_cell 
    origin = emmap1.map_origin
    nx,ny,nz = fo_lst[0].shape
    frt_lst = []
    frt_lst.append(fo_lst[0]) 
    ert_lst = []
    ert_lst.append(emmap1.eo_lst[0])  
    cov_lst = []
    #fsc12_lst = []
    i = 0
    var_covar_mat = np.zeros((nmaps, nmaps, nbin), dtype='float')
    # static map
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(fo_lst[0]))))
    iotools.write_mrc(data2write,'static_map.mrc',cell,origin)
    del data2write
    fobj.write('\n *** Map averaging using un-normalized structure factors *** \n')
    maps2average_lst = []
    maps2average_lst.append(emmap1.fo_lst[0])
    #maps2average_lst.append(emmap1.eo_lst[0])
    signalvar_static = emmap1.signalvar_lst[0]
    signalvar_lst.append(mapaverage.set_array(signalvar_static))
    var_covar_mat[0,0,:] = mapaverage.set_array(signalvar_static)
    s_grid = (emmap1.s_grid**2)/4.0

    for fo in fo_lst[1:]:
        i = i + 1
        '''# estimating covaraince between current map vs. static map
        _,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],fo,bin_idx,nbin,0,nx,ny,nz)
        # get resolution from fsc
        resol = mapaverage.get_resol(f1f2_fsc,res_arr)
        # apply transformation on half maps of moving map
        fhf1 = emmap1.fhf_lst[2*i] 
        fhf2 = emmap1.fhf_lst[2*i+1]
        frt_lst.append((fhf1+fhf2)/2.0)
        # find the relative B-factor between static and frt_full
        scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],fhf1],cell,resol)
        scale2, bfac_relative2 = mapaverage.estimate_scale_and_relativebfac([frt_lst[0],fhf2],cell,resol)
        fhf1_scaled = scale1 * fhf1 * np.exp(bfac_relative1 * s_grid)
        fhf2_scaled = scale2 * fhf2 * np.exp(bfac_relative2 * s_grid)
        # re-estimate sigvar and totvar
        bin_fsc,noisevar,signalvar,totalvar,frt_full_scaled,ert_full_scaled = fsc.halfmaps_fsc_variance(
                                            fhf1_scaled,fhf2_scaled,bin_idx,nbin)
        #maps2average_lst.append(frt_full_scaled)
        maps2average_lst.append(ert_full_scaled)
        ert_lst.append(ert_full_scaled)
        # smoothen signal
        bin_fsc_full = 2*bin_fsc/(1+bin_fsc)      
        signalvar_lst.append(mapaverage.set_array(signalvar))
        #var_covar_mat[i,i,:] = mapaverage.set_array(signalvar)
        var_covar_mat[i,i,:] = mapaverage.set_array(emmap1.signalvar_lst[1]) # test

        #f1f2_covar,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
        #                        fo_lst[0],frt_full_scaled,bin_idx,nbin,0,nx,ny,nz)

        f1f2_covar,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],fo_lst[1],bin_idx,nbin,0,nx,ny,nz)
        # smoothen signal
        # mask f1f2_covar so that anything beyond zero get zeros
        cov_lst.append(mapaverage.set_array(f1f2_covar))
        var_covar_mat[0,i,:] = mapaverage.set_array(f1f2_covar)
        var_covar_mat[i,0,:] = mapaverage.set_array(f1f2_covar)
        #totalvar_lst.append(mapaverage.set_array(totalvar))
        totalvar_lst.append(mapaverage.set_array(emmap1.totalvar_lst[1])) # test'''
        # scaling using power spectrum
        from emda.scale_maps import scale_twomaps_by_power, transfer_power
        scale_grid = transfer_power(bin_idx=bin_idx, res_arr=res_arr, \
            scale=scale_twomaps_by_power(f1=emmap1.fo_lst[0], \
                f2=emmap1.fo_lst[i], bin_idx=bin_idx, uc=cell, \
                    res_arr=res_arr))
        fhf1 = emmap1.fhf_lst[2*i] * scale_grid
        fhf2 = emmap1.fhf_lst[2*i+1] * scale_grid
        bin_fsc,noisevar,signalvar,totalvar,fo_scaled,eo_scaled = \
            fsc.halfmaps_fsc_variance(fhf1,fhf2,bin_idx,nbin)
        #fo_scaled = emmap1.fo_lst[i] * scale_grid
        var_covar_mat[i,i,:] = mapaverage.set_array(signalvar)
        totalvar_lst.append(mapaverage.set_array(totalvar))
        maps2average_lst.append(fo_scaled)
        f1f2_covar,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                fo_lst[0],fo_scaled,bin_idx,nbin,0,nx,ny,nz)
        var_covar_mat[0,i,:] = mapaverage.set_array(f1f2_covar)
        var_covar_mat[i,0,:] = mapaverage.set_array(f1f2_covar)

        ## this is assuming already scaled.
        #maps2average_lst.append(emmap1.fo_lst[i])
        ##maps2average_lst.append(emmap1.eo_lst[i])
        #var_covar_mat[i,i,:] = mapaverage.set_array(emmap1.signalvar_lst[1])
        #f1f2_covar,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
        #                        fo_lst[0],fo_lst[1],bin_idx,nbin,0,nx,ny,nz)
        #var_covar_mat[0,i,:] = mapaverage.set_array(f1f2_covar)
        #var_covar_mat[i,0,:] = mapaverage.set_array(f1f2_covar)
        #totalvar_lst.append(mapaverage.set_array(emmap1.totalvar_lst[1]))


    '''# calculating covariances between fullmaps
    fobj.write('\nCovariance between moving maps after fitting! \n')
    for i in range(len(frt_lst)):
        for j in range(len(frt_lst)):
            if i == j: continue
            if i == 0 or j == 0: continue
            if j > i:
                frt_i = frt_lst[i]
                frt_j = frt_lst[j]
                # scale 1st map to 2nd map
                _,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                        frt_i,frt_j,bin_idx,nbin,0,nx,ny,nz)
                # get resolution from fsc
                resol = mapaverage.get_resol(f1f2_fsc,res_arr)
                scale1, bfac_relative1 = mapaverage.estimate_scale_and_relativebfac([frt_i,frt_j],cell,resol)
                frt_j_scaled = scale1 * frt_j * np.exp(bfac_relative1 * s_grid)
                # calc covariance 
                f1f2_covar,f1f2_fsc,_ = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(
                                        frt_i,frt_j_scaled,bin_idx,nbin,0,nx,ny,nz)
                # smoothen signal
                fobj.write('\nCovariance between map_i and map_j maps \n')
                f1f2_covar = get_extended_signal(res_arr,mapaverage.set_array(f1f2_covar),f1f2_fsc,fobj=fobj)
                # mask f1f2_covar so that anything beyond zero get zeros
                cov_lst.append(mapaverage.set_array(f1f2_covar))
                var_covar_mat[i,j,:] = mapaverage.set_array(f1f2_covar)'''
    # map averaging with un-normalised maps
    averagemaps = utils.calc_averagemaps_simple_allmaps(maps2average_lst,
                                                  cell, nbin,
                                                  var_covar_mat,
                                                  totalvar_lst,
                                                  bin_idx,
                                                  emmap1.s_grid,
                                                  res_arr,
                                                  Bf_arr)
    return averagemaps



def get_extended_signal(res_arr,signal,bin_fsc,fobj=None):
    from scipy import stats
    import numpy as np
    from emda.plotter import plot_nlines_log
    fobj.write('\nSignal extended using B-factor \n')
    fsc_mask = bin_fsc > 0.0
    bin_fsc_masked = bin_fsc[fsc_mask]
    #dist1 = np.sqrt((bin_fsc - 0.143)**2)
    dist2 = np.sqrt((bin_fsc - 0.3)**2)
    res_arr_masked = res_arr[fsc_mask]
    uplim = np.argmin(dist2)
    #lower_lim = np.argmin(dist1)
    '''if upper_lim < len(signal):
        lower_lim = upper_lim - 10
    else:
        lower_lim = 0'''
    # Need to determine the low and high limits
    from emda.mapfit import rdp_algo
    x = 1.0/res_arr_masked[:uplim]
    y = signal[:uplim]
    xc, yc = rdp_algo.run_rdp(x, np.log(y), epsilon=0.01)
    upl = np.argmin(np.sqrt((x - xc[-1])**2))
    lwl = np.argmin(np.sqrt((x - xc[-2])**2))
    #
    #s1 = 1.0/res_arr[signal > 0.0][60:75]
    s1 = x[lwl:upl]
    sig1 = y[lwl:upl]
    fobj.write('low resol. limit: '+ str(res_arr[signal > 0.0][lwl])+' (A) \n')
    fobj.write('high resol. limit: '+ str(res_arr[signal > 0.0][upl])+' (A) \n')
    #slope, intercept,_,_,_ = stats.linregress((s1*s1)/4.0, 
    #                                           np.log(signal[signal > 0.0][60:75]))
    slope, intercept,_,_,_ = stats.linregress((s1*s1)/4.0, np.log(sig1))
    bfac = slope
    print('B factor = ', bfac, 'Intercept = ', intercept)  
    fobj.write('B-factor: '+ str(bfac)+' \n')
    fobj.write('Intercept: '+ str(intercept)+' \n')
    s = 1.0/res_arr
    signal_pred = np.exp((bfac/4)*s**2 + intercept)
    new_signal = signal
    #new_signal[75:] = signal_pred[75:]
    new_signal[upl:] = signal_pred[upl:]
    '''plot_nlines_log(res_arr,
                    [signal, signal_pred, new_signal])'''
    return new_signal