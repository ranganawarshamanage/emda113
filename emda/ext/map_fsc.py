# This module calculates FSC between maps and model
# Author: Rangana Warshamanage
# Created: 2019.06.11

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

def calc_fsc_mrc(hf1,hf2,bin_idx,nbin):
    from emda import fsc
    bin_fsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(hf1,hf2,bin_idx,nbin)
    return bin_fsc

def calculate_modelmap(uc,model,dim,resol,bfac=0.0,lig=False,lgf=None):
    import emda.emda_methods as em
    import numpy as np
    modelmap = em.model2map(modelxyz=model,
                            dim=dim,resol=resol,
                            cell=uc,lig=lig,
                            ligfile=lgf,bfac=bfac)
    f_model = np.fft.fftshift(np.fft.fftn(modelmap))
    return f_model

def pass_mtz(mtzfile,dim):
    from emda.maptools import mtz2map
    import numpy as np
    arr = mtz2map(mtzfile,dim)
    f_map = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr)))
    return f_map

def map_model_fsc(half1_map, half2_map, modelf_pdb, \
        bfac=0.0, lig=False, norm_mask=False, model1_pdb=None, \
            mask_map=None, map_size=None, model_resol=None, lgf=None):
    from emda.plotter import plot_nlines,plot_nlines2
    from emda import iotools, restools, fsc
    import emda.maskmap_class
    import numpy as np
    fsc_list = []
    ########## FSC between half maps ##########
    # if maps are in MRC format
    if half1_map.endswith(('.mrc','.mrcs','.map')):
        uc,arr1, _ = iotools.read_map(half1_map)
        uc,arr2, _ = iotools.read_map(half2_map)
    # mask taking into account
    if mask_map is not None:
         uc,msk, _ = iotools.read_map(mask_map)
    if mask_map is None:
        from emda.realsp_corr_3d import NormalizedMaps
        if norm_mask:
            nm = NormalizedMaps(hf1=np.fft.fftshift(np.fft.fftn(arr1)),
                                hf2=np.fft.fftshift(np.fft.fftn(arr2)),
                                cell=uc)
            nm.get_normdata()
            obj_maskmap = emda.maskmap_class.MaskedMaps()
            obj_maskmap.generate_mask(nm.normmap1, nm.normmap2)
            msk = obj_maskmap.mask
        else:
            # creating ccmask from half data
            obj_maskmap = emda.maskmap_class.MaskedMaps()
            obj_maskmap.generate_mask(arr1, arr2)
            msk = obj_maskmap.mask
        iotools.write_mrc(msk,'ccmask.mrc',uc)
    f_hf1 = np.fft.fftshift(np.fft.fftn(arr1 * msk))
    f_hf2 = np.fft.fftshift(np.fft.fftn(arr2 * msk))
    f_ful = (f_hf1 + f_hf2)/2.0 
    # if maps are in MTZ format
    if half1_map.endswith(('.mtz')):
        if map_size is None:
            print('Need map dimensions.')
            exit()
        dim = map_size
        if len(dim) < 3: 
            print('Need three values space delimited')
            exit()
        if len(dim) > 3:
            dim = dim[:3]
        f_hf1 = pass_mtz(half1_map, dim)
        f_hf2 = pass_mtz(half2_map, dim)  
    # making resolution grid
    dim = f_hf1.shape
    nx,ny,nz = f_hf1.shape
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f_hf1)
    # FSC between halfmaps
    bin_fsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(f_hf1,f_hf2,bin_idx,nbin)
    fsc_list.append(bin_fsc)
    #fullmap_fsc = 2.0*bin_fsc/(bin_fsc + 1.0)
    #fsc_list.append(fullmap_fsc)
    ##########################################
    if model_resol is None:
        # determine map resolution using hfmap FSC
        dist = np.sqrt((bin_fsc - 0.143)**2)
        map_resol = res_arr[np.argmin(dist)]
    else:
        map_resol = model_resol
    ########## Calculate maps from models ##########
    # if model is suppied as coordinates
    dim = [nx, ny, nz]
    if modelf_pdb.endswith(('.pdb','.ent','.cif')):
        #f_modelf = calculate_modelmap(modelf_pdb,dim,map_resol)
        f_modelf = calculate_modelmap(uc=uc,model=modelf_pdb,dim=dim,\
                                      resol=map_resol,bfac=bfac,lig=lig,lgf=lgf)
    if model1_pdb is not None:
        f_model1 = calculate_modelmap(uc=uc,model=model1_pdb,dim=dim,\
                                      resol=map_resol,bfac=bfac,lig=lig,lgf=lgf)
        ########## FSC between maps and models #########
        # FSC between halfmaps and model1
        for imap in [f_hf1,f_hf2]:
            bin_fsc = calc_fsc_mrc(imap,f_model1,bin_idx,nbin)
            fsc_list.append(bin_fsc)
    # FSC between fullmap and modelfull
    bin_fsc = calc_fsc_mrc(f_ful,f_modelf,bin_idx,nbin)
    fsc_list.append(bin_fsc)
    ################################################
    # output plots
    if len(fsc_list) == 4:
        plot_nlines(res_arr,fsc_list,
                    'allmap_fsc_modelvsmap.eps',
                    ["hf1-hf2","half1-model1","half2-model1","fullmap-model"])
        plot_nlines2(1/res_arr,fsc_list,
                    'allmap_fsc_modelvsmap-2.eps',
                    ["hf1-hf2","half1-model1","half2-model1","fullmap-model"])
    elif len(fsc_list) == 2:
        plot_nlines(res_arr,fsc_list,
                    'fsc_modelvsmap.eps',
                    ["hf1-hf2","fullmap-model"])
        plot_nlines2(1/res_arr,fsc_list,
                    'fsc_modelvsmap-2.eps',
                    ["hf1-hf2","fullmap-model"])        




    
    


