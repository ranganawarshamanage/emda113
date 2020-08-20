"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import sys
import argparse
from emda.config import *
from emda.iotools import read_map,write_mrc,read_mtz
from emda.restools import get_resArr,remove_edge,create_soft_edged_kernel_pxl
from emda.maptools import mtz2map

cmdl_parser = argparse.ArgumentParser(
                                      description='Computes 3D local correlation using maps\n')
cmdl_parser.add_argument('-h1', '--half1_map', required=True, help='Input filename for hfmap1')
cmdl_parser.add_argument('-h2', '--half2_map', required=True, help='Input filename for hfmap2')
cmdl_parser.add_argument('-ml', '--model', required=False, help='Input model MAP/MTZ file')
#cmdl_parser.add_argument('-m', '--mask_map', required=False, help='Input filename for mask')
cmdl_parser.add_argument('-k', '--kernel_size', type=np.int32, required=False, default=5, help='Kernel size (Pixels)')
cmdl_parser.add_argument('-v', '--verbose', default=False,
                         help='Verbose output')

def get_fullmapcor(hfmap_corr):
    import numpy.ma as ma
    hfmap_corr_ma = ma.masked_less_equal(hfmap_corr,0.0)
    ccf_ma = 2 * hfmap_corr_ma / (1.0 + hfmap_corr_ma)
    ccf = ccf_ma.filled(0.0)
    return ccf

def cc_twosimilarmaps(ccmap12,ccmap1,ccmap2,uc,origin):
    import numpy.ma as ma
    ccmap1_ma = ma.masked_less_equal(ccmap1,0.0)
    ccmap2_ma = ma.masked_less_equal(ccmap2,0.0)
    ccmap12_ma = ma.masked_less_equal(ccmap12,0.0)
    cc12_ma = ccmap12_ma * np.sqrt(ccmap1_ma) * np.sqrt(ccmap2_ma)
    cc12 = cc12_ma.filled(0.0)
    #write_mrc(cc12 * cc_mask,'map12_realspacecc.mrc',uc,origin)
    return cc12

def truemap_model_cc(mapmodelcc,fullmapcc):
    import numpy.ma as ma
    from emda.iotools import mask_by_value_greater
    mapmodelcc_ma = ma.masked_less_equal(mapmodelcc,0.3) # To prevent large numbers in truemapmodelcc
    fullmapcc_ma = ma.masked_less_equal(fullmapcc,0.3) # To prevent large numbers in truemapmodelcc
    #truemapmodelcc = mask_by_value_greater(mapmodelcc / np.sqrt(fullmapcc),masking_value=1.0)
    truemapmodelcc_ma = mapmodelcc_ma / np.sqrt(fullmapcc_ma)
    truemapmodelcc = mask_by_value_greater(truemapmodelcc_ma.filled(0.0),masking_value=1.0)
    #write_mrc(truemapmodelcc * cc_mask,'truemapmodel_fancg.mrc',uc,origin)
    return truemapmodelcc

def get_3d_realspcorrelation(half1,half2,kern):
    # Full map correlation using FFT convolve
    import scipy.signal
    import numpy.ma as ma
    import numpy as np
    loc3_A = scipy.signal.fftconvolve(half1, kern, "same")
    loc3_A2 = scipy.signal.fftconvolve(half1 * half1, kern, "same")
    loc3_B = scipy.signal.fftconvolve(half2, kern, "same")
    loc3_B2 = scipy.signal.fftconvolve(half2 * half2, kern, "same")
    loc3_AB = scipy.signal.fftconvolve(half1 * half2, kern, "same")

    cov3_AB = loc3_AB - loc3_A * loc3_B
    var3_A = loc3_A2 - loc3_A**2
    var3_B = loc3_B2 - loc3_B**2
    '''cov3_AB_ma = ma.masked_less_equal(cov3_AB,0.0)
    var3_A_ma = ma.masked_less_equal(var3_A,0.0)
    var3_B_ma = ma.masked_less_equal(var3_B,0.0)
    cc_realsp_ma = cov3_AB_ma / np.sqrt(var3_A_ma * var3_B_ma)
    halfmaps_cc = cc_realsp_ma.filled(0.0)'''
    cov3_AB = np.where(var3_A * var3_B <= 0., 0., cov3_AB)
    var3_A = np.where(var3_A <= 0., 0.1, var3_A)
    var3_B = np.where(var3_B <= 0., 0.1, var3_B)
    halfmaps_cc = cov3_AB / np.sqrt(var3_A * var3_B)
    halfmaps_cc = np.where(cov3_AB <= 0., 0.,halfmaps_cc)

    fullmap_cc = 2 * halfmaps_cc / (1.0 + halfmaps_cc)
    return halfmaps_cc,fullmap_cc

def get_3d_realspmapmodelcorrelation(map,model,kern_sphere):
    mapmodel_cc,_ = get_3d_realspcorrelation(map,model,kern_sphere)
    return mapmodel_cc

def calculate_modelmap(modelxyz,dim,resol,uc,bfac=0.0,lig=False,lgf=None):
    import emda.emda_methods as em
    modelmap = em.model2map(modelxyz=modelxyz,
                            dim=dim,resol=resol,
                            cell=uc,lig=lig,
                            ligfile=lgf,bfac=bfac)
    return modelmap

def rcc(half1_map, half2_map, kernel_size, norm=False, lig=False, \
        model=None, model_resol=None, mask_map=None, lgf=None):
    from emda.iotools import read_map,write_mrc
    from emda.maptools import mtz2map
    from emda import restools
    import emda.maskmap_class
    import emda.fsc
    # parameters
    hf1, hf2 = None, None
    bin_idx = None
    print('Calculating 3D correlation between half maps and fullmap. \
            Please wait...')
    uc,arr1,origin = read_map(half1_map)
    uc,arr2,origin = read_map(half2_map)
    # mask taking into account
    if mask_map is not None:
         _,cc_mask, _ = read_map(mask_map)
    if mask_map is None:
        # creating ccmask from half data
        obj_maskmap = emda.maskmap_class.MaskedMaps()
        obj_maskmap.generate_mask(arr1, arr2)
        cc_mask = obj_maskmap.mask
    #
    nx,ny,nz = arr1.shape
    # Creating soft-edged mask
    kern_sphere_soft = restools.create_soft_edged_kernel_pxl(kernel_size)
    write_mrc(kern_sphere_soft,
              'kern_sphere_soft_smax'+str(kernel_size)+'.mrc',uc,origin)
    # full map from half maps
    f_hf1 = np.fft.fftn(arr1)
    f_hf2 = np.fft.fftn(arr2)
    fullmap = np.real(np.fft.ifftn((f_hf1 + f_hf2)/2.0))
    if norm:
        hf1, hf2 = np.fft.fftshift(f_hf1), np.fft.fftshift(f_hf2)
        print('Normalising maps...')
        nm = NormalizedMaps(hf1=hf1,hf2=hf2,cell=uc)
        nm.get_normdata()
        nbin, bin_idx = nm.nbin, nm.bin_idx, 
        halffsc, res_arr = nm.binfsc, nm.res_arr
        normfull = nm.normfull
        # Real space correlation maps
        print('Calculating 3D correlation using normalized maps...')
        halfmapscc, fullmapcc = get_3d_realspcorrelation(half1=nm.normmap1, \
            half2=nm.normmap2, kern=kern_sphere_soft)
    if not norm:
        # Real space correlation maps
        print('Calculating 3D correlation...')
        halfmapscc, fullmapcc = get_3d_realspcorrelation(half1=arr1, \
            half2=arr2,kern=kern_sphere_soft)
    print('Writing out correlation maps')
    #write_mrc(mapdata=fullmapcc, \
    #    filename='rcc_fullmapxx_smax'+str(kernel_size)+'.mrc', \
    #        unit_cell=uc, map_origin=origin, label=True)
    write_mrc(mapdata=fullmapcc * cc_mask, \
        filename='rcc_fullmap_smax'+str(kernel_size)+'.mrc', \
            unit_cell=uc, map_origin=origin, label=True)
    # Map-model correlation
    if model is not None:
        print('\nMap-model correlation!\n')
        dim = [nx, ny, nz]
        if model.endswith(('.pdb','.ent','.cif')):
            if model_resol is None:
                if norm:
                    # determine map resolution using hfmap FSC
                    dist = np.sqrt((halffsc - 0.143)**2)
                    map_resol = res_arr[np.argmin(dist)]
                    print('map will be calculated upto '+ str(map_resol) + 'A')
                if not norm:
                    hf1, hf2 = np.fft.fftshift(f_hf1), np.fft.fftshift(f_hf2)           
                    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,hf1)
                    halffsc,_,_,_,_,_ = emda.fsc.halfmaps_fsc_variance(hf1, \
                        hf2,bin_idx,nbin)
                    dist = np.sqrt((halffsc - 0.143)**2)
                    map_resol = res_arr[np.argmin(dist)]
                    print('map will be calculated upto '+ str(map_resol) + 'A')
            else:
                map_resol = model_resol
                print('map will be calculated up to '+ str(map_resol) + 'A')
            print('Calculating map from model using REFMAC!')
            model_arr = calculate_modelmap(modelxyz=model,dim=dim,
                                           resol=map_resol,uc=uc,
                                           lig=lig,lgf=lgf)
        elif model.lower().endswith(('.mrcs', '.mrc', '.map')):
            _,model_arr,_ = read_map(model)
        elif model.lower().endswith('.mtz'):
            model_arr = mtz2map(model,(nx,ny,nz))
        # need to scale model_arr to fullmap
        #scaled_model = scale_model2map(fullmap, model_arr, uc)
        if norm:
            # normalisation
            normmodel = normalized(map=model_arr,bin_idx=bin_idx,nbin=nbin)
            print('Calculating model-map correlation...\n')
            mapmodelcc = get_3d_realspmapmodelcorrelation(map=normfull, \
                model=normmodel, kern_sphere=kern_sphere_soft)
        if not norm:
            print('Calculating model-map correlation...\n')
            mapmodelcc = get_3d_realspmapmodelcorrelation(map=fullmap, \
                model=model_arr, kern_sphere=kern_sphere_soft)
        #write_mrc(mapmodelcc * cc_mask,'rcc_mapmodel.mrc',uc,origin)
        write_mrc(mapdata=mapmodelcc * cc_mask, \
            filename='rcc_mapmodel_smax'+str(kernel_size)+'.mrc', \
                unit_cell=uc, map_origin=origin, label=True)
        print('Calculating truemap-model correlation...')
        # truemap-model correlation
        truemapmodelcc = truemap_model_cc(mapmodelcc,fullmapcc)
        #write_mrc(truemapmodelcc * cc_mask,'rcc_truemapmodel.mrc',uc,origin)
        write_mrc(mapdata=truemapmodelcc * cc_mask, \
            filename='rcc_truemapmodel_smax'+str(kernel_size)+'.mrc', \
                unit_cell=uc, map_origin=origin, label=True)
        print('Map-model correlation calculated!')

def scale_model2map(fullmap, model_arr, uc):
    import numpy as np
    from emda import restools
    from emda.scale_maps import scale_twomaps_by_power,transfer_power
    f1 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(fullmap)))
    f2 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(model_arr)))
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f1)   
    scale_grid = transfer_power(bin_idx,res_arr,
                                scale_twomaps_by_power(f1,
                                                       f2,
                                                       bin_idx=bin_idx,
                                                       res_arr=res_arr))
    scaled_model = np.real(np.fft.ifftn(np.fft.ifftshift(
                                                       f2 * scale_grid)))
    return scaled_model

def mapmodel_rcc(fullmap, model, resol, kernel_size=9, \
        lig=False, trim_px=1, mask_map=None, lgf=None):
    from emda.iotools import read_map,write_mrc
    from emda.maptools import mtz2map
    from emda import restools
    import emda.maskmap_class

    print('Calculating 3D correlation between half maps and fullmap. \
            Please wait...')
    uc,arr1,origin = read_map(fullmap)
    nx, ny, nz = arr1.shape
    if mask_map is not None:
        _,cc_mask, _ = read_map(mask_map)
    else:
        nbin = arr1.shape[0] // 2 - trim_px
        obj_maskmap = emda.maskmap_class.MaskedMaps()
        edge_mask = obj_maskmap.create_edgemask(nbin)
        cc_mask = np.zeros(shape=(nx,ny,nz),dtype='bool')
        cx, cy, cz = edge_mask.shape
        dx = (nx - cx)//2
        dy = (ny - cy)//2
        dz = (nz - cz)//2
        print(dx,dy,dz)
        cc_mask[dx:dx+cx, dy:dy+cy, dz:dz+cz] = edge_mask

    assert arr1.shape == cc_mask.shape
    print('\nMap-model correlation!\n')
    dim = [nx, ny, nz]
    if model.endswith(('.pdb','.ent','.cif')):
        print('map will be calculated up to '+ str(resol) + 'A')
        print('Calculating map from model using REFMAC!')
        model_arr = calculate_modelmap(modelxyz=model,dim=dim,
                                       resol=resol,uc=uc,
                                       lig=lig,lgf=lgf)
    elif model.lower().endswith(('.mrcs', '.mrc', '.map')):
        _,model_arr,_ = read_map(model)
    elif model.lower().endswith('.mtz'):
        model_arr = mtz2map(model,(nx,ny,nz))
    #write_mrc(model_arr,'modelmap.mrc',uc,origin)
    write_mrc(mapdata=model_arr, filename='modelmap.mrc', \
        unit_cell=uc, map_origin=origin)
    print('Calculating model-map correlation...\n')
    kern_sphere_soft = restools.create_soft_edged_kernel_pxl(kernel_size)
    write_mrc(kern_sphere_soft,
              'kern_sphere_soft_smax'+str(kernel_size)+'.mrc',uc,origin)
    mapmodelcc = get_3d_realspmapmodelcorrelation(arr1,model_arr, \
                                                  kern_sphere_soft)
    #write_mrc(mapmodelcc * cc_mask,'rcc_mapmodel.mrc',uc,origin)
    write_mrc(mapdata=mapmodelcc * cc_mask, filename='rcc_mapmodel.mrc', \
        unit_cell=uc, map_origin=origin, label=True)    
    print('Calculating truemap-model correlation...')
    print('Map-model correlation calculated!')

def normalized(map,bin_idx=None,nbin=None,uc=None):
    # normalise in resol bins
    from emda import fsc,restools
    if bin_idx is None:
        if uc is not None:
            nbin,res_arr,bin_idx = restools.get_resolution_array(uc,map)
    hf1 = np.fft.fftshift(np.fft.fftn(map))
    _,_,_,totalvar,_,eo = fsc.halfmaps_fsc_variance(hf1,hf1,bin_idx,nbin)
    norm_map = np.real(np.fft.ifftn(np.fft.ifftshift(eo)))
    return norm_map

def hfdata_normalized(hf1,hf2,bin_idx=None,nbin=None,uc=None):
    # normalise in resol bins
    import numpy as np
    from emda import fsc,restools
    import fcodes_fast as fc
    if bin_idx is None:
        if uc is not None:
            nbin,res_arr,bin_idx = restools.get_resolution_array(uc,hf1)
    #hf1 = np.fft.fftshift(np.fft.fftn(map))
    binfsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(hf1,hf2,bin_idx,nbin)
    _,_,_,_,_,e1 = fsc.halfmaps_fsc_variance(hf1,hf1,bin_idx,nbin)
    _,_,_,_,_,e2 = fsc.halfmaps_fsc_variance(hf2,hf2,bin_idx,nbin)
    nx,ny,nz = hf1.shape
    fsc3d = fc.read_into_grid(bin_idx,binfsc,nbin,nx,ny,nz)
    norm_map1 = np.real(np.fft.ifftn(np.fft.ifftshift(e1 * fsc3d)))
    norm_map2 = np.real(np.fft.ifftn(np.fft.ifftshift(e2 * fsc3d)))
    return norm_map1, norm_map2

class NormalizedMaps:
    def __init__(self,hf1=None,hf2=None,bin_idx=None, nbin=None, cell=None):
        self.hf1 = hf1
        self.hf2 = hf2
        self.e0 = None
        self.e1 = None
        self.e2 = None
        self.normmap1 = None
        self.normmap2 = None
        self.normfull = None
        self.bin_idx = bin_idx
        self.cell = cell
        self.nbin = nbin
        self.binfsc = None
        self.res_arr = None

    def get_parameters(self):
        from emda import restools
        if self.bin_idx is None:
            if self.cell is not None:
                self.nbin,self.res_arr,self.bin_idx = \
                    restools.get_resolution_array(self.cell,self.hf1)
    def get_normdata(self):
        import numpy as np
        from emda import fsc
        import fcodes_fast as fc
        if self.bin_idx is None:
            self.get_parameters()
        self.binfsc,_,_,_,_,self.e0 = fsc.halfmaps_fsc_variance(self.hf1, \
            self.hf2,self.bin_idx,self.nbin)
        _,_,_,_,_,self.e1 = fsc.halfmaps_fsc_variance(self.hf1,self.hf1,
                                                 self.bin_idx,self.nbin)
        _,_,_,_,_,self.e2 = fsc.halfmaps_fsc_variance(self.hf2,self.hf2,
                                                 self.bin_idx,self.nbin)
        nx,ny,nz = self.hf1.shape
        fsc3d = fc.read_into_grid(self.bin_idx,self.binfsc,self.nbin,nx,ny,nz)
        self.normmap1 = np.real(np.fft.ifftn(np.fft.ifftshift(self.e1 * fsc3d)))
        self.normmap2 = np.real(np.fft.ifftn(np.fft.ifftshift(self.e2 * fsc3d)))
        self.normfull = np.real(np.fft.ifftn(np.fft.ifftshift(self.e0 * fsc3d)))





def main():
    args = cmdl_parser.parse_args()
    # Halfmaps and fullmap correlation calculation
    print('Calculating 3D correlation between half maps and fullmap. Please wait...')
    uc,half1,origin = read_map(args.half1_map)
    uc,half2,origin = read_map(args.half2_map)
    nx,ny,nz = half1.shape
    fResArr = get_resArr(uc,nx)

    cut_mask = remove_edge(fResArr,fResArr[-1])
    cc_mask = np.zeros(shape=(nx,ny,nz),dtype='int')
    cx, cy, cz = cut_mask.shape
    dx = (nx - cx)//2
    dy = (ny - cy)//2
    dz = (nz - cz)//2
    print(dx,dy,dz)
    cc_mask[dx:dx+cx, dy:dy+cy, dz:dz+cz] = cut_mask

    # Creating soft-edged mask
    kern_sphere_soft = create_soft_edged_kernel_pxl(args.kernel_size) # sphere with radius of n pixles
    write_mrc(kern_sphere_soft,'kern_sphere_soft_smax'+str(args.kernel_size)+'.mrc',uc,origin)

    # Real space correlation maps
    halfmapscc, fullmapcc = get_3d_realspcorrelation(half1,half2,kern_sphere_soft)
    write_mrc(halfmapscc * cc_mask,'hfmaps_3dcc_smax'+str(args.kernel_size)+'.mrc',uc,origin)
    write_mrc(fullmapcc * cc_mask,'fullmap_3dcc_smax'+str(args.kernel_size)+'.mrc',uc,origin)
    
    # Create mask using correlation
    #mask = create_mask_from_rscc(fullmapcc*cc_mask,0.91)
    #write_mrc(mask * cc_mask,'mask.mrc',uc,origin)

    # Map-model correlation
    if args.model is not None:
        fullmap = (half1 + half2) / 2.0
        if args.model.lower().endswith(('.mrcs', '.mrc', '.map')):
            _,model,_ = read_map(args.model)
        elif args.model.lower().endswith('.mtz'):
            model = mtz2map(args.model,(nx,ny,nz))
        mapmodelcc = get_3d_realspmapmodelcorrelation(fullmap,model,kern_sphere_soft)
        write_mrc(mapmodelcc * cc_mask,'fullmapModel_3dcc.mrc',uc,origin)
        print('Map-model correlation calculated!')
        # truemap-model correlation
        truemapmodelcc = truemap_model_cc(mapmodelcc,fullmapcc)
        write_mrc(truemapmodelcc * cc_mask,'truemapModel_3dcc.mrc',uc,origin)

    # Two maps with and without ligand (or similar situation e.g. two conformations)
    #ccmap1 = loccc_realspace(half1,half2,kern_sphere)
    #ccmap2 = loccc_realspace(half3,half4,kern_sphere)
    #full1 = (half1 + half2) / 2.0
    #full2 = (half3 + half4) / 2.0
    #ccmap12 = loccc_realspace(full1,full2,kern_sphere)
    #cc12 = cc_twosimilarmaps(ccmap12,ccmap1,ccmap2,uc,origin)'''

if(__name__ == "__main__"):
    main()
