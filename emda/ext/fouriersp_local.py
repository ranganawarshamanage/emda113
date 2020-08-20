"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from timeit import default_timer as timer
import sys
from emda.iotools import read_map,write_mrc
from emda.restools import get_resArr,remove_edge,create_soft_edged_kernel_pxl
import argparse
from emda.config import *

cmdl_parser = argparse.ArgumentParser(
                                      description='Computes Fourier space correlation\n')
cmdl_parser.add_argument('-h1', '--half1_map', required=True, help='Input filename for hfmap1')
cmdl_parser.add_argument('-h2', '--half2_map', required=True, help='Input filename for hfmap2')
cmdl_parser.add_argument('-k', '--kernel_size', type=np.int32, required=False, default=5, help='Kernel size (Pixels)')
cmdl_parser.add_argument('-v', '--verbose', default=False,
                         help='Verbose output')

#debug_mode = 1 # 0: no debug info, 1: debug

def newvariance(hfmap,kern_3d):
    import scipy.signal
    convmap  = scipy.signal.fftconvolve(hfmap, kern_3d, "same")
    return convmap

def get_3dnoise_variance(hf1,hf2,kern_sphere):#,uc,origin):
    noisevar = (hf1 - hf2) / 2
    loc3_A = newvariance(noisevar, kern_sphere)
    loc3_AA = newvariance(noisevar * np.conjugate(noisevar), kern_sphere)
    noise3d = (loc3_AA - loc3_A * np.conjugate(loc3_A))
    return noise3d

def get_3dtotal_variance(hf1,hf2,kern_sphere):
    fullmap = (hf1 + hf2) / 2
    loc3_A = newvariance(fullmap, kern_sphere)
    loc3_AA = newvariance(fullmap * np.conjugate(fullmap), kern_sphere)
    totalvar_3d = loc3_AA - loc3_A * np.conjugate(loc3_A)
    return totalvar_3d

def get_3d_covariance(map_list,kern_sphere):
    # function for EMDA 3D map calculation
    loc3_A = newvariance(map_list[0], kern_sphere)
    loc3_B = newvariance(map_list[1], kern_sphere)
    loc3_AB = newvariance(map_list[0] * np.conjugate(map_list[1]), kern_sphere)
    cov3_AB = loc3_AB - loc3_A * loc3_B.conj()
    return cov3_AB

def get_3d_fouriercorrelation(hf1,hf2,kern_sphere):
    # halfmaps correlation in fourier space
    import numpy as np
    from emda.iotools import mask_by_value
    loc3_A = newvariance(hf1, kern_sphere)
    loc3_B = newvariance(hf2, kern_sphere)
    loc3_AA = newvariance(hf1 * np.conjugate(hf1), kern_sphere)
    loc3_BB = newvariance(hf2 * np.conjugate(hf2), kern_sphere)
    loc3_AB = newvariance(hf1 * np.conjugate(hf2), kern_sphere)
    cov3_AB = loc3_AB - loc3_A * loc3_B.conj()
    var3_A = loc3_AA - loc3_A * loc3_A.conj()
    var3_B = loc3_BB - loc3_B * loc3_B.conj()
    #
    '''var3_A = np.where(var3_A <= 0., 0.1, var3_A)
    var3_B = np.where(var3_B <= 0., 0.1, var3_B)
    f_halfmaps_corr = mask_by_value(cov3_AB / np.sqrt(var3_A * var3_B))
    f_halfmaps_corr = cov3_AB / np.sqrt(var3_A * var3_B)
    f_halfmaps_corr = np.where(cov3_AB <= 0., 0.,f_halfmaps_corr)'''
    #
    f_halfmaps_corr = mask_by_value(cov3_AB / np.sqrt(var3_A * var3_B))
    f_fullmap_corr = 2 * f_halfmaps_corr / (1.0 + f_halfmaps_corr)
    #f_truemap_corr = np.sqrt(f_fullmap_corr)
    return f_halfmaps_corr, f_fullmap_corr, np.sqrt(f_fullmap_corr)    

def fcc(half1_map, half2_map, kernel_size):
    from emda.iotools import read_map,write_mrc
    from emda.restools import get_resArr,remove_edge,create_soft_edged_kernel_pxl
    uc,half1,origin = read_map(half1_map)
    uc,half2,origin = read_map(half2_map)    
    nx,ny,nz = half1.shape
    fResArr = get_resArr(uc,nx)
    cut_mask = remove_edge(fResArr,fResArr[-1])
    cc_mask = np.zeros(shape=(nx,ny,nz),dtype='int')
    cx, cy, cz = cut_mask.shape
    dx = (nx - cx)//2
    dy = (ny - cy)//2
    dz = (nz - cz)//2
    cc_mask[dx:dx+cx, dy:dy+cy, dz:dz+cz] = cut_mask
    # Creating soft-edged mask
    kern_sphere_soft = create_soft_edged_kernel_pxl(kernel_size) 
    write_mrc(kern_sphere_soft,'kern_sphere_soft.mrc',uc,origin)
    f_halfmaps,f_fullmap,f_truemap = get_3d_fouriercorrelation(
                                        np.fft.fftshift(np.fft.fftn(half1)),
                                        np.fft.fftshift(np.fft.fftn(half2)),
                                        kern_sphere_soft)
    write_mrc(np.real(f_halfmaps) * cc_mask,'fouriercorr3d_halfmaps.mrc',uc,origin)
    write_mrc(np.real(f_fullmap) * cc_mask,'fouriercorr3d_fullmap.mrc',uc,origin)
    write_mrc(np.real(f_truemap) * cc_mask,'fouriercorr3d_truemap.mrc',uc,origin)


def main():
    args = cmdl_parser.parse_args()
    #
    uc,half1,origin = read_map(args.half1_map)
    uc,half2,origin = read_map(args.half2_map)
    #
    nx,ny,nz = half1.shape
    fResArr = get_resArr(uc,nx)
    cut_mask = remove_edge(fResArr,fResArr[-1])
    cc_mask = np.zeros(shape=(nx,ny,nz),dtype='int')
    cx, cy, cz = cut_mask.shape
    dx = (nx - cx)//2
    dy = (ny - cy)//2
    dz = (nz - cz)//2
    cc_mask[dx:dx+cx, dy:dy+cy, dz:dz+cz] = cut_mask
    
    # Creating soft-edged mask
    kern_sphere_soft = create_soft_edged_kernel_pxl(args.kernel_size) # sphere with radius of n pixles
    write_mrc(kern_sphere_soft,'kern_sphere_soft.mrc',uc,origin)

    f_halfmaps,f_fullmap,f_truemap = get_3d_fouriercorrelation(np.fft.fftshift(np.fft.fftn(half1)),
                              np.fft.fftshift(np.fft.fftn(half2)),
                              kern_sphere_soft)
    #write_mrc(np.real(f_halfmaps) * cc_mask,'fouriercorr3d_halfmaps.mrc',uc,origin)
    #write_mrc(np.real(f_fullmap) * cc_mask,'fouriercorr3d_fullmap.mrc',uc,origin)
    #write_mrc(np.real(f_truemap) * cc_mask,'fouriercorr3d_truemap.mrc',uc,origin)
    write_mrc(np.real(f_halfmaps),'fouriercorr3d_halfmaps.mrc',uc,origin)
    write_mrc(np.real(f_fullmap) ,'fouriercorr3d_fullmap.mrc',uc,origin)
    write_mrc(np.real(f_truemap) ,'fouriercorr3d_truemap.mrc',uc,origin)

if(__name__ == "__main__"):
    main()
