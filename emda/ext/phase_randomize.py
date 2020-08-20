"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import argparse
from emda.iotools import read_map
from emda.plotter import plot_nlines
import emda.fsc
import fcodes_fast
from emda.restools import get_resolution_array
from emda.config import *

# Calculate FSC_true by removing mask effects using phase randomization 
# Date: 2019.07.08
# Modi: 2019.08.13
# Author: Rangana Warshamanage

cmdl_parser = argparse.ArgumentParser(
description='Computes fsc_true in resolutuon bins using phase randomisation\n')
cmdl_parser.add_argument('-h1', '--half1_map', required=True, help='Input filename for hfmap1')
cmdl_parser.add_argument('-h2', '--half2_map', required=True, help='Input filename for hfmap2')
cmdl_parser.add_argument('-m', '--mask_map', required=True, help='Input filename for mask')
cmdl_parser.add_argument('-r', '--res_rand', type=np.float32, required=False, help='Resol.(A) for phase randomization')
cmdl_parser.add_argument('-v', '--verbose', default=False,
                         help='Verbose output')

#debug_mode = 1

def randomize_phases(sf):
    import numpy as np
    random_phases = np.random.uniform(-np.pi, np.pi, sf.shape)
    amplitudes = np.sqrt(np.real(sf * sf.conj()))
    sf_phase_randomized = (amplitudes * 
                          np.cos(random_phases) + 
                          1j * amplitudes * 
                          np.sin(random_phases))
    return sf_phase_randomized

def get_randomized_sf(resol_grid,rho,resol_randomize=10.0):
    import fcodes_fast
    print('SF phase randomization')
    f_ori = np.fft.fftshift(np.fft.fftn(rho))
    f_all_random = randomize_phases(f_ori)
    nx,ny,nz = f_ori.shape
    f_beyond_random = fcodes_fast.add_random_phase_beyond(f_ori,f_all_random,
                                                          resol_grid,
                                                          resol_randomize,
                                                          nx,ny,nz)
    return f_beyond_random

def main():
    args = cmdl_parser.parse_args()
    # Read half maps
    uc,half1,origin = read_map(args.half1_map)
    uc,half2,origin = read_map(args.half2_map)
    uc,mask,origin  = read_map(args.mask_map)
    nbin,res_arr,bin_idx = get_resolution_array(uc,half1)
    # calculate fsc
    fsc_list = []
    # unmasked fsc
    umbin_fsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(np.fft.fftshift(np.fft.fftn(half1)),
                                                    np.fft.fftshift(np.fft.fftn(half2)),
                                                    bin_idx,
                                                    nbin)
    full_fsc_unmasked = 2.0 * umbin_fsc / (1.0 + umbin_fsc)
    fsc_list.append(full_fsc_unmasked)
    # masked fsc
    bin_fsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(np.fft.fftshift(np.fft.fftn(half1*mask)),
                                                  np.fft.fftshift(np.fft.fftn(half2*mask)),
                                                  bin_idx,
                                                  nbin)
    full_fsc_t = 2.0 * bin_fsc / (1.0 + bin_fsc)
    fsc_list.append(full_fsc_t)
    
    if args.res_rand is None: 
        # resol_rand is taken as the point where half-fsc_masked falls below 0.8 - RELION
        xx = (bin_fsc - 0.8)**2
        idx = xx.argmin()
        resol_rand = res_arr[idx]
    else:
        resol_rand = args.res_rand
        idx = np.argmin((res_arr - resol_rand)**2)

    # randomize phases 
    hf1_randomized = get_randomized_sf(uc,half1,resol_rand)
    hf2_randomized = get_randomized_sf(uc,half2,resol_rand)
    # get phase randomized maps
    randhalf1 = np.real(np.fft.ifftn(np.fft.ifftshift(hf1_randomized)))
    randhalf2 = np.real(np.fft.ifftn(np.fft.ifftshift(hf2_randomized)))
    rbin_fsc,_,_,_,_,_ = fsc.halfmaps_fsc_variance(np.fft.fftshift(np.fft.fftn(randhalf1*mask)),
                                                   np.fft.fftshift(np.fft.fftn(randhalf2*mask)),
                                                   bin_idx,
                                                   nbin)
    full_fsc_n = 2.0 * rbin_fsc / (1.0 + rbin_fsc)
    fsc_list.append(full_fsc_n)
    # fsc_true from Richard's formula
    fsc_true = (full_fsc_t - full_fsc_n) / (1 - full_fsc_n)
    fsc_true[:idx+2] = full_fsc_t[:idx+2] # replace fsc_true with fsc_masked_full upto resol_rand_idx + 2 (RELION uses 2)
    fsc_list.append(fsc_true)
    # writing fsc_true into a file
    f = open("fsc_true.txt", "w")
    f.write("Resol  full_fsc_t   full_fsc_true\n")
    for i in range(len(res_arr)):
        f.write("{} {} {}\n".format(res_arr[i],full_fsc_t[i],fsc_true[i]))
    f.close()
    plot_nlines(res_arr,fsc_list,'fsc.eps',["unmasked","fsc_t","fsc_n","fsc_true"])

if(__name__ == "__main__"):
    main()
