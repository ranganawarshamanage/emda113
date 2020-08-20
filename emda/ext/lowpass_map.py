# cut map at different resolution and output

from __future__ import absolute_import, division, print_function, unicode_literals
from emda import iotools
from emda import restools 
import argparse
import numpy as np
import fcodes_fast
from emda.config import *

cmdl_parser = argparse.ArgumentParser(
description='Lowpass filter the map to given resolution\n')
cmdl_parser.add_argument('-m', '--input_map', required=True, help='Input filename for hfmap1')
cmdl_parser.add_argument('-r', '--res_lst', nargs='+', type=np.float32, required=True, help='Resol.(A) ')
cmdl_parser.add_argument('-v', '--verbose', default=False,
                         help='Verbose output')


def lowpass_map(uc, arr, resol):
    from emda import iotools, restools
    import numpy as np
    fc = np.fft.fftshift(np.fft.fftn(arr))
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,fc)
    dist = np.sqrt((res_arr - resol)**2)
    cbin = np.argmin(dist) + 1 # adding 1 because fResArr starts with zero
    fout = restools.cut_resolution(fc,bin_idx,res_arr,cbin)
    map_lwp = np.real(np.fft.ifftn(np.fft.ifftshift(fout)))
    return fout, map_lwp  

def butterworth(uc, arr, smax, order=4):
    import numpy as np
    #import fcodes_fast
    nx,ny,nz = arr.shape
    maxbin = np.amax(np.array([nx//2,ny//2,nz//2]))
    nbin,res_arr,bin_idx,sgrid = fcodes_fast.resolution_grid(\
        uc,0,maxbin,nx,ny,nz)
    fmap = np.fft.fftshift(np.fft.fftn(arr))
    dist = np.sqrt((res_arr[:nbin] - smax)**2)
    cbin = np.argmin(dist) + 1
    cutoffres = res_arr[cbin]
    order = 4 # order of the butterworth filter
    B = 1.0; D = sgrid; d = 1./cutoffres
    bwfilter = 1./(1+B*((D/d)**(2*order))) # calculate the butterworth filter
    fmap_filtered = fmap * bwfilter 
    map_lwp = np.real(np.fft.ifftn(np.fft.ifftshift(fmap_filtered)))
    return fmap_filtered, map_lwp

def main():
    args = cmdl_parser.parse_args()
    uc,arr,origin = iotools.read_map(args.input_map)
    fc = np.fft.fftshift(np.fft.fftn(arr))
    #np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr_mask)))
    nx,ny,nz = fc.shape
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,fc)
    for resol in args.res_lst:
        #fout = restools.cut_resolution(fc,bin_idx,res_arr,resol)
        fout = restools.cut_resolution(fc,bin_idx,res_arr,resol)
        np.real(np.fft.ifftn(fout))
        iotools.write_mrc(np.real(np.fft.ifftn(np.fft.ifftshift(fout))),
                          "{0}_{1}.{2}".format('lowpass',str(resol),'mrc'),
                          uc,
                          origin)

if(__name__ == "__main__"):
    main()