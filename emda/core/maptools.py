"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.config import *

def test():
    print('maptools test ... Passed')

def Bf_linefit(res_arr,
               binSgnlVar,
               low_res_cutoff,
               high_res_cutoff,
               scale=1):
    from scipy import stats
    import numpy as np
    import numpy.ma as ma
    prediSigVar = np.zeros(len(res_arr), dtype='float')
    dist_low = np.sqrt((res_arr - low_res_cutoff)**2)
    lb = np.argmin(dist_low) + 1
    dist_high = np.sqrt((res_arr - high_res_cutoff)**2)
    ub = np.argmin(dist_high) + 1
    res_arr_ma = ma.masked_equal(res_arr,0.0)
    s = 1/res_arr_ma
    x_full = (s * s) / 2.0
    x_full = x_full.filled(0.0)
    slope, intercept, rv, pv, se = stats.linregress(x_full[lb:ub],
                                                    np.log(binSgnlVar[lb:ub]))
    Bfac = slope
    #print('B factor = ', Bfac, 'Intercept = ', intercept)
    ln_s = -1.0*abs(Bfac) * x_full + intercept
    prediSigVar = scale * np.exp(ln_s)
    return Bfac,prediSigVar

def estimate_map_resol(hfmap1,hfmap2):
    from emda.iotools import read_map
    from emda import restools
    import emda.fsc
    import numpy as np
    uc,arr1,origin = read_map(hfmap1)
    uc,arr2,origin = read_map(hfmap2)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc, hf1)
    bin_fsc,_,_,_,_,_ = emda.fsc.halfmaps_fsc_variance(hf1,hf2,bin_idx,nbin)
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    dist = np.sqrt((bin_fsc - 0.143)**2)
    map_resol = res_arr[np.argmin(dist)]
    return map_resol

def get_map_power(mapin):
    from emda.iotools import read_map
    from emda import restools
    import fcodes_fast
    import numpy as np
    uc,arr,_ = read_map(mapin)
    hf = np.fft.fftshift(np.fft.fftn(arr))
    nx, ny, nz = hf.shape
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc, hf)
    power_spectrum = fcodes_fast.calc_power_spectrum(hf,bin_idx,nbin,debug_mode,nx,ny,nz)
    print('Resolution   bin     Power')
    for i in range(len(res_arr)):
        print("{:.2f} {:.4f}".format(res_arr[i], power_spectrum[i]))
    return res_arr,power_spectrum

def get_biso_from_model(mmcif_file):
    from emda.iotools import read_mmcif
    _,_,_,_,Biso_np = read_mmcif(mmcif_file)
    Biso = np.median(Biso_np)
    return Biso

def get_biso_from_map(halfmap1,halfmap2):
    from emda.plotter import plot_nlines_log
    import numpy as np
    map_resol,sig_var,res_arr = estimate_map_resol(halfmap1,halfmap2)
    plot_nlines_log(res_arr,[sig_var],["Signal Variance"])
    low_res_cutoff = np.float32(input('Enter low resolution cutoff: '))
    high_res_cutoff = np.float32(input('Enter high resolution cutoff: '))
    Biso,preditcted_signal = Bf_linefit(res_arr,sig_var,low_res_cutoff,high_res_cutoff)
    plot_nlines_log(res_arr,[sig_var,preditcted_signal],["Signal Variance","Predicted SV"],'Predicted.eps')
    return Biso

def apply_bfactor_to_map(mapname,bf_arr,mapout):
    from emda import iotools
    import fcodes_fast
    import numpy as np
    uc,ar1,origin = iotools.read_map(mapname)
    hf1 = np.fft.fftshift(np.fft.fftn(ar1))
    nx, ny, nz = ar1.shape
    nbf = len(bf_arr)
    all_mapout = fcodes_fast.apply_bfactor_to_map(hf1,bf_arr,uc,debug_mode,nx,ny,nz,nbf)
    if mapout:
        for i in range(all_mapout.shape[3]):
            if bf_arr[i] < 0.0:
                Bcode = '_blur'+str(abs(bf_arr[i]))
            elif bf_arr[i] > 0.0:
                Bcode = '_sharp'+str(abs(bf_arr[i]))
            elif bf_arr[i] == 0.0:
                Bcode = '_unsharpened'
            filename_mrc = mapname[:-4]+Bcode+'.mrc'
            data2write = np.real(np.fft.ifftn(np.fft.ifftshift(all_mapout[:,:,:,i])))
            iotools.write_mrc(data2write,filename_mrc,uc,origin)
            print('writing done')
    return all_mapout

def map2mtz(mapname,mtzname='map2mtz.mtz',factor=1.0):
    from emda.iotools import read_map,write_3d2mtz
    import numpy as np
    uc,ar1,_ = read_map(mapname)
    uc[0] = uc[0] * factor
    uc[1] = uc[1] * factor
    uc[2] = uc[2] * factor
    hf1 = np.fft.fftshift(np.fft.fftn(ar1))
    write_3d2mtz(uc,hf1,outfile=mtzname)

def mtz2map(mtzname,map_size):
    import numpy as np
    from numpy.fft import ifftn,ifftshift
    from emda.iotools import read_mtz
    import fcodes_fast
    _,dataframe = read_mtz(mtzname)
    h = dataframe['H'].astype('int')
    k = dataframe['K'].astype('int')
    l = dataframe['L'].astype('int')
    f = dataframe['Fout0'] * np.exp(np.pi * 1j * dataframe['Pout0']/180.0)
    nx, ny, nz = map_size
    f3d = fcodes_fast.mtz2_3d(h,k,l,f,nx,ny,nz,len(f))
    data2write = np.real((ifftn(ifftshift(f3d))))
    return data2write




