# curve fit
''' 
Insights: 
https://stackoverflow.com/questions/16745588/least-squares-minimization-complex-numbers/20104454#20104454
Method of minimization: scipy.optimize.leastsq; Accoding to scipy documentation method os Levenberg–Marquardt'''

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.config import *

def model(x, f2, s):
    import numpy as np
    k = x[0]; b = x[1]
    return k * f2 * np.exp(-b * s)

def fun(x, f1, f2, s):
    import numpy as np
    diff = (model(x, f2, s) - f1)
    z1d = np.zeros(f1.size*2, dtype = np.float64)
    z1d[0:z1d.size:2] = diff.real
    z1d[1:z1d.size:2] = diff.imag
    return z1d

def jac(x, f1, f2, s):
    import numpy as n
    J = np.zeros((f1.size*2, 2), dtype=np.float64)
    bs = x[1] * s
    ks = x[0] * s
    j0 = (np.exp(-bs) * f2)
    J[0:f1.size*2:2,0] = j0.real
    J[1:f1.size*2:2,0] = j0.imag
    j1 = (-ks * np.exp(-bs) * f2)
    J[0:f1.size*2:2,1] = j1.real
    J[1:f1.size*2:2,1] = j1.imag
    return J

def lsq(f1, f2, s, x0):
    from scipy.optimize import leastsq
    params, cov, infodict, mesg, ier = leastsq(fun, x0, Dfun=jac, args=(f1, f2, s),full_output=True)
    #print(params, mesg, ier)
    return params

def get_resolution(fhf_lst,uc):
    import fcodes_fast
    import numpy as np
    from emda import fsc, restools
    assert fhf_lst[0].shape == fhf_lst[1].shape
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,fhf_lst[0])
    bin_fsc,f1f2_covar = fsc.anytwomaps_fsc_covariance(fhf_lst[0],
                                                   fhf_lst[1],
                                                   bin_idx,nbin)
    '''f1f2_covar,bin_fsc = fcodes_fast.calc_covar_and_fsc_betwn_anytwomaps(fhf_lst[0],
                                fhf_lst[1],
                                bin_idx,
                                nbin,
                                debug_mode,
                                nx,ny,nz)'''
    bin_fsc = bin_fsc[bin_fsc > 0.1]
    if len(bin_fsc) > 0: dist = np.sqrt((bin_fsc - 0.143)**2)
    resol = res_arr[np.argmin(dist)]
    return resol

def main(fhf_lst,uc,resol=5.0):
    import fcodes_fast
    import numpy as np
    assert fhf_lst[0].shape == fhf_lst[1].shape
    nx,ny,nz = fhf_lst[0].shape
    print(nx,ny,nz)
    maxbin = np.amax(np.array([nx//2,ny//2,nz//2])) 
    _, s_grid, mask = fcodes_fast.resolution_grid_full(uc,resol,1,maxbin,nx,ny,nz)
    # optimization
    from scipy.optimize import leastsq
    x0 = np.array([1., 10.])
    f1 = (fhf_lst[0]*mask).flatten()
    f2 = (fhf_lst[1]*mask).flatten()
    s = s_grid.flatten()
    s = (s**2)/4
    params, cov, infodict, mesg, ier = leastsq(fun, x0, Dfun=jac, args=(f1, f2, s),full_output=True)
    print(params, mesg, ier)
    return params

if (__name__ == "__main__"):
    from emda.iotools import read_map
    import numpy as np
    maplist = ['/Users/ranganaw/MRC/Map_superposition/testmaps/str11.map',
            '/Users/ranganaw/MRC/Map_superposition/testmaps/emda/hf1_rotated_5deg.mrc']
            #'/Users/ranganaw/MRC/Map_superposition/testmaps/apply_bfac/avgmap__sharp50.mrc']
    fhf_lst = []
    for imap in maplist:
        uc,arr,_ = read_map(imap)
        fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr))))
    resol = get_resolution(fhf_lst,uc)
    main(fhf_lst,uc,resol)