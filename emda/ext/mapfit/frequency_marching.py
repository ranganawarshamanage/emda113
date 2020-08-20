
from __future__ import absolute_import, division, print_function, unicode_literals

def frequency_marching(fo,bin_idx,res_arr,bmax=None,fobj=None):
    import numpy as np
    from emda.restools import cut_resolution

    '''if bmax is not None: 
        cbin = bmax
        smax = res_arr[cbin]
    else:
        smax = 5.0 # Angstrom
        dist = np.sqrt((res_arr - smax)**2)
        cbin = np.argmin(dist) + 1 # adding 1 because fResArr starts with zero
    if cbin%2 != 0: 
        cx = cbin + 1
        #cbin = cx
    else: cx = cbin'''
    cbin = cx = bmax
    print('cbin=', cbin)
    print('fit resolution:', res_arr[cbin])
    if fobj is not None:
        fobj.write('fit resolution: '+str(res_arr[cbin])+' (A) \n')
    print('**')
    dx = int((fo.shape[0] - 2*cx)/2)
    dy = int((fo.shape[1] - 2*cx)/2)
    dz = int((fo.shape[2] - 2*cx)/2)
    cBIdx = bin_idx[dx:dx+2*cx,dy:dy+2*cx,dz:dz+2*cx]
    cutmap = cut_resolution(fo,bin_idx,res_arr,cbin)[dx:dx+2*cx,dy:dy+2*cx,dz:dz+2*cx]
    return cutmap,cBIdx,cbin