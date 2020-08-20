"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.config import *
#debug_mode = 1

def test():
    print('restools test ... Passed')

def get_resolution_array(uc,hf1):
    import fcodes_fast
    nx,ny,nz = hf1.shape
    maxbin = np.amax(np.array([nx//2,ny//2,nz//2]))
    nbin,res_arr,bin_idx,sgrid = fcodes_fast.resolution_grid(uc,debug_mode,maxbin,nx,ny,nz)
    return nbin,res_arr[:nbin],bin_idx

def get_resArr(unit_cell,nx):
    import numpy as np
    # Generate resolution array
    a,b,c = unit_cell[:3]
    narr = np.arange(nx//2)
    narr[0] = 1.0
    fResArr = a / narr
    fResArr[0] = 5000
    #for i in range(nx//2):
    #    print(i,fResArr[i])
    return fResArr

def create_kernel(fResArr, smax):
    # Create kernel. smax is resolution to which the kernel
    # is defined.
    import numpy as np
    dist = np.sqrt((fResArr - smax)**2)
    cbin = np.argmin(dist)
    print('cbin',cbin)
    box_size = 2 * cbin + 1
    box_radius = cbin + 1
    # Box mask
    #kern_3d = np.ones((box_size, box_size, box_size)) / (box_size ** 3)
    #mask = kern_3d
    # Creating a sphere mask (binary mask)
    center = [box_radius-1, box_radius-1, box_radius-1]
    print('boxsize: ',box_size,'boxradius: ',box_radius,'center:',center)
    radius = box_radius
    X, Y, Z = np.ogrid[:box_size, :box_size, :box_size]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2 + (Z-center[2])**2)
    mask = dist_from_center <= radius
    return mask

def create_soft_edged_kernel(fResArr,smax):
    # Create soft-edged-kernel. smax is resolution to which the kernel size
    # is defined
    import numpy as np
    mask = create_kernel(fResArr,smax)
    kern_sphere = mask/np.sum(mask) # Normalized mask
    #norm = 1/np.sum(mask)
    print('Number of data points inside kernel: ', np.count_nonzero(mask))
    # Creating soft-edged mask - Naive approach
    from math import sqrt,cos,sin
    kern_sphere_soft = np.zeros(shape=(kern_sphere.shape),dtype='float')
    r1 = kern_sphere.shape[0]//2 + 1
    r0 = r1 - int(round(r1 * 0.3))
    print('r1: ', r1, 'r0: ', r0)
    #print('norm: ', norm)
    kx = kern_sphere.shape[0]
    ky = kern_sphere.shape[1]
    kz = kern_sphere.shape[2]
    center = [r1-1, r1-1, r1-1]
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center[0])**2 + (j - center[0])**2 + (k - center[0])**2)
                if dist <= r1:
                    if dist < r0:
                        #kern_sphere_soft[i,j,k] = norm
                        kern_sphere_soft[i,j,k] = 1
                    else:
                        #kern_sphere_soft[i,j,k] = (norm * (1 + cos(np.pi * (dist - r0)/(r1 - r0))))/2.0
                        kern_sphere_soft[i,j,k] = ((1 + cos(np.pi * (dist - r0)/(r1 - r0))))/2.0
    kern_sphere_soft = kern_sphere_soft / np.sum(kern_sphere_soft)
    return kern_sphere_soft

def create_soft_edged_kernel_pxl(r1):
    # Create soft-edged-kernel. r1 is the radius of kernel in pixels
    from math import sqrt,cos,sin
    import numpy as np
    #from plotter import contour_nplot2
    if r1 < 3: r1 = 3
    if r1%2 == 0: r1 = r1 - 1
    boxsize = 2 * r1 + 1
    kern_sphere_soft = np.zeros(shape=(boxsize,boxsize,boxsize),dtype='float')
    kx = ky = kz = boxsize
    center = boxsize//2
    print('center: ', center)
    r1 = center
    r0 = r1 - 2
    print('r1: ', r1, 'r0: ', r0)
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center)**2 + (j - center)**2 + (k - center)**2)
                if dist < r1:
                    if dist < r0:
                        kern_sphere_soft[i,j,k] = 1
                    else:
                        kern_sphere_soft[i,j,k] = ((1 + cos(np.pi * (dist - r0)/(r1 - r0))))/2.0
                #print(i,j,k,r1,r0,dist,kern_sphere_soft[i,j,k])
    kern_sphere_soft = kern_sphere_soft / np.sum(kern_sphere_soft)
    #contour_nplot2(kern_sphere_soft)
    return kern_sphere_soft

def remove_edge(fResArr, smax):
    # Remove everything outside smax-resolution.
    import numpy as np
    dist = np.sqrt((fResArr - smax)**2)
    cbin = np.argmin(dist)
    box_radius = cbin + 1
    box_size = cbin * 2 + 1
    # Box mask
    #kern_3d = np.ones((box_size, box_size, box_size)) / (box_size ** 3)
    #mask = kern_3d
    # Creating a sphere mask
    center = [box_radius, box_radius, box_radius]
    print('boxsize: ',box_size,'boxradius: ',box_radius,'center:',center)
    radius = box_radius
    X, Y, Z = np.ogrid[:box_size, :box_size, :box_size]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2 + (Z-center[2])**2)
    mask = dist_from_center <= radius
    return mask

def cut_resolution(f,bin_idx,res_arr,cbin):
    import fcodes_fast
    # Making data for map fitting
    nx,ny,nz = f.shape
    '''dist = np.sqrt((res_arr - smax)**2)
    cbin = np.argmin(dist) + 1 # adding 1 because fResArr starts with zero
    print('cnbin=', cbin)'''
    # fcodes.cutmap cuts f according to resolution defined by smax and output same size
    # map as f but padded with zeros outside smax
    fout = fcodes_fast.cutmap(f,bin_idx,cbin,0,len(res_arr),nx,ny,nz)
    # cutmapresize.cutmap_resize imposes smax on f and resizes mapsize
    #import cutmapresize
    #fout = cutmapresize.cutmap_resize(f,bin_idx,cbin,0,len(res_arr),nx,ny,nz)
    return fout

def cut_resolution_for_linefit(f,bin_idx,res_arr,smax):
    import fcodes_fast
    import numpy as np
    # Making data for map fitting
    nx,ny,nz = f.shape
    '''dist = np.sqrt((res_arr - smax)**2)
    cbin = np.argmin(dist) #+ 1 # adding 1 because fResArr starts with zero
    if cbin%2 != 0: cx = cbin + 1
    else: cx = cbin
    if not cbin <= bin_idx.shape[0]//2:
        print('Max allowed dim: ', bin_idx.shape[0]//2)
        print('Current dim: ', cbin)
        exit()'''
    cbin = cx = smax
    dx = int((nx - 2*cx)/2)
    dy = int((ny - 2*cx)/2)
    dz = int((nz - 2*cx)/2)
    cBIdx = bin_idx[dx:dx+2*cx,dy:dy+2*cx,dz:dz+2*cx]
    # fcodes.cutmap cuts f according to resolution defined by smax and output same size
    # map as f but padded with zeros outside smax
    fout = fcodes_fast.cutmap(f,bin_idx,cbin,0,len(res_arr),nx,ny,nz)[dx:dx+2*cx,dy:dy+2*cx,dz:dz+2*cx]
    return fout,cBIdx,cbin

