from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.iotools import write_mrc,read_map
from emda.quaternions import get_quaternion,get_RM
from emda.mapfit.utils import double_the_axes, get_FRS
from emda.config import *

def write_mrc2(mapdata,filename,unit_cell,map_origin):
    import mrcfile as mrc
    import numpy as np
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(mapdata))))
    # removing outer regions
    nx, ny, nz = data2write.shape
    assert nx == ny == nz
    dx = int(nx/4); xx = int(nx/2)
    newdata = data2write[dx:dx+xx,dx:dx+xx,dx:dx+xx]
    file = mrc.new(name=filename, data=np.float32(newdata), compression=None, overwrite=True)
    file.header.cella.x = unit_cell[0]
    file.header.cella.y = unit_cell[1]
    file.header.cella.z = unit_cell[2]
    file.header.nxstart = map_origin[0]
    file.header.nystart = map_origin[1]
    file.header.nzstart = map_origin[2]
    file.close()

def map_transform(mapname,t,r,ax,outname='translformed.mrc'):
    import emda.iotools
    import numpy as np
    from emda.quaternions import get_quaternion,get_RM, rotationMatrixToEulerAngles
    import fcodes_fast
    from emda.mapfit import utils
    from scipy.ndimage.interpolation import shift
    uc, arr, origin = emda.iotools.read_map(mapname)
    '''theta = [ax,r]
    q = get_quaternion(theta)
    q = q/np.sqrt(np.dot(q,q))
    rotmat = get_RM(q)
    rotmat = np.transpose(rotmat)
    print('rotation matrix:', rotmat)
    print('Applied rotation in degrees [Euler angles]: ', 
           rotationMatrixToEulerAngles(rotmat) * 180./np.pi)
    print('Applied rotation in degrees [Overall]: ',
           np.arccos((np.trace(rotmat) - 1)/2) * 180./np.pi)
    print('Applied translation in Angstrom: ', t)
    hf = np.fft.fftshift(np.fft.fftn(arr))
    cx,cy,cz = hf.shape
    if np.sqrt(np.dot(t, t)) < 1.0e-3:
        st = 1.0
    else:
        t = (np.asarray(t)) / uc[:3]
        st,_,_,_ = fcodes_fast.get_st(cx,cy,cz,t)
    print('Applied translation in Angstrom: ', t * uc[:3])
    transformed_map = np.real(
                              np.fft.ifftn(
                              np.fft.ifftshift(
                              get_FRS(uc,rotmat,hf*st)[:,:,:,0])))
    '''
    # real space map interpolation - fcodes_fast
    t_frac = utils.ang2frac(t,uc[0]/arr.shape[0])
    theta = [ax,r]
    q = get_quaternion(theta)
    q = q/np.sqrt(np.dot(q,q))
    rotmat = get_RM(q)
    #rotmat = np.transpose(rotmat)
    print('rotation matrix:', rotmat)
    print('Applied rotation in degrees [Euler angles]: ', 
           rotationMatrixToEulerAngles(rotmat) * 180./np.pi)
    print('Applied rotation in degrees [Overall]: ',
           np.arccos((np.trace(rotmat) - 1)/2) * 180./np.pi)
    print('Applied translation in Angstrom: ', t)
    nx, ny,nz = arr.shape
    transformed_map = fcodes_fast.trilinear_map(rotmat,arr,nx,ny,nz)
    transformed_map = shift(transformed_map, t_frac)
    '''# using scipy - rotate function
    from scipy import ndimage
    f = np.fft.fftshift(np.fft.fftn(arr))
    frs_real = ndimage.rotate(f.real, r, axes=(1,2), reshape=False)
    frs_imag = ndimage.rotate(f.imag, r, axes=(1,2), reshape=False)
    frs = frs_real + 1j*frs_imag
    transformed_map = np.real(np.fft.ifftn(np.fft.ifftshift(frs)))
    angle_rad = np.deg2rad(r)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    rot_matrix = np.array([[c, s],
                              [-s, c]])
    print('rotation matrix:', rot_matrix)'''
    emda.iotools.write_mrc(transformed_map,outname,uc,origin)
    return transformed_map


def map_transform_using_rotmat(arr,t,rotmat):
    # rotmat output by EMDA overlay and final translation(in fractions) are
    # directly used.
    import numpy as np
    import fcodes_fast as fc
    #uc, arr, origin = emda.iotools.read_map(mapname)
    fmap = np.fft.fftshift(np.fft.fftn(arr))
    nx, ny, nz = fmap.shape
    st,_,_,_ = fc.get_st(nx,ny,nz,t)
    transformed_map = np.real(np.fft.ifftn(
                              np.fft.ifftshift(
                              get_FRS(rotmat,fmap*st,interp='cubic')[:,:,:,0])))
    return transformed_map


if (__name__ == "__main__"):
    mapname = '/Users/ranganaw/MRC/REFMAC/Bianka/fit/com_at_box_center/using_normalised_sf/shifted_to/static_map.mrc'
    uc,arr,origin = read_map(mapname)
    theta_init=[(0,0,1),180.0]
    q = get_quaternion(theta_init)
    q = q/np.sqrt(np.dot(q,q))
    rotmat = get_RM(q)
    arr = double_the_axes(arr)
    hf = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr)))
    frt_full = get_FRS(uc,rotmat,hf)[:,:,:,0]
    write_mrc2(frt_full,'staticmap_z_180deg.mrc',uc,origin)