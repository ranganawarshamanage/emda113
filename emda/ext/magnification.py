from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from timeit import default_timer as timer
from emda import restools
from emda import fsc as emfsc
from emda import emda_methods as em
import fcodes_fast as fcodes
#from matplotlib import pyplot as plt
from emda.restools import cut_resolution_for_linefit
from emda.config import *
 

class LineFit:
    def __init__(self):
        self.e0_lf          = None
        self.step           = None
        self.e1_lf          = None
        self.cbin_idx_lf    = None
        self.cbin_lf        = None
        self.w_grid         = None
        self.res_arr        = None
        self.mode           = 'model'


    def get_linefit_static_data(self,e0,cbin_idx,res_arr,smax):
        self.e0_lf,self.cbin_idx_lf,self.cbin_lf = cut_resolution_for_linefit(e0,
                                                    cbin_idx,
                                                    res_arr,
                                                    smax)
        self.res_arr = res_arr[:self.cbin_lf]

    def f(self, kini):
        from emda.mapfit.utils import get_fsc_wght
        #import fcmag
        nx,ny,nz = self.e0_lf.shape
        ncopies = 1
        k = self.k + kini[0] * self.step
        #print('kini[0], k: ', kini[0], k)
        eck = fcodes.tricubic_zoom(float(k), \
            self.e1_lf,0,ncopies,nx,ny,nz) 
        eckt = eck[:,:,:,0] 
        w_grid = get_fsc_wght(self.e0_lf, eckt,self.cbin_idx_lf,self.cbin_lf)
        fval = np.sum(w_grid * self.e0_lf * np.conjugate(eckt))
        #print('k: ', k, fval.real)
        return -fval.real

    def fn2(self,kini):
        #import fcmag
        nx,ny,nz = self.e0_lf.shape
        ncopies = 1
        k = self.k + kini[0] * self.step
        fcks = fcodes.tricubic_zoom(float(k), \
            self.e1_lf,0,ncopies,nx,ny,nz)
        fckt = fcks[:,:,:,0] 
        _,_,_,fval = fcodes.scalesigmafval_full( \
            self.e0_lf,fckt,self.cbin_idx_lf,self.res_arr,0, \
                self.cbin_lf,nx,ny,nz)
        #print('k, fval: ', k, fval.real)
        return fval.real

    def calc_fval_for_different_kvalues_at_this_step(self,k=1.0,e1=None):
        from scipy import optimize
        start = timer()
        nx,ny,nz = e1.shape
        cx = self.e0_lf.shape[0] // 2
        cy = cz = cx
        dx = int((nx - 2*cx)/2)
        dy = int((ny - 2*cy)/2)
        dz = int((nz - 2*cz)/2)
        self.e1_lf = e1[dx:dx+2*cx,dy:dy+2*cy,dz:dz+2*cz]
        assert self.e1_lf.shape == self.e0_lf.shape
        init_guess = [k]
        self.k = k
        if self.mode == 'map': # for map-map
            minimum = optimize.minimize(self.f, init_guess, method='Powell')
        if self.mode == 'model': # for map-model
            minimum = optimize.minimize(self.fn2, init_guess, method='Powell')
        end = timer()
        #print(' time for line search: ', end-start)
        return minimum.x



def create_xyz_grid(uc,nxyz):
    x = np.fft.fftfreq(nxyz[0]) #* uc[0]
    y = np.fft.fftfreq(nxyz[1]) #* uc[1]
    z = np.fft.fftfreq(nxyz[2]) #* uc[2]
    xv, yv, zv = np.meshgrid(x,y,z)
    xyz = [yv,xv,zv]
    for i in range(3):
        xyz[i] = np.fft.ifftshift(xyz[i])
    return xyz

def get_xyz_sum(xyz):
    xyz_sum = np.zeros(shape=(6),dtype='float')
    n = -1
    for i in range(3):  
        for j in range(3):
            if i == 0:
                sumxyz = np.sum(xyz[i] * xyz[j])
            elif i > 0 and j >= i:
                sumxyz = np.sum(xyz[i] * xyz[j])
            else:
                continue
            n = n + 1
            xyz_sum[n] = sumxyz   
    return xyz_sum

def dFks(mapin,uc):
    xyz = create_xyz_grid(uc,mapin.shape)
    vol = uc[0] * uc[1] * uc[2]

    # Calculating dFC(ks)/dk using FFT 
    dfk = np.zeros(mapin.shape, np.complex64)
    xyz_sum = 0.0
    for i in range(3):  
        xyz_sum = xyz_sum + np.sum(xyz[i])
    #dfk = np.fft.fftshift((-1/vol) * 2j * np.pi * np.fft.fftn(mapin * xyz_sum))
    dfk = np.fft.fftshift(-2j * np.pi * np.fft.fftn(mapin * xyz_sum))

    # second derivative
    xyz_sum = 0.0
    tp2 = (2.0 * np.pi)**2 
    """ for i in range(3):
        for j in range(3):
            xyz_sum = xyz_sum + np.sum(xyz[i] * xyz[j]) """
    xyz_sum = np.sum(get_xyz_sum(xyz))
    #ddfk = -(tp2/vol) * np.fft.fftshift(np.fft.fftn(mapin * xyz_sum))
    ddfk = -(tp2) * np.fft.fftshift(np.fft.fftn(mapin * xyz_sum))
    return dfk, ddfk

def dFts(Fc,sv):
    nx,ny,nz = Fc.shape
    tp2 = (2.0 * np.pi)**2
    dt_arr = np.zeros(shape=(nx,ny,nz,3),dtype='complex')
    ddt_arr = np.zeros(shape=(nx,ny,nz,3,3),dtype='complex')
    for i in range(3): 
        # 1st derivative 
        dfs = 2.0 * 1j * np.pi * sv[i] * Fc
        dt_arr[:,:,:,i] = dfs
        # 2nd derivative
        for j in range(3):
            if i == 0:
                ddfs = -tp2 * sv[i] * sv[j] * Fc
            elif i > 0 and j >= i:
                ddfs = -tp2 * sv[i] * sv[j] * Fc
            else:
                ddfs = ddt_arr[:,:,:,j,i] 
            ddt_arr[:,:,:,i,j] = ddfs
    return dt_arr, ddt_arr

def get_ll(fmaplist,bin_idx,res_arr,nbin,k,t):
    #k = 0.9
    nx,ny,nz = fmaplist[-1].shape #model
    st,s1,s2,s3 = fcodes.get_st(nx,ny,nz,t)
    sv = np.array([s1,s2,s3])
    ncopies = 1
    fcks = fcodes.tricubic_zoom(k,fmaplist[-1],0,ncopies,nx,ny,nz) 
    fckt = fcks[:,:,:,0] * st
    """ if len(fmaplist) == 3:
        fo,scale_d,sigma,noisevar,fval = fcmag.scalesigmafval( \
            fmaplist[0],fmaplist[1],fckt,bin_idx,res_arr,0,nbin,nx,ny,nz)
        return fo,fckt,scale_d,sigma,noisevar,fval,sv """
    if len(fmaplist) == 2:
        scale_d,sigma,totalvar,fval = fcodes.scalesigmafval_full( \
            fmaplist[0],fckt,bin_idx,res_arr,0,nbin,nx,ny,nz)
        return fmaplist[0],fckt,scale_d,sigma,totalvar,fval,sv

def derivatives_mapmodel(fo,fc,bin_idx,sv,D,totalvar,uc):
    # 1. Calculate dFc/dk and dFc/dT
    mapin = np.fft.ifftn(np.fft.ifftshift(fc))
    #em.write_mrc(mapin.real,"mapin.mrc",uc)
    dk, ddk = dFks(mapin,uc)
    #dt, ddt = dFts(fc,sv)

    nbin = len(totalvar)
    nx,ny,nz = fc.shape
    dll,ddll = fcodes.ll_derivatives(fo,fc,bin_idx,D,totalvar,1,nbin,nx,ny,nz)

    # 1st derivatives
    df_val = np.zeros(shape=(4), dtype='float')
    df_val[0] = np.sum(np.real(dll * np.conjugate(dk)))
    """ for i in range(3):
        df_val[i+1] = np.sum(np.real(dll * np.conjugate(2j * np.pi * sv[i] * fc))) """
    
    # 2nd derivatives
    ddf_val = np.zeros(shape=(4,4), dtype='float')
    ddf_val[0,0] = np.sum(np.real(ddll * np.conjugate(ddk)))
    """ tp2 = (2.0 * np.pi)**2
    for i in range(3):
        for j in range(3):
            ddf_val[i+1,j+1] = -tp2 * np.sum(np.real(ddll * \
                np.conjugate(fc * sv[i] * sv[j]))) """
    ddf_val_inv = np.linalg.pinv(ddf_val)
    step = ddf_val_inv.dot(-df_val)
    return step


def minimizer_mapmodel(hfmaplist,uc,ncycles=4):
    # get bin_idx
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,hfmaplist[0])
    tol = 1e-2
    if nbin <= 50:
        smax = nbin
    else:
        smax = 50
    for ifit in range(len(hfmaplist) - 1):
        f_list = []
        k_list = []
        t_list = []
        t=[0.0,0.0,0.0]; k=1.0
        for i in range(ncycles):
            start = timer()
            maplist = [hfmaplist[ifit], hfmaplist[-1]]

            fo,fckt,scale_d,sigma,totalvar,fval,sv = get_ll( \
                maplist,bin_idx,res_arr,nbin,k,t)
            #maplist = [fo,fckt]
            f_list.append(fval) 
            k_list.append(k)
            t_list.append(t)
            #print(i, fval, k)
            if i == 0:
                fval_previous = fval
                k_previous = k
                print()
                print('ifit    cycle#    func val.   magnification')
                print("{:5d} {:5d} {:8.4f} {:6.4f} {:6.4f}".format(
                    ifit, i, fval, k, k_previous))
            if i > 0 and (fval_previous - fval) > tol or i == ncycles-1:
                print("{:5d} {:5d} {:8.4f} {:6.4f} {:6.4f}".format(
                    ifit, i, fval, k, k_previous))
                k_previous = k
            if i > 0 and (fval - fval_previous) > tol or i == ncycles-1:
                #parameter output
                nx,ny,nz = hfmaplist[1].shape
                #st,_,_,_ = fcodes.get_st(nx,ny,nz,t)
                magerror = abs(k_previous - 1.0)*100.0
                print('magnification error (%): ', \
                    "{:.2f}".format(magerror))
                fck = fcodes.tricubic_zoom(float(1/k_previous), \
                    maplist[0],0,1,nx,ny,nz) 
                fckt = fck[:,:,:,0] #* st
                zoomedmap = np.fft.fftshift((np.fft.ifftn( \
                    np.fft.ifftshift(fckt))).real)
                mapname = 'magcorretedmap_'+str(ifit)+'.mrc'
                em.write_mrc(zoomedmap,mapname,uc)
                break
            step = derivatives_mapmodel(fo,fckt,bin_idx,sv,scale_d, \
                sigma,uc)
            #print(step)
            if i == 0:
                linefit = LineFit()
                linefit.get_linefit_static_data(e0=fo, \
                    cbin_idx=bin_idx,res_arr=res_arr,smax=smax)
            linefit.step = step[0]
            alpha = linefit.calc_fval_for_different_kvalues_at_this_step(\
                k=k,e1=fckt)
            k = k + (alpha) * step[0]
            t = t #+ step[1:] # t ignored as maps are analysed for 
                            # magnification at the origin of the box
            fval_previous = fval
            end = timer()
            #print('per cycle time: ', end-start)

def get_fmaplist(maplist):
    from emda.iotools import resample2staticmap
    hfmaplist = []
    # reference map
    uc_target, arr_ref, orig = em.read_map(maplist[-1])
    em.write_mrc(arr_ref,'reference.mrc',uc_target,orig)
    target_dim = list(arr_ref.shape)
    target_pix_size = uc_target[0] / target_dim[0]
    for imap in maplist[:-1]:
        uc, arr, orig = em.read_map(imap)
        curnt_pix_size = uc[0] / arr.shape[0]
        arr = resample2staticmap(curnt_pix=curnt_pix_size, \
                targt_pix=target_pix_size, targt_dim=target_dim, \
                    arr=arr)
        hfmaplist.append(np.fft.fftshift(
            np.fft.fftn(np.fft.fftshift(arr))))
    hfmaplist.append(np.fft.fftshift(
            np.fft.fftn(np.fft.fftshift(arr_ref))))
    return hfmaplist, uc_target

def main(maplist, mode='model'):
    hfmaplist, uc = get_fmaplist(maplist)
    """ if mode == 'map':
        minimizer_twomaps(hfmaplist,uc) """
    if mode == 'model':
        minimizer_mapmodel(hfmaplist,uc)