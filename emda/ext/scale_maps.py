from emda.config import *

# calculate scale between two maps in resolution bins
def scale_onemap2another(f1, f2, bin_idx, res_arr):
    import numpy as np
    import fcodes_fast
    assert f1.shape == f2.shape == bin_idx.shape
    nx, ny, nz = bin_idx.shape
    scale_np = np.zeros(len(res_arr), dtype='float')
    s_grid = fcodes_fast.read_into_grid(bin_idx,
                                1.0/res_arr,
                                len(res_arr),
                                nx,ny,nz)
    print('resolution       scale')
    for i, res in enumerate(res_arr):
        slope = estimate_scale(f1[bin_idx==i], 
                                f2[bin_idx==i], 
                                s_grid[bin_idx==i])
        #scale_np[i] = params[0]
        scale_np[i] = slope
        print("{:8.4f} {:6.2f}".format(res, slope))   
    return scale_np

def estimate_scale(f1, f2, s):
    import numpy as np
    from emda.mapfit import curve_fit_3
    from scipy import stats
    x0 = np.array([1., 10.])
    if f1.ndim > 1:
        f1 = f1.flatten()
    if f2.ndim > 1:
        f1 = f2.flatten()
    if s.ndim > 1:
        f1 = s.flatten()            
    s = (s**2)/4
    #params = curve_fit_3.lsq(f1, f2, s, x0)
    # just scaling
    slope, intercept,_,_,_ = stats.linregress(np.real(f1*f1.conj()), 
                                              np.real(f2*f2.conj()))
    return slope

def scale_twomaps_by_power(f1, f2, bin_idx=None, uc=None, res_arr=None):
    from emda import restools, fsc, plotter
    import fcodes_fast
    import numpy as np
    from emda.mapfit import mapaverage
    assert f1.shape == f2.shape == bin_idx.shape
    nx, ny, nz = f1.shape
    if bin_idx is None:
        nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f1)
    else:
        nbin = np.max(bin_idx) + 1
    # find how far the signal is
    f1f2_fsc,f1f2_covar = fsc.anytwomaps_fsc_covariance(f1=f1,
                                                        f2=f2,
                                                        bin_idx=bin_idx,
                                                        nbin=nbin)
    mask = (mapaverage.set_array(f1f2_covar, 0.1) > 0.0).astype(int)
    inverse_mask = (mask < 1).astype(int)
    power_1 = fcodes_fast.calc_power_spectrum(f1,bin_idx,nbin,debug_mode,nx,ny,nz)
    power_2 = fcodes_fast.calc_power_spectrum(f2,bin_idx,nbin,debug_mode,nx,ny,nz)
    scale_np = power_2/power_1
    scale_np = scale_np * mask + inverse_mask.astype(float)
    scale_np = 1.0 / scale_np
    plotter.plot_nlines_log(res_arr,
                        [power_1,power_2,power_2*scale_np],
                        ["power1","power2","scaled2to1"],
                        'log_totalvariances.eps')
    for i, res in enumerate(res_arr):
        print("{:8.4f} {:6.2f}".format(res, scale_np[i]))
    return scale_np

def transfer_power(bin_idx,res_arr,scale):
    import fcodes_fast
    nx, ny, nz = bin_idx.shape
    scale_grid = fcodes_fast.read_into_grid(bin_idx,
                                scale,
                                len(res_arr),
                                nx,ny,nz)
    return scale_grid

if (__name__ == "__main__"):
    from emda.iotools import read_map
    import numpy as np
    from emda.restools import get_resolution_array
    maplist = ['/Users/ranganaw/MRC/REFMAC/Bianka/emda_overlay_test/emd_4572_rot180_axr010.mrc',
                #'/Users/ranganaw/MRC/REFMAC/Bianka/emda_overlay_test/emd_4572_rot180_axr010.mrc']
                #'/Users/ranganaw/MRC/REFMAC/Bianka/emda_overlay_test/emd_4572_rot180_axr010_tra200.mrc']
            '/Users/ranganaw/MRC/REFMAC/Bianka/emda_overlay_test/emd_4572_rot180_axr010_rotated20_axr100.mrc']
        
    fhf_lst = []
    for imap in maplist:
        uc,arr,_ = read_map(imap)
        fhf_lst.append(np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr))))
    nbin,res_arr,bin_idx = get_resolution_array(uc,fhf_lst[0])
    '''scale_np = scale_onemap2another(fhf_lst[0],
                         fhf_lst[1],
                         bin_idx,
                         res_arr)'''
    scale_np = scale_twomaps_by_power(fhf_lst[0],
                                      fhf_lst[1],
                                      bin_idx=bin_idx,
                                      res_arr=res_arr)
