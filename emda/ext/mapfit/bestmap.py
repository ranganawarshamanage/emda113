# calcuating bestmap - normalized and FSC weighted map
def bestmap(f1, f2, mode=1, kernel_size=5, bin_idx=None, nbin=None):
    if mode==1:
        print('bestmap is calculated using FSC in resolution bins')
        bestmap = bestmap_1dfsc(f1, f2, bin_idx, nbin)
    elif mode==2:
        print('bestmap is calculated using local FSC')
        bestmap = bestmap_3dfsc(f1, f2, kernel_size)
    return bestmap

def bestmap_1dfsc(f1, f2, bin_idx, nbin):
    import fcodes_fast
    from emda import fsc
    import numpy as np
    cx,cy,cz = f1.shape
    fsc,_,_,_,_,eo = fsc.halfmaps_fsc_variance(f1, f2, bin_idx, nbin)
    fsc_grid = fcodes_fast.read_into_grid(bin_idx,
                            fsc,
                            nbin,
                            cx,cy,cz)
    fsc_grid_filtered = np.where(fsc_grid < 0.0, 0.0, fsc_grid)
    return np.sqrt(fsc_grid_filtered) * eo

def bestmap_3dfsc(f1, f2, kernel_size=5):
    import numpy as np
    from emda.restools import get_resArr,remove_edge,create_soft_edged_kernel_pxl
    from emda.fouriersp_corr_3d import get_3dnoise_variance,get_3dtotal_variance
    kernel = create_soft_edged_kernel_pxl(kernel_size)
    # calculate 3D fourier correlation
    noise = get_3dnoise_variance(f1, f2, kernel)
    total = get_3dtotal_variance(f1, f2, kernel)
    total_pos = np.where(total <= 0., 0.1, total)
    fsc = (total - noise) / total_pos
    fsc_grid_filtered = np.where((total - noise) <= 0., 0., fsc)
    eo = (f1 + f2) / (2 * np.sqrt(total_pos))
    eo = np.where(total <= 0., 0., eo)
    return np.sqrt(fsc_grid_filtered) * eo