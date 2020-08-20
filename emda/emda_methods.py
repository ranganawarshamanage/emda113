# Caller script
from __future__ import absolute_import, division, print_function, unicode_literals
from emda.config import *

def read_map(mapname):
    import emda.iotools
    uc,arr,origin = emda.iotools.read_map(mapname)
    return uc,arr,origin

def read_mtz(mtzfile):
    import emda.iotools
    uc,df = emda.iotools.read_mtz(mtzfile)
    return uc,df

def write_mrc(mapdata,filename,unit_cell,map_origin=[0.0,0.0,0.0]):
    import emda.iotools
    emda.iotools.write_mrc(mapdata,filename,unit_cell)

def write_mtz(uc,arr,outfile='output.mtz'):
    import emda.iotools
    emda.iotools.write_3d2mtz(uc,
                              arr,
                              outfile='output.mtz')

def resample_data(curnt_pix,targt_pix,targt_dim,arr):
    import emda.iotools
    new_arr = emda.iotools.resample2staticmap(curnt_pix=curnt_pix,
                targt_pix=targt_pix,targt_dim=targt_dim,
                arr=arr)
    return new_arr

def estimate_map_resol(hfmap1name,hfmap2name):
    import emda.maptools
    map_resol = emda.maptools.estimate_map_resol(hfmap1name,hfmap2name)
    return map_resol

def get_map_power(mapname):
    import emda.maptools
    res_arr,power_spectrum = emda.maptools.get_map_power(mapname)
    return res_arr,power_spectrum

def get_biso_from_model(mmcif_file):
    import emda.maptools
    biso = emda.maptools.get_biso_from_model(mmcif_file)
    return biso

def get_biso_from_map(halfmap1,halfmap2):
    import emda.maptools
    biso = emda.maptools.get_biso_from_map(halfmap1,halfmap2)
    return biso

def apply_bfactor_to_map(mapname,bf_arr,mapout):
    import emda.maptools
    all_mapout = emda.maptools.apply_bfactor_to_map(mapname,
                                                    bf_arr,
                                                    mapout)
    return all_mapout
    
def map2mtz(mapname,mtzname):
    import emda.maptools
    emda.maptools.map2mtz(mapname,mtzname=mtzname)

def map2mtzfull(uc,arr1,arr2,mtzname):
    from emda import iotools
    import numpy as np
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    iotools.write_3d2mtz_x(uc=uc,hf1=hf1,hf2=hf2,outfile=mtzname)

def mtz2map(mtzname,map_size):
    import emda.maptools
    data2write = emda.maptools.mtz2map(mtzname,map_size)
    return data2write

def lowpass_map(uc,arr1,resol, filter='ideal', order=4):
    import emda.lowpass_map
    if filter == 'ideal':
        fmap1, map1 = emda.lowpass_map.lowpass_map(uc, arr1, resol)
    if filter == 'butterworth':
        fmap1, map1 = emda.lowpass_map.butterworth(uc, arr1, resol, order)
    return fmap1, map1

def half2full(half1name,half2name,outfile):
    import emda.half2full
    uc,arr1,origin = read_map(half1name)
    uc,arr2,origin = read_map(half2name)
    fullmap = emda.half2full.half2full(arr1, arr2)
    write_mrc(fullmap,outfile,uc,origin)
    return fullmap

def map_transform(mapname,t,r,ax,outname):
    import emda.transform_map
    transformed_map = emda.transform_map.map_transform(mapname,
                                                       t,
                                                       r,
                                                       tuple(ax),
                                                       outname)
    return transformed_map

def halfmap_fsc(half1name,half2name,filename=None, maskname=None):
    import emda.restools
    import emda.fsc
    import fcodes_fast
    import numpy as np
    import os
    uc,arr1,origin = read_map(half1name)
    uc,arr2,origin = read_map(half2name)
    if maskname is not None:
        _, mask, _ = read_map(maskname)
        arr1 = arr1 * mask
        arr2 = arr2 * mask
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin,res_arr,bin_idx = emda.restools.get_resolution_array(uc,hf1)
    _,_,noisevar,signalvar,totalvar,bin_fsc,bincount = fcodes_fast.calc_fsc_using_halfmaps(hf1,
                                                    hf2,
                                                    bin_idx,
                                                    nbin,
                                                    debug_mode,
                                                    hf1.shape[0],
                                                    hf1.shape[1],
                                                    hf1.shape[2])
    if filename is not None: 
        tdata = open(filename, "w") 
        tdata.write('halfmap1 file: %s\n' %os.path.abspath(half1name))
        tdata.write('halfmap2 file: %s\n' %os.path.abspath(half2name))
        tdata.write('\n')
        tdata.write('bin # \n')  
        tdata.write('resolution (Ang.) \n') 
        tdata.write('signal variance \n')  
        tdata.write('noise variance \n')  
        tdata.write('total variance \n')  
        tdata.write('halfmap fsc \n')   
        tdata.write('# reflx \n')                                                
        for i in range(len(res_arr)):
            sv = signalvar[i]
            nv = noisevar[i]
            tv = totalvar[i]
            fsc = bin_fsc[i]
            nfc = bincount[i]
            tdata.write("{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(i, 
                                                                      res_arr[i], 
                                                                      sv, nv, tv, 
                                                                      fsc,nfc))
    print('Bin    Resolution     FSC')
    for i in range(len(res_arr)):
        print("{:5d} {:6.2f} {:14.4f}".format(i, res_arr[i], bin_fsc[i]))
    return res_arr,bin_fsc

def get_variance(half1name,half2name,filename=None, maskname=None):
    import emda.restools
    import fcodes_fast
    import numpy as np
    uc,arr1,origin = read_map(half1name)
    uc,arr2,origin = read_map(half2name)
    if maskname is not None:
        _, mask, _ = read_map(maskname)
        arr1 = arr1 * mask
        arr2 = arr2 * mask
    hf1 = np.fft.fftshift(np.fft.fftn(arr1))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin,res_arr,bin_idx = emda.restools.get_resolution_array(uc,hf1)
    _,_,noisevar,signalvar,totalvar,bin_fsc,bincount = \
        fcodes_fast.calc_fsc_using_halfmaps(hf1,
                                            hf2,
                                            bin_idx,
                                            nbin,
                                            debug_mode,
                                            hf1.shape[0],
                                            hf1.shape[1],
                                            hf1.shape[2])   
    return res_arr, noisevar, signalvar

def twomap_fsc(map1name,map2name,fobj=None,xmlobj=None):
    import emda.restools
    import emda.fsc
    import numpy as np
    import os
    uc,arr1,origin = read_map(map1name)
    uc,arr2,origin = read_map(map2name)
    f1 = np.fft.fftshift(np.fft.fftn(arr1))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin,res_arr,bin_idx = emda.restools.get_resolution_array(uc,f1)
    bin_fsc,f1f2_covar = emda.fsc.anytwomaps_fsc_covariance(f1,f2,bin_idx,nbin)
    if xmlobj is not None:
        xmlobj.map1path = os.path.abspath(map1name)
        xmlobj.map2path = os.path.abspath(map2name)
        xmlobj.res_arr = res_arr
        xmlobj.fsc = bin_fsc
        xmlobj.outmap = 'fullmap.mtz'
        xmlobj.write_xml()
    if fobj is not None: 
        fobj.write('map1 file: %s\n' %os.path.abspath(map1name))
        fobj.write('map2 file: %s\n' %os.path.abspath(map2name))
        fobj.write('\n')
        fobj.write('bin # \n')  
        fobj.write('resolution (Ang.) \n')  
        fobj.write('fsc \n')                                                   
        for bin, fsc in enumerate(bin_fsc):
            fobj.write("{:5d} {:6.2f} {:6.3f}\n".format(bin, 
                                                         res_arr[bin], 
                                                         fsc))
    print('Bin      Resolution     FSC')
    for bin, fsc in enumerate(bin_fsc):
        print("{:5d} {:6.2f} {:6.3f}".format(bin, res_arr[bin], fsc))
    return res_arr,bin_fsc

def balbes_data(map1name,map2name,fsccutoff=0.5,mode=None):
    from emda import xmlclass
    if mode=='half':
        prepare_refmac_data(hf1name=map1name, hf2name=map2name, fsccutoff=fsccutoff)
    else:
        xml = xmlclass.Xml()
        res_arr,bin_fsc = twomap_fsc(map1name=map1name, map2name=map2name, xmlobj=xml)


def singlemap_fsc(map1name):
    import emda.restools
    import emda.fsc
    import numpy as np
    from emda.fakehalf import fakehalf
    uc,arr1,origin = read_map(map1name)
    f1 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr1)))
    f2 = fakehalf(f1)
    data2write = np.real(np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(f2))))
    write_mrc(data2write,'fakehalf.mrc',uc,origin)
    nbin,res_arr,bin_idx = emda.restools.get_resolution_array(uc,f1)
    bin_fsc,noisevar,signalvar,totalvar,fo,eo = emda.fsc.halfmaps_fsc_variance(f1,
                                                            f2,bin_idx,nbin)
    print('Resolution bin     FSC')
    for i in range(len(res_arr)):
        print("{:.2f} {:.4f}".format(res_arr[i], bin_fsc[i]))
    return res_arr,bin_fsc

def mask_from_halfmaps(uc, half1, half2, radius, norm=False, iter=1, thresh=None):
    import numpy as np
    import emda.maskmap_class
    from emda.realsp_corr_3d import normalized, hfdata_normalized
    arr1, arr2 = half1, half2
    # normalized maps
    if norm:
        #arr1 = normalized(map=half1, uc=uc)
        #arr2 = normalized(map=half2, uc=uc)
        hf1 = np.fft.fftshift(np.fft.fftn(arr1))
        hf2 = np.fft.fftshift(np.fft.fftn(arr2))
        arr1,arr2 = hfdata_normalized(hf1=hf1,hf2=hf2,uc=uc)
        write_mrc(arr1,'normarr1.mrc',uc)
    obj_maskmap = emda.maskmap_class.MaskedMaps()
    obj_maskmap.generate_mask(arr1, arr2, smax=radius,iter=iter,threshold=thresh)
    mask = obj_maskmap.mask
    return mask

def mask_from_map(arr, kern, uc, prob=0.99, itr=3):
    import emda.maskmap_class as mc
    mask = mc.mapmask(arr=arr, uc=uc, kern_rad=kern, prob=prob, itr=itr)
    return mask

def sphere_kernel_softedge(radius=5):
    import emda.restools
    kernel = emda.restools.create_soft_edged_kernel_pxl(radius)
    return kernel

def overlay_maps(maplist,masklist,tra,rot,ax,ncy,res,fobj,interp,hfm):
    from emda.mapfit import mapoverlay
    theta_init = [tuple(ax),rot]
    mapoverlay.main(maplist=maplist,
                    masklist=masklist,
                    ncycles=ncy,
                    t_init=tra,
                    theta_init=theta_init,
                    smax=res,
                    fobj=fobj,
                    interp=interp,
                    halfmaps=hfm)

def average_maps(maplist,masklist,tra,rot,ax,ncy,res,fobj,interp):
    from emda.mapfit import mapaverage
    theta_init = [tuple(ax),rot]
    mapaverage.main(maplist=maplist,
                    masklist=masklist,
                    ncycles=ncy,
                    t_init=tra,
                    theta_init=theta_init,
                    smax=res,
                    fobj=fobj,
                    interp=interp)

def realsp_correlation(half1map, half2map, kernel_size, norm=False, \
        lig=False, model=None, model_resol=None, mask_map=None, lgf=None):
    from emda.realsp_corr_3d import rcc
    rcc(half1_map=half1map, half2_map=half2map, kernel_size=kernel_size, \
        norm=norm, lig=lig, model=model, model_resol=model_resol, \
        mask_map=mask_map, lgf=lgf)

def realsp_correlation_mapmodel(fullmap, model, resol, \
        kernel_size, lig=False, trimpx=1, mask_map=None, lgf=None):
    from emda.realsp_corr_3d import mapmodel_rcc as frcc
    frcc(fullmap=fullmap, model=model, kernel_size=kernel_size, \
        lig=lig, resol=resol, mask_map=mask_map, lgf=lgf, trim_px=trimpx)

def fouriersp_correlation(half1_map, half2_map, kernel_size):
    from emda.fouriersp_corr_3d import fcc
    fcc(half1_map, half2_map, kernel_size)

def map_model_validate(half1map, half2map, modelfpdb, model1pdb, mask, \
        mapsize, modelresol, bfac=0.0, lig=False, lgf=None ):
    from emda.map_fsc import map_model_fsc
    #map_model_fsc(half1map, half2map, modelfpdb, model1pdb, mask, mapsize, modelresol)
    map_model_fsc(half1_map=half1map,half2_map=half2map,modelf_pdb=modelfpdb, \
        bfac=bfac,lig=lig,model1_pdb=model1pdb,mask_map=mask, \
            map_size=mapsize, model_resol=modelresol,lgf=lgf)

def difference_map(maplist, masklist, smax, mode='norm'):
    import numpy as np
    from emda.mapfit import difference
    #for imap, imsk in zip(maplist, masklist):
    uc,arr1,origin = read_map(maplist[0])
    _,arr2,_ = read_map(maplist[1])
    _,msk1,_ = read_map(masklist[0])
    _,msk2,_ = read_map(masklist[1])
    f1 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr1)))# * msk1)))
    f2 = np.fft.fftshift(np.fft.fftn(np.fft.fftshift(arr2)))# * msk2)))
    #difference.difference_map(f1, f2, uc,origin)
    if mode=='power':
        dm1_dm2, dm2_dm1 = difference.diffmap_scalebypower(\
            f1=f1, f2=f2, cell=uc,origin=origin,smax=smax)
        # calculate map rmsd
        masked_mean = np.sum(dm1_dm2 * msk1) / np.sum(msk1)
        diff = (dm1_dm2 - masked_mean) * msk1
        rmsd = np.sqrt(np.sum(diff * diff) / np.sum(msk1))
        print('rmsd: ', rmsd)

    if mode=='norm':
        diffmap = difference.diffmap_normalisedsf(f1=f1, f2=f2, \
            cell=uc,origin=origin,smax=smax)
        list_maps = []
        for i in range(diffmap.shape[3]):
            map = np.real(np.fft.ifftshift(
                            np.fft.ifftn(
                            np.fft.ifftshift(diffmap[:,:,:,i]))))
            list_maps.append(map)

        # calculate map rmsd
        masked_mean = np.sum(list_maps[0] * msk1) / np.sum(msk1)
        diff = (list_maps[0] - masked_mean) * msk1
        rmsd = np.sqrt(np.sum(diff * diff) / np.sum(msk1))
        print('rmsd: ', rmsd)
        masked_mean = np.sum(list_maps[1] * msk2) / np.sum(msk2)
        diff = (list_maps[1] - masked_mean) * msk2
        rmsd = np.sqrt(np.sum(diff * diff) / np.sum(msk2))
        print('rmsd of diffmap_m1-m2_amp: ', rmsd)

        write_mrc(list_maps[0]*msk1, 'diffmap_m1-m2_nrm.mrc',uc, origin)
        write_mrc(list_maps[1]*msk2, 'diffmap_m2-m1_nrm.mrc',uc, origin)
        write_mrc(list_maps[2]*msk1, 'map1.mrc',uc, origin)
        write_mrc(list_maps[3]*msk2, 'map2.mrc',uc, origin)
    
    if mode=='ampli':
        diffmap = difference.diffmap_scalebyampli(f1=f1, f2=f2, \
            cell=uc,origin=origin,smax=smax)
        list_maps = []
        for i in range(diffmap.shape[3]):
            map = np.real(np.fft.ifftshift(
                            np.fft.ifftn(
                            np.fft.ifftshift(diffmap[:,:,:,i]))))
            list_maps.append(map)
        # calculate map rmsd
        masked_mean = np.sum(list_maps[0] * msk1) / np.sum(msk1)
        diff = (list_maps[0] - masked_mean) * msk1
        rmsd = np.sqrt(np.sum(diff * diff) / np.sum(msk1))
        print('rmsd of diffmap_m1-m2_amp: ', rmsd)
        masked_mean = np.sum(list_maps[1] * msk2) / np.sum(msk2)
        diff = (list_maps[1] - masked_mean) * msk2
        rmsd = np.sqrt(np.sum(diff * diff) / np.sum(msk2))
        print('rmsd of diffmap_m2-m1_amp: ', rmsd)
        # difference map output
        write_mrc(list_maps[0]*msk1, 'diffmap_m1-m2_amp.mrc',uc, origin)
        write_mrc(list_maps[1]*msk2, 'diffmap_m2-m1_amp.mrc',uc, origin)
        write_mrc(list_maps[2]*msk1, 'map1.mrc',uc, origin)
        write_mrc(list_maps[3]*msk2, 'map2.mrc',uc, origin)
        """ write_mrc(list_maps[0]*msk1, 'diffmap_m1-m2_amp.mrc',uc, origin)
        write_mrc(list_maps[1]*msk2, 'diffmap_m2-m1_amp.mrc',uc, origin)
        write_mrc(list_maps[2], 'map1.mrc',uc, origin)
        write_mrc(list_maps[3], 'map2.mrc',uc, origin) """


def applymask(mapname, maskname, outname):
    uc,arr1,origin = read_map(mapname)
    _,mask,_ = read_map(maskname)
    write_mrc(arr1*mask,outname,uc,origin)

def scale_map2map(staticmap,map2scale,outfile):
    import numpy as np
    from emda.scale_maps import scale_twomaps_by_power,transfer_power
    from emda import restools
    uc,arr1,origin = read_map(staticmap)
    uc,arr2,origin = read_map(map2scale)
    f1 = np.fft.fftshift(np.fft.fftn(arr1))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f1)
    scale_grid = transfer_power(bin_idx,res_arr,
                                scale_twomaps_by_power(f1,
                                                       f2,
                                                       bin_idx=bin_idx,
                                                       res_arr=res_arr))
    data2write = np.real(np.fft.ifftn(np.fft.ifftshift(
                                                       f2 * scale_grid)))
    write_mrc(data2write,outfile,uc,origin)
    return data2write

def bestmap(hf1name, hf2name, outfile, mode=1, knl=5, mask=None):
    import numpy as np
    from emda import restools
    from emda.mapfit import bestmap
    if mask is None: 
        msk = 1.0
    else:
        _,msk,_ = read_map(mask)
    uc,arr1,origin = read_map(hf1name)
    uc,arr2,origin = read_map(hf2name)
    if mask: print('mask is not included in FSC calculation')
    f1 = np.fft.fftshift(np.fft.fftn(arr1))# * msk))
    f2 = np.fft.fftshift(np.fft.fftn(arr2))# * msk))
    if mode == 1:
        nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f1)
        f_map = bestmap.bestmap(f1=f1, f2=f2, bin_idx=bin_idx, nbin=nbin,
                                mode=mode)
    elif mode == 2:
        f_map = bestmap.bestmap(f1=f1, f2=f2, mode=mode, kernel_size=knl)        
    data2write = np.real(np.fft.ifftn(np.fft.ifftshift(f_map))) * msk
    write_mrc(data2write,outfile,uc,origin)

def predict_fsc(hf1name, hf2name, nparticles=None, bfac=None, mask=None):
    import numpy as np
    from emda import restools,fsc,plotter
    uc,arr1,_ = read_map(hf1name)
    uc,arr2,_ = read_map(hf2name)
    if mask is not None:
        _,msk,_ = read_map(mask)
    else:
        msk = 1.0
    if nparticles is None and bfac is None:
        print('Either nparticles or bfac needs to be given!')
        exit()
    f1 = np.fft.fftshift(np.fft.fftn(arr1 * msk))
    f2 = np.fft.fftshift(np.fft.fftn(arr2 * msk))
    nbin,res_arr,bin_idx = restools.get_resolution_array(uc,f1)  
    if nparticles is not None: 
        bfac = None
        #hf1,hf2,bin_idx,nbin,res_arr,nparticles=None, bfac=None
        #np.asarray(nparticles, dtype='float')
        nparticles = 1.0 / np.asarray(nparticles, dtype='float')
        fsc_lst = fsc.predict_fsc(hf1=f1,hf2=f2,bin_idx=bin_idx,nbin=nbin, \
            nparticles=nparticles, res_arr=res_arr)
        labels = [str(i) for i in nparticles]
    if bfac is not None:
        fsc_lst = fsc.predict_fsc(hf1=f1,hf2=f2,bin_idx=bin_idx,nbin=nbin, \
            bfac=bfac, res_arr=res_arr)
        labels = [str(i) for i in bfac]
    labels.append('reference')
    plotter.plot_nlines(res_arr, fsc_lst, 
                    curve_label=labels, mapname='fsc_predicted.eps')
    return fsc_lst,res_arr,bin_idx,nbin

def prepare_refmac_data(hf1name,hf2name,outfile='fullmap.mtz', \
    bfac=None,maskname=None, xmlobj=None, fsccutoff=None):
    import numpy as np
    from emda import iotools,restools
    import fcodes_fast
    import os
    from emda import xmlclass
    xml = xmlclass.Xml()
    uc,arr1,origin = read_map(hf1name)
    uc,arr2,origin = read_map(hf2name)
    xml.map1path = os.path.abspath(hf1name)
    xml.map2path = os.path.abspath(hf2name)
    if maskname is None:
        msk = 1.0
    if maskname is not None:
        _,msk,_ = read_map(maskname)
        if arr1.shape != msk.shape:
            print('mask dim and map dim do not match!')
            print('map dim:', arr1.shape, 'mask dim:', msk.shape)
            exit()
    if bfac is None:
        bfac = np.array([0.0], dtype='float')
    else:
        bfac = np.asarray(bfac, dtype='float')
        bfac = np.insert(bfac,0,0.0)
    hf1 = np.fft.fftshift(np.fft.fftn(arr1 * msk))
    hf2 = np.fft.fftshift(np.fft.fftn(arr2 * msk))
    # stats from half maps
    nx,ny,nz = hf1.shape
    maxbin = np.amax(np.array([nx//2,ny//2,nz//2]))
    nbin,res_arr,bin_idx,sgrid = fcodes_fast.resolution_grid(uc,debug_mode,maxbin,nx,ny,nz)
    res_arr = res_arr[:nbin]
    _,_,nvar,svar,tvar,binfsc,bincount = fcodes_fast.calc_fsc_using_halfmaps(hf1,
                                                    hf2,
                                                    bin_idx,
                                                    nbin,
                                                    debug_mode,
                                                    hf1.shape[0],
                                                    hf1.shape[1],
                                                    hf1.shape[2])
    xml.res_arr = res_arr; xml.fsc = binfsc
    xml.outmap = outfile
    # determine resolution
    bin_fsc = binfsc[binfsc > 0.1]
    if len(bin_fsc) > 0:
        if fsccutoff is None:
            fsccutoff = 0.5
        dist500 = np.sqrt((bin_fsc - float(fsccutoff))**2)
        dist143 = np.sqrt((bin_fsc - 0.143)**2)
        xml.fsccutoff1 = float(0.143)
        xml.fsccutoff2 = float(fsccutoff)
        xml.mapres1 = res_arr[np.argmin(dist143)]
        xml.mapres2 = res_arr[np.argmin(dist500)]
    xml.write_xml()

    tdata = open('table_variances.txt', "w") 
    tdata.write('halfmap1 file: %s\n' %os.path.abspath(hf1name))
    tdata.write('halfmap2 file: %s\n' %os.path.abspath(hf2name))
    tdata.write('\n')
    tdata.write('bin # \n')  
    tdata.write('resolution (Ang.) \n') 
    tdata.write('signal variance \n')  
    tdata.write('noise variance \n')  
    tdata.write('total variance \n')  
    tdata.write('halfmap fsc \n')   
    tdata.write('# reflx \n')                                                
    for i in range(len(res_arr)):
        sv = svar[i]
        nv = nvar[i]
        tv = tvar[i]
        fsc = binfsc[i]
        nfc = bincount[i]
        tdata.write("{:-3d} {:-6.2f} {:-14.4f} {:-14.4f} {:-14.4f} {:-14.4f} {:-10d}\n".format(i, 
                                                                  res_arr[i], 
                                                                  sv, nv, tv, 
                                                                  fsc,nfc))
    print('Bin    Resolution     FSC')
    for i in range(len(res_arr)):
        print("{:5d} {:6.2f} {:14.4f}".format(i, res_arr[i], binfsc[i]))
    # output mtz file    
    iotools.write_3d2mtz_refmac(uc,
                                sgrid,
                                (hf1+hf2)/2.0,
                                (hf1-hf2),
                                bfac,
                                outfile=outfile)
    
def overall_cc(map1name, map2name, space='real', maskname=None):
    from emda import cc
    import numpy as np
    uc, arr1, origin = read_map(map1name)
    uc, arr2, origin = read_map(map2name)
    if maskname is not None:
        uc, msk, origin = read_map(maskname)
        arr1 = arr1 * msk
        arr2 = arr2 * msk
    if space == 'fourier':
        print('Overall CC calculation in Fourier space')
        f1 = np.fft.fftn(arr1)
        f2 = np.fft.fftn(arr2)
        occ = cc.cc_overall_fouriersp(f1=f1, f2=f2)
        print('Overall Correlation in Fourier space= ', occ)
    else:
        print('Overall CC calculation in Real/Image space')
        occ = cc.cc_overall_realsp(map1=arr1,map2=arr2)
        print('Overall Correlation in real space= ', occ)
    return occ

def mirror_map(mapname):
    # gives the inverted copy of the map
    import numpy as np
    uc, arr, origin = read_map(mapname)
    data = np.real(np.fft.ifftn(np.conjugate(np.fft.fftn(arr))))
    write_mrc(data,'mirror.mrc',uc,origin)

def model2map(modelxyz,dim,resol,cell,bfac=0.0,lig=False,ligfile=None):
    import numpy as np
    from emda.iotools import run_refmac_sfcalc
    from emda.maptools import mtz2map
    import gemmi as gm
    # check for valid sampling: 
    if np.any(np.mod(dim, 2)) != 0:
        dim = dim + 1
    # check for minimum sampling
    min_pix_size = resol / 2 # in Angstrom
    min_dim = np.asarray(cell[:3], dtype='float') / min_pix_size
    min_dim = np.ceil(min_dim).astype(int)
    if np.any(np.mod(min_dim, 2)) != 0:
        min_dim = min_dim + 1
    if min_dim[0] > dim[0]:
        print('Minimum dim should be: ', min_dim)
        exit()
    # replace/add cell and write model.cif
    a, b, c = cell[:3]
    structure = gm.read_structure(modelxyz)
    structure.cell.set(a, b, c, 90.0, 90.0, 90.0)
    structure.make_mmcif_document().write_file('model.cif')
    # run refmac using model.cif just created
    run_refmac_sfcalc('./model.cif',resol,bfac,lig=lig,ligfile=ligfile)
    modelmap = mtz2map('./sfcalc_from_crd.mtz',dim)
    return modelmap

def read_atomsf(atom, fpath=None):
    from emda import iotools
    import subprocess
    if fpath is None:
        CMD = 'echo $CLIBD'
        p = subprocess.Popen(CMD, stdout=subprocess.PIPE, shell=True, \
            executable='/bin/bash')
        list_of_strings = [x.decode('utf-8').rstrip('\n') \
            for x in iter(p.stdout.readlines())]
        fpath = list_of_strings[0]+'/atomsf.lib'
    z, ne,a, b, ier = iotools.read_atomsf(atom, fpath=fpath)
    return z, ne, a, b, ier

def compositemap(maps, masks):
    from emda import composite
    composite.main(mapslist=maps, masklist=masks)

def mapmagnification(maplist, rmap):
    from emda import magnification
    maplist.append(rmap)
    magnification.main(maplist=maplist)

    

