from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from emda.quaternions import derivatives_wrt_q
from timeit import default_timer as timer

debug_mode = 0
timeit = False

def dFs2(mapin,xyz,vol):
    import numpy as np
    nx,ny,nz = mapin.shape
    #start = timer()
    #tp2 = (2.0 * np.pi)**2 
    #n = -1
    dfs_arr = np.zeros(shape=(nx,ny,nz,3),dtype=np.complex64)
    '''ddfs_arr = np.zeros(shape=(nx,ny,nz,6),dtype=np.complex64)'''
    # Calculating dFC/ds using FFT - 1st derivative using eq. 15
    for i in range(3):  
        dfs = (1/vol) * 2j * np.pi * np.fft.fftn(mapin * xyz[i])
        dfs_arr[:,:,:,i] = np.fft.fftshift(dfs)
        # Calculating ddFC/ds using FFT - 2nd derivative using eq. 17
        '''for j in range(3):
            if i == 0:
                ddfs = -(tp2/vol) * np.fft.fftn(mapin * xyz[i] * xyz[j])
            elif i > 0 and j >= i:
                ddfs = -(tp2/vol) * np.fft.fftn(mapin * xyz[i] * xyz[j])
            else:
                continue
            n = n + 1
            ddfs_arr[:,:,:,n] = np.fft.fftshift(ddfs) '''   
    #end = timer() 
    #print('time for dfs and ddfs calculation: ', end-start)
    return dfs_arr#,ddfs_arr

def derivatives(e0,e1,w_grid,w2_grid,q,sv,xyz,xyz_sum,vol):
    import numpy as np
    import sys
    df_val = np.zeros(6,dtype=np.float64)
    ddf_val = np.zeros(shape=(6,6),dtype=np.float64)
    tp2 = (2.0 * np.pi)**2
    # Translation derivatives
    start = timer()
    for i in range(3):
        # 1st derivative. eq. 33
        df_tmp = np.sum(w_grid * e0 * 
                np.conjugate(e1 * (2.0 * np.pi * 1j) * sv[i])) 
        df_val[i] = np.real(df_tmp)
        for j in range(3):
            # estimated 2nd derivative
            if i == 0:
                ddf_tmp = -tp2 * np.sum(w2_grid * sv[i] * sv[j]) 
            elif i > 0 and j >= i:
                ddf_tmp = -tp2 * np.sum(w2_grid * sv[i] * sv[j]) 
            else:
                ddf_val[i,j] = ddf_val[j,i]
            ddf_val[i,j] = np.real(ddf_tmp)
    
    # Rotation derivatives
    dRdq = derivatives_wrt_q(q) # give current q - vector
    # impose current rotation and translation on F2
    start_dfs2 = timer()
    dFRs = dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))),xyz,vol)
    end_dfs2 = timer()
    print('     time for DFS calc: ', end_dfs2-start_dfs2)
    a = np.zeros(shape=(3,3),dtype=np.float64)
    b = np.zeros(shape=(3,3),dtype=np.float64)
    start_dfr = timer()
    for i in range(3):
        a[:,:] = 0.0
        for k in range(3):
            for l in range(3):
                if k == 0:
                    tmp1 = np.sum(w_grid * np.conjugate(e0) * 
                    (dFRs[:,:,:,k] * sv[l] * dRdq[i,k,l]))
                elif k > 0 and l >= k:
                    tmp1 = np.sum(w_grid * np.conjugate(e0) * 
                    (dFRs[:,:,:,k] * sv[l] * dRdq[i,k,l]))
                else:
                    a[k,l] = a[l,k]
                a[k,l] = tmp1.real
        df_val[i+3] = np.sum(a) # df_val[3] to df_val[5]
    end_dfr = timer()
    print('     time for dfr: ', end_dfr-start_dfr)
    wfsc = w_grid * np.conjugate(e0) * e1 # THIS WORKS
    start_ddfr = timer()
    for i in range(3):
        for j in range(3):
            if i == 0:
                b[:,:] = 0.0
                n = -1
                for k in range(3):
                    for l in range(3): 
                        if k == 0:
                            n = n + 1
                            #tmp2 = np.sum(w_grid * np.conjugate(e0) * 
                            #    (ddfrs_all[:,:,:,n] * sv[k] * sv[l] * 
                            #    dRdq[i,k,l] * dRdq[j,k,l]))
                            tmp2 = -(tp2/vol) * xyz_sum[n] * np.sum(wfsc * sv[k] * sv[l] *
                                dRdq[i,k,l] * dRdq[j,k,l])
                        elif k > 0 and l >= k:
                            n = n + 1
                            #tmp2 = np.sum(w_grid * np.conjugate(e0) * 
                            #    (ddfrs_all[:,:,:,n] * sv[k] * sv[l] * 
                            #    dRdq[i,k,l] * dRdq[j,k,l]))
                            tmp2 = -(tp2/vol) * xyz_sum[n] * np.sum(wfsc * sv[k] * sv[l] *
                                dRdq[i,k,l] * dRdq[j,k,l])
                        else:
                            b[k,l] = b[l,k]
                        b[k,l] = tmp2.real
                ddf_val[i+3,j+3] = np.sum(b) # ddf_val[3] to [5] in i and j indices
            elif i > 0 and j >= i:
                b[:,:] = 0.0
                n = -1
                for k in range(3):
                    for l in range(3): 
                        if k == 0:
                            n = n + 1
                            #tmp2 = np.sum(w_grid * np.conjugate(e0) * 
                            #    (ddfrs_all[:,:,:,n] * sv[k] * sv[l] * 
                            #    dRdq[i,k,l] * dRdq[j,k,l]))
                            tmp2 = -(tp2/vol) * xyz_sum[n] * np.sum(wfsc * sv[k] * sv[l] *
                                dRdq[i,k,l] * dRdq[j,k,l])
                        elif k > 0 and l >= k:
                            n = n + 1
                            #tmp2 = np.sum(w_grid * np.conjugate(e0) * 
                            #    (ddfrs_all[:,:,:,n] * sv[k] * sv[l] * 
                            #    dRdq[i,k,l] * dRdq[j,k,l]))
                            tmp2 = -(tp2/vol) * xyz_sum[n] * np.sum(wfsc * sv[k] * sv[l] *
                                dRdq[i,k,l] * dRdq[j,k,l])
                        else:
                            b[k,l] = b[l,k]
                        b[k,l] = tmp2.real
                ddf_val[i+3,j+3] = np.sum(b) # ddf_val[3] to [5] in i and j indices
            else:
                ddf_val[i+3,j+3] = ddf_val[j+3,i+3]    
    #print(ddf_val)
    # Mixed derivatives
    # NEED A REVIEW
    #for i in range(3):
    #    for j in range(3):
    #        a[:,:] = 0.0
    #        for k in range(3):
    #            for l in range(3):
    #                tmp1 = np.sum(w_grid * np.conjugate(e0) *
    #                    (dFRs[:,:,:,k] * sv[l] * sv[j] * (2.0 * np.pi * 1j) * dRdq[i,k,l]))
    #                a[k,l] = tmp1.real
    #                #print(tmp1)
    #        ddf_val[i,j+3] = np.sum(a) 
    #        ddf_val[i+3,j] = np.sum(a)
    #print(np.linalg.det(ddf_val))
    end_ddfr = timer()
    print('     time for ddfr: ', start_ddfr-end_ddfr)
    end = timer()
    print(' time for derivative calculation: ', end-start)
    #print('df: ', df_val)
    #print('ddf: ', ddf_val)
    if np.linalg.cond(ddf_val) < 1/sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
    else:
        print('Derivative matrix is non-invertible! Stopping now...')
        exit()
    step = ddf_val_inv.dot(-df_val)
    return step,df_val,e1

def new_dFs2(mapin,xyz,vol):
    import numpy as np
    nx,ny,nz = mapin.shape
    dfs = np.zeros(shape=(nx,ny,nz,3),dtype=np.complex64)
    for i in range(3):  
        dfs[:,:,:,i] = np.fft.fftshift(
                                      (1/vol) * 2j * np.pi * 
                                      np.fft.fftn(mapin * xyz[i]))
    return dfs

def new_derivatives(e0,e1,w_grid,w2_grid,q,sv,xyz,xyz_sum,vol,dfrs=None):
    import numpy as np
    import sys
    import fcodes_fast
    nx, ny, nz = e0.shape
    start = timer()
    sv_np = np.zeros((nx,ny,nz,3), dtype=np.float64)
    for i in range(3):
        sv_np[:,:,:,i] = sv[i]
    dRdq = derivatives_wrt_q(q)
    if dfrs is None:
        #print('Calculating derivatives...')
        start = timer()
        dFRs = new_dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))),xyz,vol) 
        end = timer()
        if timeit: print(' time for dFRs calculation: ', end-start)    
    if dfrs is not None:
        print('Interpolated DFRs are used!')
        dFRs = dfrs
    start = timer()
    df_val,ddf_val = fcodes_fast.calc_derivatives(e0,e1,w_grid,w2_grid,sv_np,dFRs,dRdq,xyz_sum,vol,nx,ny,nz)
    end = timer()
    if timeit: print(' time for derivative calculation: ', end-start)
    if np.linalg.cond(ddf_val) < 1/sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
    else:
        print('Derivative matrix is non-invertible! Stopping now...')
        exit()
    step = ddf_val_inv.dot(-df_val)
    return step,df_val

def optimized_derivative_calc(e0,e1,w_grid,w2_grid,q,sv,xyz,xyz_sum,vol):
    import numpy as np
    import sys
    df_val = np.zeros(6,dtype=np.float64)
    ddf_val = np.zeros(shape=(6,6),dtype=np.float64)
    tp2 = (2.0 * np.pi)**2
    # Translation derivatives
    start = timer()
    for i in range(3):
        df_val[i] = np.real(np.sum(w_grid * e0 * 
                np.conjugate(e1 * (2.0 * np.pi * 1j) * sv[i]))) 
        for j in range(3):
            if i == 0 or (i > 0 and j >= i):
                ddf_val[i,j] = -tp2 * np.sum(w2_grid * sv[i] * sv[j])
            else:
                ddf_val[i,j] = ddf_val[j,i]
    
    # Rotation derivatives
    dRdq = derivatives_wrt_q(q) # give current q - vector
    start_dfs2 = timer()
    dFRs = dFs2(np.real(np.fft.ifftn(np.fft.ifftshift(e1))),xyz,vol)
    end_dfs2 = timer()
    print('     time for DFS calc: ', end_dfs2-start_dfs2)
    a = np.zeros(shape=(3,3),dtype=np.float64)
    b = np.zeros(shape=(3,3),dtype=np.float64)
    start_dfr = timer()
    for i in range(3):
        a[:,:] = 0.0
        for k in range(3):
            for l in range(3):
                if k == 0 or (k > 0 and l >= k):
                    a[k,l] = np.real(np.sum(w_grid * np.conjugate(e0) * 
                    (dFRs[:,:,:,k] * sv[l] * dRdq[i,k,l])))
                else:
                    a[k,l] = a[l,k]
        df_val[i+3] = np.sum(a) # df_val[3] to df_val[5]
    end_dfr = timer()
    print('     time for dfr: ', end_dfr-start_dfr)
    wfsc = np.real(w_grid * np.conjugate(e0) * e1) # THIS WORKS
    start_ddfr = timer()
    for i in range(3):
        for j in range(3):
            if i == 0 or (i > 0 and j >= i):
                b[:,:] = 0.0
                n = -1
                for k in range(3):
                    for l in range(3): 
                        if k == 0 or (k > 0 and l >= k):
                            n = n + 1
                            b[k,l] = -(tp2/vol) * xyz_sum[n] * np.sum(wfsc * sv[k] * sv[l] *
                                dRdq[i,k,l] * dRdq[j,k,l])
                        else:
                            b[k,l] = b[l,k]
                ddf_val[i+3,j+3] = np.sum(b) # ddf_val[3] to [5] in i and j indices
            else:
                ddf_val[i+3,j+3] = ddf_val[j+3,i+3]    
    end_ddfr = timer()
    print('     time for ddfr: ', start_ddfr-end_ddfr)
    end = timer()
    print(' time for derivative calculation: ', end-start)
    if np.linalg.cond(ddf_val) < 1/sys.float_info.epsilon:
        ddf_val_inv = np.linalg.pinv(ddf_val)
    else:
        print('Derivative matrix is non-invertible! Stopping now...')
        exit()
    step = ddf_val_inv.dot(-df_val)
    return step,df_val,e1