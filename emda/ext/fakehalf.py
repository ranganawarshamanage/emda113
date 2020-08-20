# Neighbor averaging method based fake half map
def fakehalf(f_map, knl=3):
    import numpy as np
    import scipy.signal
    #kernel = np.ones((knl, knl, knl), dtype=int)
    #kernel[1,1,1] = 0
    #kernel = kernel/np.sum(kernel)
    box_size = knl
    kernel = np.ones((box_size, box_size, box_size)) / (box_size ** 3)
    kernel[1,1,1] = 0
    kernel = kernel/np.sum(kernel)
    fakehalf = scipy.signal.fftconvolve(f_map, kernel, "same")
    return fakehalf