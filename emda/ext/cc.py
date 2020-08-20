

def cc_overall_realsp(map1, map2):
    # Overall correlation - real space
    import numpy as np
    import math
    mean1 = np.mean(map1)
    mean2 = np.mean(map2)
    covar = np.sum(map1 * map2 - (mean1 * mean2))
    var1 = np.sum(map1 * map1 - (mean1 * mean1))
    var2 = np.sum(map2 * map2 - (mean2 * mean2))
    occ = covar / math.sqrt(var1 * var2)
    return 2.0 * occ / (1.0 + occ)

def cc_overall_fouriersp(f1, f2):
    import numpy as np
    import math
    f1_mean = np.mean(f1)
    f2_mean = np.mean(f2)
    covar = np.sum(f1 * np.conjugate(f2) - f1_mean * np.conjugate(f2_mean))
    var1 = np.sum(f1 * np.conjugate(f1) - f1_mean * np.conjugate(f1_mean))
    var2 = np.sum(f2* np.conjugate(f2) - f2_mean * np.conjugate(f2_mean))
    fsc = (covar.real / math.sqrt(var1.real * var2.real))
    return 2.0 * fsc / (1.0 + fsc)
