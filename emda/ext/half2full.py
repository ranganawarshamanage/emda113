"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import argparse
import sys
from emda import iotools, restools, fsc,plotter
from emda.config import *

# Fullmap from halfmaps

cmdl_parser = argparse.ArgumentParser(description='Fullmap from halfmaps\n')
cmdl_parser.add_argument('-h1', '--half1', required=True, help='input halfmap1')
cmdl_parser.add_argument('-h2', '--half2', required=True, help='input halfmap2')
cmdl_parser.add_argument('-o', '--outfull', required=False, help='output fullmap')

def half2full(ar1, ar2):
    import numpy as np
    assert ar1.shape == ar2.shape
    hf1 = np.fft.fftn(ar1)
    hf2 = np.fft.fftn(ar2)    
    return np.real(np.fft.ifftn((hf1 + hf2)/2.0))

def main():
    args = cmdl_parser.parse_args()
    uc, ar, org = iotools.read_map(args.half1)
    uc, ar2, org = iotools.read_map(args.half2)
    assert ar.shape == ar2.shape
    outfile = 'fullmap.mrc'
    if args.outfull is not None: outfile=args.outfull
    hf1 = np.fft.fftn(ar)
    hf2 = np.fft.fftn(ar2)
    iotools.write_mrc(np.real(np.fft.ifftn((hf1 + hf2)/2.0)),outfile,uc)
    print('Fullmap was written! ', outfile)


if(__name__ == "__main__"):
    main()

