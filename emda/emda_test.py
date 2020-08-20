"""
Author: "Rangana Warshamanage, Garib N. Murshudov"
MRC Laboratory of Molecular Biology

This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""

from __future__ import absolute_import, division, print_function, unicode_literals

from emda import iotools, maptools, restools
import fcodes_fast

def main():
    iotools.test()
    maptools.test()
    restools.test()
    fcodes_fast.test()

if __name__ == '__main__':
    main()
