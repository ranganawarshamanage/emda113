from __future__ import division, absolute_import, print_function
import setuptools
from numpy.distutils.core import setup, Extension

ex1 = Extension(name = 'fcodes_fast',
                sources = ['emda/fcodes_fast.f90'])

setup(name= 'emda',
    version= '1.1.3',
    description= 'Electron Microscopy map and model manipulation tools',
    long_description='''\
        Python library for Electron Microscopy map and model
        manipulations. To work with MRC/CCP4.MAP and MTZ files. Map to map
        fitting and average/difference map calculation. Map and map-model validation
        using 3D correlations.''',
    url= None,
    author= 'Rangana Warshamanage, Garib N. Murshudov',
    author_email= 'ranganaw@mrc-lmb.cam.ac.uk, garib@mrc-lmb.cam.ac.uk',
    license= 'MPL-2.0',
    packages= ['emda','emda.mapfit'],
    ext_modules = [ex1],
    install_requires=['pandas==0.23.4','mrcfile','matplotlib','numpy','scipy','gemmi'],
    test_suite='emda.tests',
    entry_points={
      'console_scripts': [
          'emda_test = emda.emda_test:main',
          'emda = emda.emda_cmd_caller:main',
                          ],
      },
    zip_safe= False)
