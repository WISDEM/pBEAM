#!/usr/bin/env python
# encoding: utf-8

import setuptools
from numpy.distutils.core import setup, Extension
from os.path import join
import platform

path = 'src'
src = ['Poly.cpp', 'myMath.cpp', 'BeamFEA.cpp', 'Beam.cpp', 'pyBEAM.cpp']

for i in range(4):
    src[i] = join(path, 'pBEAM', src[i])
src[4] = join(path, 'pyBEAM', src[4])


# f = open('MANIFEST.in', 'a')
# f.write('recursive-include ' + path + ' * \n')
# f.write('recursive-exclude ' + join(path, 'pBEAM.xcodeproj') + ' * \n')
# f.write('exclude ' + join(path, 'pBEAM', 'main.cpp') + '\n')
# f.close()

if platform.system() == 'Windows':
    setup(
        name='pBEAM',
        version='0.1.0',
        description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
        author='S. Andrew Ning',
        author_email='andrew.ning@nrel.gov',
        # install_requires=['numpy', 'scipy'],
        license='Apache License, Version 2.0',
        # Windows
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-O2'],
                               include_dirs=[join(path, 'pBEAM'), 'C:/boost_1_55_0'],
                               library_dirs=['C:/boost_1_55_0/stage/lib', 'C:/lapack'],
                               libraries=['boost_python-mgw46-mt-1_55', 'lapack'])])
else:
    setup(
        name='pBEAM',
        version='0.1.0',
        description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
        author='S. Andrew Ning',
        author_email='andrew.ning@nrel.gov',
        # install_requires=['numpy', 'scipy'],
        license='Apache License, Version 2.0',
        # OS X, Linux
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-O2'])])