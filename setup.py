#!/usr/bin/env python
# encoding: utf-8

import setuptools
from numpy.distutils.core import setup, Extension
import os, platform

path = 'src'
src = ['Poly.cpp', 'myMath.cpp', 'BeamFEA.cpp', 'Beam.cpp', 'pyBEAM.cpp']

for i in range(4):
    src[i] = os.path.join(path, 'pBEAM', src[i])
src[4] = os.path.join(path, 'pyBEAM', src[4])


if platform.system() == 'Windows':
    setup(
        name='pBEAM',
        version='0.2.0',
        description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
        author='S. Andrew Ning and Garrett E. Barter',
        author_email='garrett.barter@nrel.gov',
        license=['Apache License, Version 2.0','Mozilla Public License (MPL) version 2.0',
                 'Boost Software License 1.0','pybind11 license'],
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-std=gnu11','-fPIC'],
                                   include_dirs=[os.path.join(path, 'pBEAM'), os.path.join(path, 'include')])],
        zip_safe=False
    )
else:
    setup(
        name='pBEAM',
        version='0.2.0',
        description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
        author='S. Andrew Ning and Garrett E. Barter',
        author_email='garrett.barter@nrel.gov',
        license=['Apache License, Version 2.0','Mozilla Public License (MPL) version 2.0',
                 'Boost Software License 1.0','pybind11 license'],
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-std=c++11','-fPIC'],
                                   include_dirs=[os.path.join(path, 'pBEAM'), os.path.join(path, 'include')])],
        zip_safe=False
    )