#!/usr/bin/env python
# encoding: utf-8

import setuptools
from numpy.distutils.core import setup, Extension
import os, platform

path = 'src'
src = ['Poly.cpp', 'myMath.cpp', 'BeamFEA.cpp', 'CurveFEM.cpp', 'Beam.cpp', 'pyBEAM.cpp']

for i in range(len(src)-1):
    src[i] = os.path.join(path, 'pBEAM', src[i])
src[-1] = os.path.join(path, 'pyBEAM', src[-1])


if platform.system() == 'Windows':
    # Note: must use mingw compiler on windows or a Visual C++ compiler version that supports std=c++11
    arglist = ['-std=gnu++11','-fPIC']
else:
    arglist = ['-std=c++11','-fPIC']
    
setup(
    name='pBEAM',
    version='0.2.0',
    description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
    author='NREL WISDEM Team and Garrett E. Barter',
    author_email='systems.engineering@nrel.gov',
    license=['Apache License, Version 2.0','Mozilla Public License (MPL) version 2.0','pybind11 license'],
    ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=arglist,
                           include_dirs=[os.path.join(path, 'pBEAM'), os.path.join(path, 'include')])],
    zip_safe=False
)
