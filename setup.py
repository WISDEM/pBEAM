#!/usr/bin/env python
# encoding: utf-8

import setuptools
from numpy.distutils.core import setup, Extension
import os
import platform

path = 'src'
src = ['Poly.cpp', 'myMath.cpp', 'BeamFEA.cpp', 'Beam.cpp', 'pyBEAM.cpp']

for i in range(4):
    src[i] = os.path.join(path, 'pBEAM', src[i])
src[4] = os.path.join(path, 'pyBEAM', src[4])


# f = open('MANIFEST.in', 'a')
# f.write('recursive-include ' + path + ' * \n')
# f.write('recursive-exclude ' + os.path.join(path, 'pBEAM.xcodeproj') + ' * \n')
# f.write('exclude ' + os.path.join(path, 'pBEAM', 'main.cpp') + '\n')
# f.close()

#os.environ["CC"] = "g++"
#os.environ["CXX"] = "g++"

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
                               include_dirs=[os.path.join(path, 'pBEAM'), 'C:/boost_1_55_0'],
                               library_dirs=['C:/boost_1_55_0/stage/lib', 'C:/lapack'],
                               libraries=['boost_python-mgw46-mt-1_55', 'lapack', 'blas'])])
elif platform.system() == 'Darwin':
    setup(
        name='pBEAM',
        version='0.1.0',
        description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
        author='S. Andrew Ning',
        author_email='andrew.ning@nrel.gov',
        # install_requires=['numpy', 'scipy'],
        license='Apache License, Version 2.0',
        # OS X, Linux
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-O2'],
                               include_dirs=[os.path.join(path, 'pBEAM'),'/opt/local/include'],
                               library_dirs=['/opt/local/lib', '/opt/local/lib/lapack/'],
                               libraries=['boost_python-mt', 'boost_numpy-mt', 'lapack', 'blas'])])
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
        ext_modules=[Extension('_pBEAM', sources=src, extra_compile_args=['-O2'],
                               include_dirs=[os.path.join(path, 'pBEAM')],
                               libraries=['boost_python','boost_numpy','lapack'])])

