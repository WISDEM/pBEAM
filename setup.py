#!/usr/bin/env python
# encoding: utf-8

import setuptools
from setuptools import setup, Extension
import os
import platform
import pybind11

# Define source files and paths
path = 'src'
src = ['Poly.cpp', 'myMath.cpp', 'BeamFEA.cpp', 'CurveFEM.cpp', 'Beam.cpp', 'pyBEAM.cpp']

# Construct full paths for source files
src = [os.path.join(path, 'pBEAM', f) if f != 'pyBEAM.cpp' else os.path.join(path, 'pyBEAM', f) for f in src]

# Platform-specific compiler arguments
if platform.system() == 'Windows':
    # Note: Use mingw or a Visual C++ compiler supporting C++11
    extra_compile_args = ['-std=gnu++11', '-fPIC']
else:
    extra_compile_args = ['-std=c++11', '-fPIC']

# Define the extension module
pbeam_extension = Extension(
    '_pBEAM',
    sources=src,
    include_dirs=[
        os.path.join(path, 'pBEAM'),
        os.path.join(path, 'include'),
        pybind11.get_include()  # Dynamically include pybind11 headers
    ],
    extra_compile_args=extra_compile_args,
    language='c++'
)

setup(
    name='pBEAM',
    version='0.2.0',
    description='Polynomial Beam Element Analysis Module. Finite element analysis for beam-like structures.',
    author='NREL WISDEM Team and Garrett E. Barter',
    author_email='systems.engineering@nrel.gov',
    license='Apache License, Version 2.0; Mozilla Public License (MPL) version 2.0; pybind11 license',
    ext_modules=[pbeam_extension],
    zip_safe=False,
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.20',
        'pybind11>=2.9.0'  # Ensure pybind11 version supports Python 3.11
    ]
)