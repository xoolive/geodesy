#! /usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

import sys

INCLUDE = []
CARGS = []

try:
    import numpy
    INCLUDE.append(numpy.get_include())
except ImportError:
    print ("numpy is required in order to use geodesy toolkit")
    raise

if not sys.platform == "win32":
    # suppress warnings for importing numpy
    CARGS.append("-Wno-unused-function")

extensions = [
    Extension("geodesy.sphere",
              ["geodesy/sphere.pyx"],
              include_dirs=INCLUDE, language="c++", extra_compile_args=CARGS),
    Extension("geodesy.wgs84",
              ["geodesy/wgs84.pyx"],
              include_dirs=INCLUDE, language="c++", extra_compile_args=CARGS),
]

setup(name="geodesy",
      version="1.1",
      author="Xavier Olive",
      author_email="xavier@xoolive.org",
      description="A pragmatic and efficient geodesy toolkit",
      ext_modules=cythonize(extensions),
      license="MIT",
      packages=['geodesy', ],
      url='https://github.com/xoolive/geodesy',
      )
