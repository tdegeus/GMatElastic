desc = '''
Linear elastic material model.
'''

from setuptools import setup, Extension

import sys,re
import setuptools
import pybind11
import pyxtensor

header = open('include/GMatElastic/config.h','r').read()
world  = re.split(r'(.*)(\#define GMATELASTIC_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split(r'(.*)(\#define GMATELASTIC_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split(r'(.*)(\#define GMATELASTIC_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

ext_modules = [
  Extension(
    'GMatElastic',
    ['python/main.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      pyxtensor.get_include(False),
      pyxtensor.get_include(True ),
      pyxtensor.find_xtensor(),
      pyxtensor.find_xtl(),
    ],
    language='c++'
  ),
]

setup(
  name             = 'GMatElastic',
  description      = 'Linear elastic material model',
  long_description = desc,
  version          = __version__,
  license          = 'MIT',
  author           = 'Tom de Geus',
  author_email     = 'tom@geus.me',
  url              = 'https://github.com/tdegeus/GMatElastic',
  ext_modules      = ext_modules,
  install_requires = ['pybind11>=2.2.0'],
  cmdclass         = {'build_ext': pyxtensor.BuildExt},
  zip_safe         = False,
)
