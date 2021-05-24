#!/usr/bin/env python
# -*- coding: utf-8 -*-
u'''
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright:
    * Until package version 1.1.23: © 2011 Daniel Müllner <http://danifold.net>
    * All changes from version 1.1.24 on: © Google Inc. <http://google.com>
'''
import os
import sys
import numpy
from setuptools import setup, Extension
from io import open

with open('fastcluster.py', encoding='utf_8') as f:
    for line in f:
        if line.find('__version_info__ =') == 0:
            version = '.'.join(line.split("'")[1:-1:2])
            break

print('Fastcluster version: ' + version)
print('Python version: ' + sys.version)

setup(name='fastcluster',
      version=version,
      py_modules=['fastcluster'],
      description='Fast hierarchical clustering routines for R and Python.',
      python_requires='>=3',
      requires=['numpy'],
      install_requires=["numpy>=1.9"],
      extras_require={'test':  ['scipy>=1.6.3']},
      provides=['fastcluster'],
      ext_modules=[Extension('_fastcluster',
                             ['src/fastcluster_python.cpp'],
                             extra_compile_args=['/EHsc'] if os.name == 'nt' else [],
                             include_dirs=[numpy.get_include()],
# Feel free to uncomment the line below if you use the GCC.
# This switches to more aggressive optimization and turns
# more warning switches on. No warning should appear in
# the compilation process.
#
# Also, the author's Python distribution generates debug
# symbols by default. This can be turned off, resulting a in
# much smaller compiled library.
#
# Optimization
#extra_compile_args=['-O2', '-g0', '-march=native', '-mtune=native', '-fno-math-errno'],
#
# List of all warning switches, somewhere from stackoverflow.com
#extra_compile_args=['-Wall', '-Weffc++', '-Wextra', '-Wall', '-Wcast-align', '-Wchar-subscripts', '-Wcomment', '-Wconversion', '-Wsign-conversion', '-Wdisabled-optimization', '-Wfloat-equal', '-Wformat', '-Wformat=2', '-Wformat-nonliteral', '-Wformat-security', '-Wformat-y2k', '-Wimport', '-Winit-self', '-Winline', '-Winvalid-pch', '-Wunsafe-loop-optimizations', '-Wmissing-braces', '-Wmissing-field-initializers', '-Wmissing-format-attribute', '-Wmissing-include-dirs', '-Wmissing-noreturn', '-Wpacked', '-Wparentheses', '-Wpointer-arith', '-Wredundant-decls', '-Wreturn-type', '-Wsequence-point', '-Wshadow', '-Wsign-compare', '-Wstack-protector', '-Wstrict-aliasing', '-Wstrict-aliasing=2', '-Wswitch', '-Wswitch-enum', '-Wtrigraphs', '-Wuninitialized', '-Wunknown-pragmas', '-Wunreachable-code', '-Wunused', '-Wunused-function', '-Wunused-label', '-Wunused-parameter', '-Wunused-value', '-Wunused-variable', '-Wvariadic-macros', '-Wvolatile-register-var', '-Wwrite-strings', '-Wlong-long', '-Wpadded', '-Wcast-qual', '-Wswitch-default', '-Wnon-virtual-dtor', '-Wold-style-cast', '-Woverloaded-virtual', '-Waggregate-return', '-Werror'],
#
# Linker optimization
#extra_link_args=['-Wl,--strip-all'],
      )],
      keywords=['dendrogram', 'linkage', 'cluster', 'agglomerative',
                'hierarchical', 'hierarchy', 'ward'],
      author=u"Daniel Müllner",
      author_email="daniel@danifold.net",
      license="BSD <http://opensource.org/licenses/BSD-2-Clause>",
      classifiers=[
          "Topic :: Scientific/Engineering :: Information Analysis",
          "Topic :: Scientific/Engineering :: Artificial Intelligence",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Mathematics",
          "Programming Language :: Python",
          "Programming Language :: Python :: 3",
          "Programming Language :: C++",
          "Operating System :: OS Independent",
          "License :: OSI Approved :: BSD License",
          "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
          "Intended Audience :: Science/Research",
          "Development Status :: 5 - Production/Stable"],
      url='http://danifold.net',
      test_suite='tests.fastcluster_test',
)
