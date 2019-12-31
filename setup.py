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
      long_description=u"""
This library provides Python functions for hierarchical clustering. It
generates hierarchical clusters from distance matrices or from vector data.

Part of this module is intended to replace the functions ::

    linkage, single, complete, average, weighted, centroid, median, ward

in the module ``scipy.cluster.hierarchy`` with the same functionality but much
faster algorithms. Moreover, the function ``linkage_vector`` provides
memory-efficient clustering for vector data.

The interface is very similar to MATLAB's Statistics Toolbox API to make code
easier to port from MATLAB to Python/NumPy. The core implementation of this
library is in C++ for efficiency.

**User manual:** `fastcluster.pdf
<https://github.com/dmuellner/fastcluster/raw/master/docs/fastcluster.pdf>`_.

Installation files for Windows are provided on `PyPI
<https://pypi.org/project/fastcluster/#files>`__ and on `Christoph Gohlke's
web page <http://www.lfd.uci.edu/~gohlke/pythonlibs/#fastcluster>`_.

**The fastcluster package is considered stable and will undergo few changes
from now on. If some years from now there have not been any updates, this
does not necessarily mean that the package is unmaintained but maybe it just
was not necessary to correct anything. Of course, please still report potential
bugs and incompatibilities to daniel@danifold.net. You may also use**
`my GitHub repository <https://github.com/dmuellner/fastcluster/>`_
**for bug reports, pull requests etc.**

Note that `PyPI <https://pypi.org/project/fastcluster/>`__ and `my GitHub
repository <https://github.com/dmuellner/fastcluster/>`_ host the source code
for the Python interface only. The archive with both the R and the Python
interface is available on `CRAN
<https://CRAN.R-project.org/package=fastcluster>`_ and the GitHub repository
`“cran/fastcluster” <https://github.com/cran/fastcluster>`_. Even though I
appear as the author also of this second GitHub repository, this is just an
automatic, read-only mirror of the CRAN archive, so please do not attempt to
report bugs or contact me via this repository.

Christoph Dalitz wrote a pure `C++ interface to fastcluster
<http://informatik.hsnr.de/~dalitz/data/hclust>`_.

Reference: Daniel Müllner, *fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python*, Journal of Statistical Software, **53**
(2013), no. 9, 1–18, http://www.jstatsoft.org/v53/i09/.
""",
      long_description_content_type='text/x-rst',
      requires=['numpy'],
      install_requires=["numpy>=1.9"],
      extras_require={'test':  ['scipy']},
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
          "Programming Language :: Python :: 2",
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
