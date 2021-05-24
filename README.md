This library provides Python functions for hierarchical clustering. It generates
hierarchical clusters from distance matrices or from vector data.

This module is intended to replace the functions

    linkage, single, complete, average, weighted, centroid, median, ward

in the module [`scipy.cluster.hierarchy`](
https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html) with the same
functionality but much faster algorithms. Moreover, the function
`linkage_vector` provides memory-efficient clustering for vector data.

The interface is very similar to MATLAB's Statistics Toolbox API to make code
easier to port from MATLAB to Python/NumPy. The core implementation of this
library is in C++ for efficiency.

**User manual:** [fastcluster.pdf](
https://github.com/dmuellner/fastcluster/raw/master/docs/fastcluster.pdf).

The “Yule” distance function changed in fastcluster version 1.2.0. This is
following a [change in SciPy 1.6.3](
https://github.com/scipy/scipy/commit/3b22d1da98dc1b5f64bc944c21f398d4ba782bce).
It is recommended to use fastcluster version 1.1.x together with SciPy versions
before 1.6.3 and fastcluster 1.2.x with SciPy ≥1.6.3.

The fastcluster package is considered stable and will undergo few changes
from now on. If some years from now there have not been any updates, this does
not necessarily mean that the package is unmaintained but maybe it just was
not necessary to correct anything. Of course, please still report potential
bugs and incompatibilities to daniel@danifold.net. You may also use
[my GitHub repository](https://github.com/dmuellner/fastcluster/)
for bug reports, pull requests etc.

Note that [PyPI](https://pypi.org/project/fastcluster/) and [my GitHub
repository](https://github.com/dmuellner/fastcluster/) host the source code
for the Python interface only. The archive with both the R and the Python
interface is available on
[CRAN](https://CRAN.R-project.org/package=fastcluster) and the GitHub repository
[“cran/fastcluster”](https://github.com/cran/fastcluster). Even though I appear
as the author also of this second GitHub repository, this is just an automatic,
read-only mirror of the CRAN archive, so please do not attempt to report bugs or
contact me via this repository.

Installation files for Windows are provided on [PyPI](
https://pypi.org/project/fastcluster/#files) and on [Christoph Gohlke's web
page](http://www.lfd.uci.edu/~gohlke/pythonlibs/#fastcluster).

Christoph Dalitz wrote a pure [C++ interface to fastcluster](
https://lionel.kr.hs-niederrhein.de/~dalitz/data/hclust/).

Reference: Daniel Müllner, *fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python*, Journal of Statistical Software, **53**
(2013), no. 9, 1–18, https://www.jstatsoft.org/v53/i09/.
