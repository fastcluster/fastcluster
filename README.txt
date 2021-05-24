fastcluster: Fast hierarchical clustering routines for R and Python

Copyright:
  * Until package version 1.1.23: © 2011 Daniel Müllner <http://danifold.net>
  * All changes from version 1.1.24 on: © Google Inc. <http://google.com>

The fastcluster package is a C++ library for hierarchical, agglomerative
clustering. It efficiently implements the seven most widely used clustering
schemes: single, complete, average, weighted/McQuitty, Ward, centroid and
median linkage. The library currently has interfaces to two languages: R and
Python/NumPy. Part of the functionality is designed as drop-in replacement for
existing routines: “linkage” in the SciPy package “scipy.cluster.hierarchy”,
“hclust” in R's “stats” package, and the “flashClust” package. Once the
fastcluster library is loaded at the beginning of the code, every program that
uses hierarchical clustering can benefit immediately and effortlessly from the
performance gain. Moreover, there are memory-saving routines for clustering of
vector data, which go beyond what the existing packages provide.

See the author's home page <http://danifold.net> for more
information, in particular a performance comparison with other clustering
packages. The User's manual is the file docs/fastcluster.pdf in the
source distribution.

The fastcluster package is distributed under the BSD license. See the file
LICENSE in the source distribution or
<http://opensource.org/licenses/BSD-2-Clause>.

The distribution on pypi.python.org contains only the files which are necessary
for the Python interface. The full source distribution with both interfaces
is available on CRAN

    https://CRAN.R-project.org/package=fastcluster

The Python package can be installed either from PyPI (conveniently with pip)
or manually from the source package at CRAN. Both distributions compile and
install identical libraries.

Christoph Dalitz wrote a pure C++ interface to fastcluster:
<https://lionel.kr.hs-niederrhein.de/~dalitz/data/hclust/>.


Installation
‾‾‾‾‾‾‾‾‾‾‾‾
See the file INSTALL.txt in the source distribution, which also explains how to
install the fastcluster package for R.


Usage
‾‾‾‾‾
1. R
‾‾‾‾
In R, load the package with the following command:

    library('fastcluster')

The package overwrites the function hclust from the “stats” package (in the
same way as the flashClust package does). Please remove any references to the
flashClust package in your R files to not accidentally overwrite the hclust
function with the flashClust version.

The new hclust function has exactly the same calling conventions as the old
one. You may just load the package and immediately and effortlessly enjoy the
performance improvements. The function is also an improvement to the flashClust
function from the “flashClust” package. Just replace every call to flashClust
by hclust and expect your code to work as before, only faster. (If you are
using flashClust prior to version 1.01, update it! See the change log for
flashClust:

    http://cran.r-project.org/web/packages/flashClust/ChangeLog )

If you need to access the old function or make sure that the right function is
called, specify the package as follows:

    fastcluster::hclust(…)
    flashClust::hclust(…)
    stats::hclust(…)

Vector data can be clustered with a memory-saving algorithm with the command

    hclust.vector(…)

See the User's manual docs/fastcluster.pdf for further details.

WARNING
‾‾‾‾‾‾‾
R and Matlab/SciPy use different conventions for the “Ward”, “centroid” and
“median” methods. R assumes that the dissimilarity matrix consists of squared
Euclidean distances, while Matlab and SciPy expect non-squared Euclidean
distances. The fastcluster package respects these conventions and uses
different formulas in the two interfaces.

If you want the same results in both interfaces, then feed the hclust function
in R with the entry-wise square of the distance matrix, D^2, for the “Ward”,
“centroid” and “median” methods and later take the square root of the height
field in the dendrogram. For the “average” and “weighted” alias “mcquitty”
methods, you must still take the same distance matrix D as in the Python
interface for the same results. The “single” and “complete” methods only depend
on the relative order of the distances, hence it does not make a difference
whether the method operates on the distances or the squared distances.

The code example in the R documentation (enter ?hclust or example(hclust) in R)
contains an instance where the squared distance matrix is generated from
Euclidean data.

2. Python
‾‾‾‾‾‾‾‾‾
The fastcluster package is imported as usual by

    import fastcluster

It provides the following functions:

    linkage(X, method='single', metric='euclidean', preserve_input=True)
    single(X)
    complete(X)
    average(X)
    weighted(X)
    ward(X)
    centroid(X)
    median(X)
    linkage_vector(X, method='single', metric='euclidean', extraarg=None)

The argument X is either a compressed distance matrix or a collection of n
observation vectors in d dimensions as an (n×d) array. Apart from the argument
preserve_input, the methods have the same input and output as the functions of
the same name in the package scipy.cluster.hierarchy.

The additional, optional argument preserve_input specifies whether the
fastcluster package first copies the distance matrix or writes into the
existing array. If the dissimilarities are generated for the clustering step
only and are not needed afterward, approximately half the memory can be saved
by specifying preserve_input=False. Note that the input array X contains
unspecified values after this procedure. You may want to write

    linkage(X, method='…', preserve_input=False)
    del X

to make sure that the matrix X is not accessed accidentally after it has been
used as scratch memory.

The method

    linkage_vector(X, method='single', metric='euclidean', extraarg=None)

provides memory-saving clustering for vector data. It also accepts a collection
of n observation vectors in d dimensions as an (n×d) array as the first parameter.
The parameter 'method' is either 'single', 'ward', 'centroid' or 'median'. The
'ward', 'centroid' and 'median' methods require the Euclidean metric. In case
of single linkage, the 'metric' parameter can be chosen from all metrics which
are implemented in scipy.spatial.dist.pdist. There may be differences between

    linkage(scipy.spatial.dist.pdist(X, metric='…'))
and
    linkage_vector(X, metric='…')

since there have been made a few corrections compared to the pdist function.
Please consult the the User's manual docs/fastcluster.pdf for
comprehensive details.
