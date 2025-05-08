# *fastcluster*: Fast hierarchical clustering routines for R and Python

The fastcluster package is a C++ library for hierarchical, agglomerative
clustering. It efficiently implements the seven most widely used clustering
schemes: single, complete, average, weighted/McQuitty, Ward, centroid and
median linkage. The library has interfaces to two languages: R and Python.

The Python module is designed to replace the functions
```
linkage, single, complete, average, weighted, centroid, median, ward
```
in the module [`scipy.cluster.hierarchy`][] with the same functionality but
faster algorithms. Moreover, the function `linkage_vector` provides
memory-efficient clustering for vector data.

The R package is meant to replace [`hclust`][] in the [`stats`][] package
and the [`flashClust`][] package.

[`linkage`]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
[`scipy.cluster.hierarchy`]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
[`hclust`]: https://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
[`stats`]: https://stat.ethz.ch/R-manual/R-patched/library/stats/html/00Index.html
[`flashClust`]: https://CRAN.R-project.org/package=flashClust

See the [author's home page][] for more information, in particular a performance
comparison with other clustering packages. The user's manual is the file
[docs/fastcluster.pdf][] in the source distribution.

[author's home page]: https://danifold.net
[docs/fastcluster.pdf]: https://raw.githubusercontent.com/fastcluster/fastcluster/master/docs/fastcluster.pdf

## Distribution

The distributions on [GitHub][] and [PyPi][] contain only the files for the
Python interface. The full source distribution with both interfaces is
available on CRAN: <https://CRAN.R-project.org/package=fastcluster>.

The Python package can be installed from PyPI (conveniently with [pip][]),
GitHub, or from the source package at CRAN. All distributions compile and
install the same Python/C++ libraries.

> [!NOTE]
> The following sections describe the Python interface. See the R package on
> CRAN for the documentation of the R interface.

[GitHub]: https://github.com/fastcluster/fastcluster/
[PyPi]: https://pypi.org/project/fastcluster/
[pip]: https://pip.pypa.io


## Quick installation

```
pip install fastcluster
```

## Usage

The fastcluster module is imported as usual by
```
import fastcluster
```
It provides the following functions:
```
linkage(X, method='single', metric='euclidean', preserve_input=True)
single(X)
complete(X)
average(X)
weighted(X)
ward(X)
centroid(X)
median(X)
linkage_vector(X, method='single', metric='euclidean', extraarg=None)
```
The argument `X` is either a compressed distance matrix or a collection of *n*
observation vectors in *d* dimensions as an (*n*×*d*) array. Apart from the
argument `preserve_input`, the methods have the same input and output as the
functions of the same name in the module `scipy.cluster.hierarchy`.

The optional argument `preserve_input` specifies whether the fastcluster package
first copies the distance matrix or writes into the existing array. If the
dissimilarities are generated for the clustering step only and are not
needed afterward, approximately half the memory can be saved by specifying
`preserve_input=False`. Note that the input array `X` contains unspecified
values after this procedure. You may want to write
```
linkage(X, method='…', preserve_input=False)
del X
```
to make sure that the matrix `X` is not accessed accidentally after it has been
used as scratch memory.

The method
```
linkage_vector(X, method='single', metric='euclidean', extraarg=None)
```
provides memory-saving clustering for vector data. It also accepts a collection
of *n* observation vectors in *d* dimensions as an (*n*×*d*) array as the first
parameter. The parameter `method` is either `single`, `ward`, `centroid` or
`median`. The `ward`, `centroid` and `median` methods require the Euclidean
metric. In case of single linkage, the `metric` parameter can be chosen from
all metrics which are implemented in [`scipy.spatial.dist.pdist`][]. There may be differences between
```
linkage(scipy.spatial.dist.pdist(X, metric='…'))
```
and
```
linkage_vector(X, metric='…')
```
since a few corrections have been made compared to the pdist function. Please
consult the [user's manual] for comprehensive details.

[`scipy.spatial.dist.pdist`]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
[user's manual]: https://raw.githubusercontent.com/fastcluster/fastcluster/master/docs/fastcluster.pdf

## Copyright

  * Until package version 1.1.23: © 2011 Daniel Müllner <https://danifold.net>
  * All changes from version 1.1.24 on: © Google Inc. <https://www.google.com>

## License

The fastcluster package is distributed under the [BSD license][]. See the file
[LICENSE][] in the source distribution.

[BSD License]: https://opensource.org/licenses/BSD-2-Clause
[LICENSE]: https://raw.githubusercontent.com/fastcluster/fastcluster/refs/heads/master/LICENSE

## Citation

To cite fastcluster in publications, please use:

Daniel Müllner, fastcluster: Fast Hierarchical, Agglomerative Clustering
Routines for R and Python, Journal of Statistical Software, 53 (2013), no. 9,
1–18, https://doi.org/10.18637/jss.v053.i09.

## Further links

* Project home page:
  [https://danifold.net/fastcluster](https://danifold.net/fastcluster)
* Christoph Dalitz wrote a pure C++ interface to fastcluster:
  <https://lionel.kr.hs-niederrhein.de/~dalitz/data/hclust/>.
