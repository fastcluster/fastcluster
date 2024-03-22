/*
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright:
    * Until package version 1.1.23: © 2011 Daniel Müllner <https://danifold.net>
    * All changes from version 1.1.24 on: © Google Inc. <https://www.google.com>
*/

// for INT32_MAX in fastcluster.cpp
// This must be defined here since Python.h loads the header file pyport.h,
// and from this stdint.h. INT32_MAX is defined in stdint.h, but only if
// __STDC_LIMIT_MACROS is defined.
#define __STDC_LIMIT_MACROS

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 6))
#define HAVE_DIAGNOSTIC 1
#endif

#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wlong-long"
#pragma GCC diagnostic ignored "-Wformat"
#endif
#include <Python.h>
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wlong-long"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif
#include <numpy/arrayobject.h>
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif

/* It's complicated, but if I do not include the C++ math headers, GCC
   will complain about conversions from 'double' to 'float', whenever 'isnan'
   is called in a templated function (but not outside templates).

   The '#include <cmath>' seems to cure the problem.
*/
// #include <cmath>
#define fc_isnan(X) ((X) != (X))

// There is Py_IS_NAN but it is so much slower on my x86_64 system with GCC!

#include <algorithm> // for std::stable_sort
#include <cmath>     // for std::abs, std::pow, std::sqrt
#include <cstddef>   // for std::ptrdiff_t
#include <exception> // for std::exception
#include <limits>    // for std::numeric_limits<...>::infinity()
#include <new>       // for std::bad_alloc

#include "fastcluster.cpp"

// backwards compatibility
#ifndef NPY_ARRAY_CARRAY_RO
#define NPY_ARRAY_CARRAY_RO NPY_CARRAY_RO
#endif

/* Since the public interface is given by the Python respectively R interface,
 * we do not want other symbols than the interface initalization routines to be
 * visible in the shared object file. The "visibility" switch is a GCC concept.
 * Hiding symbols keeps the relocation table small and decreases startup time.
 * See http://gcc.gnu.org/wiki/Visibility
 */
#if HAVE_VISIBILITY
#pragma GCC visibility push(hidden)
#endif

/*
  Convenience class for the output array: automatic counter.
*/
template <typename T> class linkage_output {
private:
  T *Z;

public:
  linkage_output(T *const Z_) : Z(Z_) {}

  void append(const t_index node1, const t_index node2, const T dist,
              const T size) {
    if (node1 < node2) {
      *(Z++) = static_cast<T>(node1);
      *(Z++) = static_cast<T>(node2);
    } else {
      *(Z++) = static_cast<T>(node2);
      *(Z++) = static_cast<T>(node1);
    }
    *(Z++) = dist;
    *(Z++) = size;
  }
};

/*
  Generate the SciPy-specific output format for a dendrogram from the
  clustering output.

  The list of merging steps can be sorted or unsorted.
*/
// The size of a node is either 1 (a single point) or is looked up from
// one of the clusters.
#define size_(r_) (((r_ < N) ? 1 : Z_(r_ - N, 3)))

template <const bool sorted, typename T>
static void generate_SciPy_dendrogram(T *const Z, cluster_result<T> &Z2,
                                      const t_index N) {
  // The array "nodes" is a union-find data structure for the cluster
  // identities (only needed for unsorted cluster_result input).
  union_find nodes(sorted ? 0 : N);
  if (!sorted) {
    std::stable_sort(Z2[0], Z2[N - 1]);
  }

  linkage_output<T> output(Z);
  t_index node1, node2;

  for (node<T> const *NN = Z2[0]; NN != Z2[N - 1]; ++NN) {
    // Get two data points whose clusters are merged in step i.
    if (sorted) {
      node1 = NN->node1;
      node2 = NN->node2;
    } else {
      // Find the cluster identifiers for these points.
      node1 = nodes.Find(NN->node1);
      node2 = nodes.Find(NN->node2);
      // Merge the nodes in the union-find data structure by making them
      // children of a new node.
      nodes.Union(node1, node2);
    }
    output.append(node1, node2, NN->dist, size_(node1) + size_(node2));
  }
}

/*
  Python interface code
*/
static PyObject *linkage_wrap(PyObject *const self, PyObject *const args);
static PyObject *linkage_vector_wrap(PyObject *const self,
                                     PyObject *const args);

// List the C++ methods that this extension provides.
static PyMethodDef _fastclusterWrapMethods[] = {
    {"linkage_wrap", linkage_wrap, METH_VARARGS, NULL},
    {"linkage_vector_wrap", linkage_vector_wrap, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL} /* Sentinel - marks the end of this structure */
};

/* Tell Python about these methods.

   Python 2.x and 3.x differ in their C APIs for this part.
*/
#if PY_VERSION_HEX >= 0x03000000

static struct PyModuleDef fastclustermodule = {
    PyModuleDef_HEAD_INIT,
    "_fastcluster",
    NULL, // no module documentation
    -1,   /* size of per-interpreter state of the module,
             or -1 if the module keeps state in global variables. */
    _fastclusterWrapMethods,
    NULL,
    NULL,
    NULL,
    NULL};

/* Make the interface initalization routines visible in the shared object
 * file.
 */
#if HAVE_VISIBILITY
#pragma GCC visibility push(default)
#endif

PyMODINIT_FUNC PyInit__fastcluster(void) {
  PyObject *m;
  m = PyModule_Create(&fastclustermodule);
  if (!m) {
    return NULL;
  }
  import_array(); // Must be present for NumPy. Called first after above line.
  return m;
}

#if HAVE_VISIBILITY
#pragma GCC visibility pop
#endif

#else // Python 2.x

#if HAVE_VISIBILITY
#pragma GCC visibility push(default)
#endif

PyMODINIT_FUNC init_fastcluster(void) {
  (void)Py_InitModule("_fastcluster", _fastclusterWrapMethods);
  import_array(); // Must be present for NumPy. Called first after above line.
}

#if HAVE_VISIBILITY
#pragma GCC visibility pop
#endif

#endif // PY_VERSION

class GIL_release {
private:
  // noncopyable
  GIL_release(GIL_release const &);
  GIL_release &operator=(GIL_release const &);

public:
  inline GIL_release(bool really = true)
      : _save(really ? PyEval_SaveThread() : NULL) {}

  inline ~GIL_release() {
    if (_save)
      PyEval_RestoreThread(_save);
  }

private:
  PyThreadState *_save;
};

/*
  Interface to Python, part 1:
  The input is a dissimilarity matrix.
*/

static PyObject *linkage_wrap(PyObject *const, PyObject *const args) {
  PyArrayObject *D, *Z;
  long int N_ = 0;
  unsigned char method;
  int high_precision_flag = 0;

  try {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
    // Parse the input arguments
    if (!PyArg_ParseTuple(args, "lO!O!bp",
                          &N_,                     // signed long integer
                          &PyArray_Type, &D,       // NumPy array
                          &PyArray_Type, &Z,       // NumPy array
                          &method,                 // unsigned char
                          &high_precision_flag)) { // int
      return NULL; // Error if the arguments have the wrong type.
    }

    bool high_precision = high_precision_flag != 0;
    int T_mant_dig = high_precision ? DBL_MANT_DIG : FLT_MANT_DIG;

#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    if (N_ < 1) {
      // N must be at least 1.
      PyErr_SetString(PyExc_ValueError,
                      "At least one element is needed for clustering.");
      return NULL;
    }

    /*
      (1)
      The biggest index used below is 4*(N-2)+3, as an index to Z. This must
      fit into the data type used for indices.
      (2)
      The largest representable integer, without loss of precision, by a
      floating point number of type T is 2^T_mant_dig. Here, we
      make sure that all cluster labels from 0 to 2N-2 in the output can be
      accurately represented by a floating point number.

      Conversion of N to 64 bits below is not really necessary but it prevents
      a warning ("shift count >= width of type") on systems where "long int"
      is 32 bits wide.
    */
    if (N_ > MAX_INDEX / 4 ||
        static_cast<int64_t>(N_ - 1) >> (T_mant_dig - 1) > 0) {
      PyErr_SetString(PyExc_ValueError, "Data is too big, index overflow.");
      return NULL;
    }

    t_index N = static_cast<t_index>(N_);

    // Allow threads!
    GIL_release G;

    if (high_precision) {
      using T = double;
      T *const D_ = reinterpret_cast<T *>(PyArray_DATA(D));
      cluster_result<T> Z2(N - 1);
      auto_array_ptr<t_index> members;
      // For these methods, the distance update formula needs the number of
      // data points in a cluster.
      if (method == METHOD_METR_AVERAGE || method == METHOD_METR_WARD ||
          method == METHOD_METR_CENTROID) {
        members.init(N, 1);
      }
      // Operate on squared distances for these methods.
      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID ||
          method == METHOD_METR_MEDIAN) {
        for (T *DD = D_;
             DD != D_ + static_cast<std::ptrdiff_t>(N) * (N - 1) / 2; ++DD)
          *DD *= *DD;
      }

      switch (method) {
      case METHOD_METR_SINGLE:
        MST_linkage_core<T>(N, D_, Z2);
        break;
      case METHOD_METR_COMPLETE:
        NN_chain_core<METHOD_METR_COMPLETE, t_index, T>(N, D_, NULL, Z2);
        break;
      case METHOD_METR_AVERAGE:
        NN_chain_core<METHOD_METR_AVERAGE, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_WEIGHTED:
        NN_chain_core<METHOD_METR_WEIGHTED, t_index, T>(N, D_, NULL, Z2);
        break;
      case METHOD_METR_WARD:
        NN_chain_core<METHOD_METR_WARD, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_CENTROID:
        generic_linkage<METHOD_METR_CENTROID, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_MEDIAN:
        generic_linkage<METHOD_METR_MEDIAN, t_index, T>(N, D_, NULL, Z2);
        break;
      default:
        throw std::runtime_error(std::string("Invalid method index."));
      }

      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID ||
          method == METHOD_METR_MEDIAN) {
        Z2.sqrt();
      }

      T *const Z_ = reinterpret_cast<T *>(PyArray_DATA(Z));
      if (method == METHOD_METR_CENTROID || method == METHOD_METR_MEDIAN) {
        generate_SciPy_dendrogram<true, T>(Z_, Z2, N);
      } else {
        generate_SciPy_dendrogram<false, T>(Z_, Z2, N);
      }
    }

    else {
      using T = float;
      T *const D_ = reinterpret_cast<T *>(PyArray_DATA(D));
      cluster_result<T> Z2(N - 1);
      auto_array_ptr<t_index> members;
      // For these methods, the distance update formula needs the number of
      // data points in a cluster.
      if (method == METHOD_METR_AVERAGE || method == METHOD_METR_WARD ||
          method == METHOD_METR_CENTROID) {
        members.init(N, 1);
      }
      // Operate on squared distances for these methods.
      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID ||
          method == METHOD_METR_MEDIAN) {
        for (T *DD = D_;
             DD != D_ + static_cast<std::ptrdiff_t>(N) * (N - 1) / 2; ++DD)
          *DD *= *DD;
      }

      switch (method) {
      case METHOD_METR_SINGLE:
        MST_linkage_core<T>(N, D_, Z2);
        break;
      case METHOD_METR_COMPLETE:
        NN_chain_core<METHOD_METR_COMPLETE, t_index, T>(N, D_, NULL, Z2);
        break;
      case METHOD_METR_AVERAGE:
        NN_chain_core<METHOD_METR_AVERAGE, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_WEIGHTED:
        NN_chain_core<METHOD_METR_WEIGHTED, t_index, T>(N, D_, NULL, Z2);
        break;
      case METHOD_METR_WARD:
        NN_chain_core<METHOD_METR_WARD, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_CENTROID:
        generic_linkage<METHOD_METR_CENTROID, t_index, T>(N, D_, members, Z2);
        break;
      case METHOD_METR_MEDIAN:
        generic_linkage<METHOD_METR_MEDIAN, t_index, T>(N, D_, NULL, Z2);
        break;
      default:
        throw std::runtime_error(std::string("Invalid method index."));
      }

      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID ||
          method == METHOD_METR_MEDIAN) {
        Z2.sqrt();
      }

      T *const Z_ = reinterpret_cast<T *>(PyArray_DATA(Z));
      if (method == METHOD_METR_CENTROID || method == METHOD_METR_MEDIAN) {
        generate_SciPy_dendrogram<true, T>(Z_, Z2, N);
      } else {
        generate_SciPy_dendrogram<false, T>(Z_, Z2, N);
      }
    }
  } // try
  catch (const std::bad_alloc &) {
    return PyErr_NoMemory();
  } catch (const std::exception &e) {
    PyErr_SetString(PyExc_EnvironmentError, e.what());
    return NULL;
  } catch (const nan_error &) {
    PyErr_SetString(PyExc_FloatingPointError, "NaN dissimilarity value.");
    return NULL;
  }
#ifdef FE_INVALID
  catch (const fenv_error &) {
    PyErr_SetString(PyExc_FloatingPointError,
                    "NaN dissimilarity value in intermediate results.");
    return NULL;
  }
#endif
  catch (...) {
    PyErr_SetString(
        PyExc_EnvironmentError,
        "C++ exception (unknown reason). Please send a bug report.");
    return NULL;
  }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
  Py_RETURN_NONE;
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
}

/*
   Part 2: Clustering on vector data
*/

/* Metric codes.

   These codes must agree with the dictionary mtridx in fastcluster.py.
*/
enum metric_codes {
  // metrics
  METRIC_EUCLIDEAN = 0,
  METRIC_MINKOWSKI = 1,
  METRIC_CITYBLOCK = 2,
  METRIC_SEUCLIDEAN = 3,
  METRIC_SQEUCLIDEAN = 4,
  METRIC_COSINE = 5,
  METRIC_HAMMING = 6,
  METRIC_JACCARD = 7,
  METRIC_CHEBYCHEV = 8,
  METRIC_CANBERRA = 9,
  METRIC_BRAYCURTIS = 10,
  METRIC_MAHALANOBIS = 11,
  METRIC_YULE = 12,
  METRIC_MATCHING = 13,
  METRIC_DICE = 14,
  METRIC_ROGERSTANIMOTO = 15,
  METRIC_RUSSELLRAO = 16,
  METRIC_SOKALSNEATH = 17,
  METRIC_KULSINSKI = 18,
  METRIC_USER = 19,
  METRIC_INVALID = 20,      // sentinel
  METRIC_JACCARD_BOOL = 21, // separate function for Jaccard metric on
};                          // Boolean input data

/*
   Helper class: Throw this if calling the Python interpreter from within
   C returned an error.
*/
class pythonerror {};

/*
  This class handles all the information about the dissimilarity
  computation.
*/

template <typename T> class python_dissimilarity {
private:
  T *Xa;
  std::ptrdiff_t dim; // size_t saves many statis_cast<> in products
  t_index N;
  auto_array_ptr<T> Xnew;
  t_index *members;
  void (cluster_result<T>::*postprocessfn)(const T) const;
  T postprocessarg;

  T (python_dissimilarity::*distfn)(const t_index, const t_index) const;

  // for user-defined metrics
  PyObject *X_Python;
  PyObject *userfn;

  auto_array_ptr<T> precomputed;
  T *precomputed2;

  PyArrayObject *V;
  const T *V_data;

  // noncopyable
  python_dissimilarity();
  python_dissimilarity(python_dissimilarity const &);
  python_dissimilarity &operator=(python_dissimilarity const &);

public:
  // Ignore warning about uninitialized member variables. I know what I am
  // doing here, and some member variables are only used for certain metrics.
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
  python_dissimilarity(PyArrayObject *const Xarg, t_index *const members_,
                       const method_codes method, const metric_codes metric,
                       PyObject *const extraarg, bool temp_point_array)
      : Xa(reinterpret_cast<T *>(PyArray_DATA(Xarg))),
        dim(PyArray_DIM(Xarg, 1)),
        N(static_cast<t_index>(PyArray_DIM(Xarg, 0))),
        Xnew(temp_point_array ? (N - 1) * dim : 0), members(members_),
        postprocessfn(NULL), V(NULL) {
    int numpyType;
    if (typeid(T) == typeid(float)) {
      numpyType = NPY_FLOAT;
    } else if (typeid(T) == typeid(double)) {
      numpyType = NPY_DOUBLE;
    } else {
      // Handle error or unsupported type, should be caught earlier
      throw std::runtime_error(
          "Unsupported data type.");
    }
    switch (method) {
    case METHOD_METR_SINGLE:
      postprocessfn = NULL; // default
      switch (metric) {
      case METRIC_EUCLIDEAN:
        set_euclidean();
        break;
      case METRIC_SEUCLIDEAN:
        if (extraarg == NULL) {
          PyErr_SetString(
              PyExc_TypeError,
              "The 'seuclidean' metric needs a variance parameter.");
          throw pythonerror();
        }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
        V = reinterpret_cast<PyArrayObject *>(
            PyArray_FromAny(extraarg, PyArray_DescrFromType(numpyType), 1, 1,
                            NPY_ARRAY_CARRAY_RO, NULL));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
        if (PyErr_Occurred()) {
          throw pythonerror();
        }
        if (PyArray_DIM(V, 0) != dim) {
          PyErr_SetString(PyExc_ValueError, "The variance vector must have the "
                                            "same dimensionality as the data.");
          throw pythonerror();
        }
        V_data = reinterpret_cast<T *>(PyArray_DATA(V));
        distfn = &python_dissimilarity::seuclidean;
        postprocessfn = &cluster_result<T>::sqrt;
        break;
      case METRIC_SQEUCLIDEAN:
        distfn = &python_dissimilarity::sqeuclidean<false>;
        break;
      case METRIC_CITYBLOCK:
        set_cityblock();
        break;
      case METRIC_CHEBYCHEV:
        set_chebychev();
        break;
      case METRIC_MINKOWSKI:
        set_minkowski(extraarg);
        break;
      case METRIC_COSINE:
        distfn = &python_dissimilarity::cosine;
        postprocessfn = &cluster_result<T>::plusone;
        // precompute norms
        precomputed.init(N);
        for (t_index i = 0; i < N; ++i) {
          T sum = 0;
          for (t_index k = 0; k < dim; ++k) {
            sum += X(i, k) * X(i, k);
          }
          precomputed[i] = 1 / std::sqrt(sum);
        }
        break;
      case METRIC_HAMMING:
        distfn = &python_dissimilarity::hamming;
        postprocessfn = &cluster_result<T>::divide;
        postprocessarg = static_cast<T>(dim);
        break;
      case METRIC_JACCARD:
        distfn = &python_dissimilarity::jaccard;
        break;
      case METRIC_CANBERRA:
        distfn = &python_dissimilarity::canberra;
        break;
      case METRIC_BRAYCURTIS:
        distfn = &python_dissimilarity::braycurtis;
        break;
      case METRIC_MAHALANOBIS:
        if (extraarg == NULL) {
          PyErr_SetString(PyExc_TypeError,
                          "The 'mahalanobis' metric needs a parameter for the "
                          "inverse covariance.");
          throw pythonerror();
        }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
        V = reinterpret_cast<PyArrayObject *>(
            PyArray_FromAny(extraarg, PyArray_DescrFromType(numpyType), 2, 2,
                            NPY_ARRAY_CARRAY_RO, NULL));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
        if (PyErr_Occurred()) {
          throw pythonerror();
        }
        if (PyArray_DIM(V, 0) != N || PyArray_DIM(V, 1) != dim) {
          PyErr_SetString(PyExc_ValueError,
                          "The inverse covariance matrix has the wrong size.");
          throw pythonerror();
        }
        V_data = reinterpret_cast<T *>(PyArray_DATA(V));
        distfn = &python_dissimilarity::mahalanobis;
        postprocessfn = &cluster_result<T>::sqrt;
        break;
      case METRIC_YULE:
        distfn = &python_dissimilarity::yule;
        break;
      case METRIC_MATCHING:
        distfn = &python_dissimilarity::matching;
        postprocessfn = &cluster_result<T>::divide;
        postprocessarg = static_cast<T>(dim);
        break;
      case METRIC_DICE:
        distfn = &python_dissimilarity::dice;
        break;
      case METRIC_ROGERSTANIMOTO:
        distfn = &python_dissimilarity::rogerstanimoto;
        break;
      case METRIC_RUSSELLRAO:
        distfn = &python_dissimilarity::russellrao;
        postprocessfn = &cluster_result<T>::divide;
        postprocessarg = static_cast<T>(dim);
        break;
      case METRIC_SOKALSNEATH:
        distfn = &python_dissimilarity::sokalsneath;
        break;
      case METRIC_KULSINSKI:
        distfn = &python_dissimilarity::kulsinski;
        postprocessfn = &cluster_result<T>::plusone;
        precomputed.init(N);
        for (t_index i = 0; i < N; ++i) {
          t_index sum = 0;
          for (t_index k = 0; k < dim; ++k) {
            sum += Xb(i, k);
          }
          precomputed[i] = -.5 / static_cast<T>(sum);
        }
        break;
      case METRIC_USER:
        X_Python = reinterpret_cast<PyObject *>(Xarg);
        this->userfn = extraarg;
        distfn = &python_dissimilarity::user;
        break;
      default: // case METRIC_JACCARD_BOOL:
        distfn = &python_dissimilarity::jaccard_bool;
      }
      break;

    case METHOD_METR_WARD:
      postprocessfn = &cluster_result<T>::sqrtdouble;
      break;

    default:
      postprocessfn = &cluster_result<T>::sqrt;
    }
  }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif

  ~python_dissimilarity() {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
    Py_XDECREF(V);
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
  }

  inline T operator()(const t_index i, const t_index j) const {
    return (this->*distfn)(i, j);
  }

  inline T X(const t_index i, const t_index j) const { return Xa[i * dim + j]; }

  inline bool Xb(const t_index i, const t_index j) const {
    return reinterpret_cast<bool *>(Xa)[i * dim + j];
  }

  inline T *Xptr(const t_index i, const t_index j) const {
    return Xa + i * dim + j;
  }

  void merge(const t_index i, const t_index j, const t_index newnode) const {
    T const *const Pi = i < N ? Xa + i * dim : Xnew + (i - N) * dim;
    T const *const Pj = j < N ? Xa + j * dim : Xnew + (j - N) * dim;
    for (t_index k = 0; k < dim; ++k) {
      Xnew[(newnode - N) * dim + k] = (Pi[k] * static_cast<T>(members[i]) +
                                       Pj[k] * static_cast<T>(members[j])) /
                                      static_cast<T>(members[i] + members[j]);
    }
    members[newnode] = members[i] + members[j];
  }

  void merge_weighted(const t_index i, const t_index j,
                      const t_index newnode) const {
    T const *const Pi = i < N ? Xa + i * dim : Xnew + (i - N) * dim;
    T const *const Pj = j < N ? Xa + j * dim : Xnew + (j - N) * dim;
    for (t_index k = 0; k < dim; ++k) {
      Xnew[(newnode - N) * dim + k] = (Pi[k] + Pj[k]) * .5;
    }
  }

  void merge_inplace(const t_index i, const t_index j) const {
    T const *const Pi = Xa + i * dim;
    T *const Pj = Xa + j * dim;
    for (t_index k = 0; k < dim; ++k) {
      Pj[k] = (Pi[k] * static_cast<T>(members[i]) +
               Pj[k] * static_cast<T>(members[j])) /
              static_cast<T>(members[i] + members[j]);
    }
    members[j] += members[i];
  }

  void merge_inplace_weighted(const t_index i, const t_index j) const {
    T const *const Pi = Xa + i * dim;
    T *const Pj = Xa + j * dim;
    for (t_index k = 0; k < dim; ++k) {
      Pj[k] = (Pi[k] + Pj[k]) * .5;
    }
  }

  void postprocess(cluster_result<T> &Z2) const {
    if (postprocessfn != NULL) {
      (Z2.*postprocessfn)(postprocessarg);
    }
  }

  inline T ward(const t_index i, const t_index j) const {
    T mi = static_cast<T>(members[i]);
    T mj = static_cast<T>(members[j]);
    return sqeuclidean<true>(i, j) * mi * mj / (mi + mj);
  }

  inline T ward_initial(const t_index i, const t_index j) const {
    // alias for sqeuclidean
    // Factor 2!!!
    return sqeuclidean<true>(i, j);
  }

  // This method must not produce NaN if the input is non-NaN.
  inline static T ward_initial_conversion(const T min) { return min * .5; }

  inline T ward_extended(const t_index i, const t_index j) const {
    T mi = static_cast<T>(members[i]);
    T mj = static_cast<T>(members[j]);
    return sqeuclidean_extended(i, j) * mi * mj / (mi + mj);
  }

  /* We need two variants of the Euclidean metric: one that does not check
     for a NaN result, which is used for the initial distances, and one which
     does, for the updated distances during the clustering procedure.
  */
  template <const bool check_NaN>
  T sqeuclidean(const t_index i, const t_index j) const {
    T sum = 0;
    /*
      for (t_index k=0; k<dim; ++k) {
      T diff = X(i,k) - X(j,k);
      sum += diff*diff;
      }
    */
    // faster
    T const *Pi = Xa + i * dim;
    T const *Pj = Xa + j * dim;
    for (t_index k = 0; k < dim; ++k) {
      T diff = Pi[k] - Pj[k];
      sum += diff * diff;
    }
    if (check_NaN) {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
      if (fc_isnan(sum))
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
        throw(nan_error());
    }
    return sum;
  }

  T sqeuclidean_extended(const t_index i, const t_index j) const {
    T sum = 0;
    T const *Pi = i < N ? Xa + i * dim : Xnew + (i - N) * dim; // TBD
    T const *Pj = j < N ? Xa + j * dim : Xnew + (j - N) * dim;
    for (t_index k = 0; k < dim; ++k) {
      T diff = Pi[k] - Pj[k];
      sum += diff * diff;
    }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
    if (fc_isnan(sum))
      throw(nan_error());
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    return sum;
  }

private:
  void set_minkowski(PyObject *extraarg) {
    if (extraarg == NULL) {
      PyErr_SetString(PyExc_TypeError,
                      "The Minkowski metric needs a parameter.");
      throw pythonerror();
    }
    postprocessarg = static_cast<T>(PyFloat_AsDouble(extraarg));
    if (PyErr_Occurred()) {
      throw pythonerror();
    }

#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
    if (postprocessarg == std::numeric_limits<T>::infinity()) {
      set_chebychev();
    } else if (postprocessarg == 1.0) {
      set_cityblock();
    } else if (postprocessarg == 2.0) {
      set_euclidean();
    } else {
      distfn = &python_dissimilarity::minkowski;
      postprocessfn = &cluster_result<T>::power;
    }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
  }

  void set_euclidean() {
    distfn = &python_dissimilarity::sqeuclidean<false>;
    postprocessfn = &cluster_result<T>::sqrt;
  }

  void set_cityblock() { distfn = &python_dissimilarity::cityblock; }

  void set_chebychev() { distfn = &python_dissimilarity::chebychev; }

  T seuclidean(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      T diff = X(i, k) - X(j, k);
      sum += diff * diff / V_data[k];
    }
    return sum;
  }

  T cityblock(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      sum += std::abs(X(i, k) - X(j, k));
    }
    return sum;
  }

  T minkowski(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      sum += std::pow(std::abs(X(i, k) - X(j, k)), postprocessarg);
    }
    return sum;
  }

  T chebychev(const t_index i, const t_index j) const {
    T max = 0;
    for (t_index k = 0; k < dim; ++k) {
      T diff = std::abs(X(i, k) - X(j, k));
      if (diff > max) {
        max = diff;
      }
    }
    return max;
  }

  T cosine(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      sum -= X(i, k) * X(j, k);
    }
    return sum * precomputed[i] * precomputed[j];
  }

  T hamming(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
      sum += (X(i, k) != X(j, k));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    }
    return sum;
  }

  // Differs from scipy.spatial.distance: equal vectors correctly
  // return distance 0.
  T jaccard(const t_index i, const t_index j) const {
    t_index sum1 = 0;
    t_index sum2 = 0;
    for (t_index k = 0; k < dim; ++k) {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
      sum1 += (X(i, k) != X(j, k));
      sum2 += ((X(i, k) != 0) || (X(j, k) != 0));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    }
    return sum1 == 0 ? 0 : static_cast<T>(sum1) / static_cast<T>(sum2);
  }

  T canberra(const t_index i, const t_index j) const {
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      T numerator = std::abs(X(i, k) - X(j, k));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
      sum += numerator == 0
                 ? 0
                 : numerator / (std::abs(X(i, k)) + std::abs(X(j, k)));
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    }
    return sum;
  }

  T user(const t_index i, const t_index j) const {
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
    PyObject *u = PySequence_ITEM(X_Python, i);
    PyObject *v = PySequence_ITEM(X_Python, j);
    PyObject *result = PyObject_CallFunctionObjArgs(userfn, u, v, NULL);
    Py_DECREF(u);
    Py_DECREF(v);
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    if (result == NULL) {
      throw pythonerror();
    }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
    const T C_result = PyFloat_AsDouble(result);
    Py_DECREF(result);
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
    if (PyErr_Occurred()) {
      throw pythonerror();
    }
    return C_result;
  }

  T braycurtis(const t_index i, const t_index j) const {
    T sum1 = 0;
    T sum2 = 0;
    for (t_index k = 0; k < dim; ++k) {
      sum1 += std::abs(X(i, k) - X(j, k));
      sum2 += std::abs(X(i, k) + X(j, k));
    }
    return sum1 / sum2;
  }

  T mahalanobis(const t_index i, const t_index j) const {
    // V_data contains the product X*VI
    T sum = 0;
    for (t_index k = 0; k < dim; ++k) {
      sum += (V_data[i * dim + k] - V_data[j * dim + k]) * (X(i, k) - X(j, k));
    }
    return sum;
  }

  t_index mutable NTT; // 'local' variables
  t_index mutable NXO;
  t_index mutable NTF;
#define NTFFT NTF
#define NFFTT NTT

  void nbool_correspond(const t_index i, const t_index j) const {
    NTT = 0;
    NXO = 0;
    for (t_index k = 0; k < dim; ++k) {
      NTT += (Xb(i, k) & Xb(j, k));
      NXO += (Xb(i, k) ^ Xb(j, k));
    }
  }

  void nbool_correspond_tfft(const t_index i, const t_index j) const {
    NTT = 0;
    NXO = 0;
    NTF = 0;
    for (t_index k = 0; k < dim; ++k) {
      NTT += (Xb(i, k) & Xb(j, k));
      NXO += (Xb(i, k) ^ Xb(j, k));
      NTF += (Xb(i, k) & !Xb(j, k));
    }
    NTF *= (NXO - NTF);                             // NTFFT
    NTT *= (static_cast<t_index>(dim) - NTT - NXO); // NFFTT
  }

  void nbool_correspond_xo(const t_index i, const t_index j) const {
    NXO = 0;
    for (t_index k = 0; k < dim; ++k) {
      NXO += (Xb(i, k) ^ Xb(j, k));
    }
  }

  void nbool_correspond_tt(const t_index i, const t_index j) const {
    NTT = 0;
    for (t_index k = 0; k < dim; ++k) {
      NTT += (Xb(i, k) & Xb(j, k));
    }
  }

  T yule(const t_index i, const t_index j) const {
    nbool_correspond_tfft(i, j);
    return (NTFFT == 0)
               ? 0
               : static_cast<T>(2 * NTFFT) / static_cast<T>(NTFFT + NFFTT);
  }

  // Prevent a zero denominator for equal vectors.
  T dice(const t_index i, const t_index j) const {
    nbool_correspond(i, j);
    return (NXO == 0) ? 0 : static_cast<T>(NXO) / static_cast<T>(NXO + 2 * NTT);
  }

  T rogerstanimoto(const t_index i, const t_index j) const {
    nbool_correspond_xo(i, j);
    return static_cast<T>(2 * NXO) / static_cast<T>(NXO + dim);
  }

  T russellrao(const t_index i, const t_index j) const {
    nbool_correspond_tt(i, j);
    return static_cast<T>(dim - NTT);
  }

  // Prevent a zero denominator for equal vectors.
  T sokalsneath(const t_index i, const t_index j) const {
    nbool_correspond(i, j);
    return (NXO == 0) ? 0
                      : static_cast<T>(2 * NXO) / static_cast<T>(NTT + 2 * NXO);
  }

  T kulsinski(const t_index i, const t_index j) const {
    nbool_correspond_tt(i, j);
    return static_cast<T>(NTT) * (precomputed[i] + precomputed[j]);
  }

  // 'matching' distance = Hamming distance
  T matching(const t_index i, const t_index j) const {
    nbool_correspond_xo(i, j);
    return static_cast<T>(NXO);
  }

  // Prevent a zero denominator for equal vectors.
  T jaccard_bool(const t_index i, const t_index j) const {
    nbool_correspond(i, j);
    return (NXO == 0) ? 0 : static_cast<T>(NXO) / static_cast<T>(NXO + NTT);
  }
};

static PyObject *linkage_vector_wrap(PyObject *const, PyObject *const args) {
  PyArrayObject *X, *Z;
  unsigned char method, metric;
  PyObject *extraarg;
  int high_precision_flag = 0;

  try {
    // Parse the input arguments
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
    if (!PyArg_ParseTuple(args, "O!O!bbOp", &PyArray_Type, &X, // NumPy array
                          &PyArray_Type, &Z,                   // NumPy array
                          &method,                             // unsigned char
                          &metric,                             // unsigned char
                          &extraarg,                           // PyObject
                          &high_precision_flag)) {             // int
      return NULL;
    }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif

    bool high_precision = high_precision_flag != 0;
    int T_mant_dig = high_precision ? DBL_MANT_DIG : FLT_MANT_DIG;

    if (PyArray_NDIM(X) != 2) {
      PyErr_SetString(PyExc_ValueError,
                      "The input array must be two-dimensional.");
    }
    npy_intp const N_ = PyArray_DIM(X, 0);
    if (N_ < 1) {
      // N must be at least 1.
      PyErr_SetString(PyExc_ValueError,
                      "At least one element is needed for clustering.");
      return NULL;
    }

    npy_intp const dim = PyArray_DIM(X, 1);
    if (dim < 1) {
      PyErr_SetString(PyExc_ValueError, "Invalid dimension of the data set.");
      return NULL;
    }

    /*
      (1)
      The biggest index used below is 4*(N-2)+3, as an index to Z. This must
      fit into the data type used for indices.
      (2)
      The largest representable integer, without loss of precision, by a
      floating point number of type t_float is 2^T_FLOAT_MANT_DIG. Here, we
      make sure that all cluster labels from 0 to 2N-2 in the output can be
      accurately represented by a floating point number.

      Conversion of N to 64 bits below is not really necessary but it prevents
      a warning ("shift count >= width of type") on systems where "int" is 32
      bits wide.
    */
    if (N_ > MAX_INDEX / 4 || dim > MAX_INDEX ||
        static_cast<int64_t>(N_ - 1) >> (T_mant_dig - 1) > 0) {
      PyErr_SetString(PyExc_ValueError, "Data is too big, index overflow.");
      return NULL;
    }
    t_index N = static_cast<t_index>(N_);

    if (high_precision) {
      using T = double;
      cluster_result<T> Z2(N - 1);

      auto_array_ptr<t_index> members;
      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID) {
        members.init(2 * N - 1, 1);
      }

      if ((method != METHOD_METR_SINGLE && metric != METRIC_EUCLIDEAN) ||
          metric >= METRIC_INVALID) {
        PyErr_SetString(PyExc_IndexError, "Invalid metric index.");
        return NULL;
      }

      if (PyArray_ISBOOL(X)) {
        if (metric == METRIC_HAMMING) {
          metric = METRIC_MATCHING; // Alias
        }
        if (metric == METRIC_JACCARD) {
          metric = METRIC_JACCARD_BOOL;
        }
      }

      if (extraarg != Py_None && metric != METRIC_MINKOWSKI &&
          metric != METRIC_SEUCLIDEAN && metric != METRIC_MAHALANOBIS &&
          metric != METRIC_USER) {
        PyErr_SetString(PyExc_TypeError,
                        "No extra parameter is allowed for this metric.");
        return NULL;
      }

      /* temp_point_array must be true if the alternative algorithm
         is used below (currently for the centroid and median methods). */
      bool temp_point_array =
          (method == METHOD_METR_CENTROID || method == METHOD_METR_MEDIAN);

      python_dissimilarity<T> dist(
          X, members, static_cast<method_codes>(method),
          static_cast<metric_codes>(metric), extraarg, temp_point_array);

      if (method != METHOD_METR_SINGLE && method != METHOD_METR_WARD &&
          method != METHOD_METR_CENTROID && method != METHOD_METR_MEDIAN) {
        PyErr_SetString(PyExc_IndexError, "Invalid method index.");
        return NULL;
      }

      // Allow threads if the metric is not "user"!
      GIL_release G(metric != METRIC_USER);

      switch (method) {
      case METHOD_METR_SINGLE:
        MST_linkage_core_vector<T>(N, dist, Z2);
        break;
      case METHOD_METR_WARD:
        generic_linkage_vector<METHOD_VECTOR_WARD, T>(N, dist, Z2);
        break;
      case METHOD_METR_CENTROID:
        generic_linkage_vector_alternative<METHOD_VECTOR_CENTROID, T>(N, dist, Z2);
        break;
      default: // case METHOD_METR_MEDIAN:
        generic_linkage_vector_alternative<METHOD_VECTOR_MEDIAN, T>(N, dist, Z2);
      }

      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID) {
        members.free();
      }

      dist.postprocess(Z2);

      T *const Z_ = reinterpret_cast<T *>(PyArray_DATA(Z));
      if (method != METHOD_METR_SINGLE) {
        generate_SciPy_dendrogram<true, T>(Z_, Z2, N);
      } else {
        generate_SciPy_dendrogram<false, T>(Z_, Z2, N);
      }
    } else {
      using T = float;
      cluster_result<T> Z2(N - 1);

      auto_array_ptr<t_index> members;
      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID) {
        members.init(2 * N - 1, 1);
      }

      if ((method != METHOD_METR_SINGLE && metric != METRIC_EUCLIDEAN) ||
          metric >= METRIC_INVALID) {
        PyErr_SetString(PyExc_IndexError, "Invalid metric index.");
        return NULL;
      }

      if (PyArray_ISBOOL(X)) {
        if (metric == METRIC_HAMMING) {
          metric = METRIC_MATCHING; // Alias
        }
        if (metric == METRIC_JACCARD) {
          metric = METRIC_JACCARD_BOOL;
        }
      }

      if (extraarg != Py_None && metric != METRIC_MINKOWSKI &&
          metric != METRIC_SEUCLIDEAN && metric != METRIC_MAHALANOBIS &&
          metric != METRIC_USER) {
        PyErr_SetString(PyExc_TypeError,
                        "No extra parameter is allowed for this metric.");
        return NULL;
      }

      /* temp_point_array must be true if the alternative algorithm
         is used below (currently for the centroid and median methods). */
      bool temp_point_array =
          (method == METHOD_METR_CENTROID || method == METHOD_METR_MEDIAN);

      python_dissimilarity<T> dist(
          X, members, static_cast<method_codes>(method),
          static_cast<metric_codes>(metric), extraarg, temp_point_array);

      if (method != METHOD_METR_SINGLE && method != METHOD_METR_WARD &&
          method != METHOD_METR_CENTROID && method != METHOD_METR_MEDIAN) {
        PyErr_SetString(PyExc_IndexError, "Invalid method index.");
        return NULL;
      }

      // Allow threads if the metric is not "user"!
      GIL_release G(metric != METRIC_USER);

      switch (method) {
      case METHOD_METR_SINGLE:
        MST_linkage_core_vector<T>(N, dist, Z2);
        break;
      case METHOD_METR_WARD:
        generic_linkage_vector<METHOD_VECTOR_WARD, T>(N, dist, Z2);
        break;
      case METHOD_METR_CENTROID:
        generic_linkage_vector_alternative<METHOD_VECTOR_CENTROID, T>(N, dist, Z2);
        break;
      default: // case METHOD_METR_MEDIAN:
        generic_linkage_vector_alternative<METHOD_VECTOR_MEDIAN, T>(N, dist, Z2);
      }

      if (method == METHOD_METR_WARD || method == METHOD_METR_CENTROID) {
        members.free();
      }

      dist.postprocess(Z2);

      T *const Z_ = reinterpret_cast<T *>(PyArray_DATA(Z));
      if (method != METHOD_METR_SINGLE) {
        generate_SciPy_dendrogram<true, T>(Z_, Z2, N);
      } else {
        generate_SciPy_dendrogram<false, T>(Z_, Z2, N);
      }
    }
  } // try
  catch (const std::bad_alloc &) {
    return PyErr_NoMemory();
  } catch (const std::exception &e) {
    PyErr_SetString(PyExc_EnvironmentError, e.what());
    return NULL;
  } catch (const nan_error &) {
    PyErr_SetString(PyExc_FloatingPointError, "NaN dissimilarity value.");
    return NULL;
  } catch (const pythonerror) {
    return NULL;
  } catch (...) {
    PyErr_SetString(
        PyExc_EnvironmentError,
        "C++ exception (unknown reason). Please send a bug report.");
    return NULL;
  }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
  Py_RETURN_NONE;
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
}

#if HAVE_VISIBILITY
#pragma GCC visibility pop
#endif