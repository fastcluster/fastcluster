#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TBD test single on integer matrices for hamming/jaccard
print('''
Test program for the 'fastcluster' package.
Copyright:
  * Until package version 1.1.23: (c) 2011 Daniel Müllner <http://danifold.net>
  * All changes from version 1.1.24 on: (c) Google Inc. <http://google.com>''')
import sys
import fastcluster as fc
import numpy as np
from scipy.spatial.distance import pdist, squareform
import math

version = '1.2.4'
if fc.__version__ != version:
    raise ValueError('Wrong module version: {} instead of {}.'.format(fc.__version__, version))

import atexit
def print_seed():
  print("Seed: {0}".format(seed))
atexit.register(print_seed)

seed  = np.random.randint(0,1e9)

print_seed()
np.random.seed(seed)
abstol = 1e-14 # absolute tolerance
rtol = 1e-13 # relative tolerance

# NaN values are used in computations. Do not warn about them.
np.seterr(invalid='ignore')

def correct_for_zero_vectors(D, pcd, metric):
    # Correct some metrics: we want the distance from the zero vector
    # to itself to be 0, not NaN.
    if metric in ('jaccard', 'dice', 'sokalsneath'):
        z = np.flatnonzero(np.all(pcd==0, axis=1))
        if len(z):
            DD = squareform(D)
            DD[np.ix_(z, z)] = 0
            D = squareform(DD)
    return D

def test_all(n,dim):
  method = 'single'

  # metrics for boolean vectors
  pcd = np.random.randint(0, 2, size=(n,dim), dtype=bool)
  pcd2 = pcd.copy()

  for metric in ('hamming', 'jaccard', 'yule', 'matching', 'dice',
                 'rogerstanimoto',
                 #'sokalmichener',
                 # exclude, bug in Scipy
                 # http://projects.scipy.org/scipy/ticket/1486
                 'russellrao', 'sokalsneath',
                 #'kulsinski'
                 # exclude, bug in Scipy
                 # http://projects.scipy.org/scipy/ticket/1484
                 ):
    sys.stdout.write("Metric: " + metric + "...")
    D = pdist(pcd, metric=metric)
    D = correct_for_zero_vectors(D, pcd, metric)

    try:
        Z2 = fc.linkage_vector(pcd, method, metric)
    except FloatingPointError:
        # If linkage_vector reported a NaN dissimilarity value,
        # check whether the distance matrix really contains NaN.
        if np.any(np.isnan(D)):
            print("Skip this test: NaN dissimilarity value.")
            continue
        else:
            raise AssertionError('"linkage_vector" erroneously reported NaN.')

    if np.any(pcd2!=pcd):
      raise AssertionError('Input array was corrupted.', pcd)
    check(Z2, method, D)

  # metrics for real vectors
  bound = math.sqrt(n)
  pcd = np.random.randint(-bound, bound + 1, (n,dim))
  for metric in ['euclidean', 'sqeuclidean', 'cityblock', 'chebychev',
                 'minkowski', 'cosine', 'correlation', 'hamming', 'jaccard',
                 'canberra',
                 # canberra: see bug in older Scipy versions
                 # http://projects.scipy.org/scipy/ticket/1430
                 'braycurtis', 'seuclidean', 'mahalanobis', 'user']:
    sys.stdout.write("Metric: " + metric + "...")
    if metric=='minkowski':
        p = np.random.uniform(1.,10.)
        sys.stdout.write("p: " + str(p) + "...")
        D = pdist(pcd, metric=metric, p=p)
        Z2 = fc.linkage_vector(pcd, method, metric, p)
    elif metric=='user':
        # Euclidean metric as a user function
        fn = (lambda u, v: np.sqrt(((u-v)*(u-v).T).sum()))
        D = pdist(pcd, metric=fn)
        Z2 = fc.linkage_vector(pcd, method, fn)
    else:
        D = pdist(pcd, metric=metric)
        D = correct_for_zero_vectors(D, pcd, metric)
        try:
            Z2 = fc.linkage_vector(pcd, method, metric)
        except FloatingPointError:
            if np.any(np.isnan(D)):
                print("Skip this test: NaN dissimilarity value.")
                continue
            else:
                raise AssertionError(
                    '"linkage_vector" erroneously reported NaN.')

    check(Z2, method, D)

  D = pdist(pcd)
  for method in ['ward', 'centroid', 'median']:
    Z2 = fc.linkage_vector(pcd, method)
    check(Z2, method, D)

def check(Z2, method, D):
    sys.stdout.write("Method: " + method + "...")
    I = np.array(Z2[:,:2], dtype=int)

    Ds = squareform(D)
    n = len(Ds)
    row_repr = np.arange(2*n-1)
    row_repr[n:] = -1
    size = np.ones(n, dtype=int)

    np.fill_diagonal(Ds, np.nan)

    mins = np.empty(n-1)

    for i in range(n-1):
      for j in range(n-1):
        mins[j] = np.nanmin(Ds[j,j+1:])
      gmin = np.nanmin(mins)
      if abs(Z2[i,2]-gmin) > max(abs(Z2[i,2]),abs(gmin))*rtol and \
            abs(Z2[i,2]-gmin)>abstol:
          raise AssertionError(
              'Not the global minimum in step {2}: {0}, {1}'.
              format(Z2[i,2], gmin,i), squareform(D))
      i1, i2 = row_repr[I[i,:]]
      if (i1<0):
        raise AssertionError('Negative index i1.', squareform(D))
      if (i2<0):
        raise AssertionError('Negative index i2.', squareform(D))
      if I[i,0]>=I[i,1]:
        raise AssertionError('Convention violated.', squareform(D))
      if i1>i2:
        i1, i2 = i2, i1
      if abs(Ds[i1,i2]-gmin) > max(abs(Ds[i1,i2]),abs(gmin))*rtol and \
            abs(Ds[i1,i2]-gmin)>abstol:
          raise AssertionError(
              'The global minimum is not at the right place in step {5}: '
              '({0}, {1}): {2} != {3}. Difference: {4}'
              .format(i1, i2, Ds[i1, i2], gmin, Ds[i1, i2]-gmin, i),
              squareform(D))

      s1 = size[i1]
      s2 = size[i2]
      S = float(s1+s2)
      if method=='single':
          if i1>0: # mostly unnecessary; workaround for a bug/feature in NumPy
          # 1.7.0.dev, see http://projects.scipy.org/numpy/ticket/2078
              Ds[:i1,i2]   = np.min( Ds[:i1,(i1,i2)],axis=1)
          Ds[i1:i2,i2] = np.minimum(Ds[i1,i1:i2],Ds[i1:i2,i2])
          Ds[i2,i2:]   = np.min( Ds[(i1,i2),i2:],axis=0)
      elif method=='complete':
          if i1>0:
              Ds[:i1,i2]   = np.max( Ds[:i1,(i1,i2)],axis=1)
          Ds[i1:i2,i2] = np.maximum(Ds[i1,i1:i2],Ds[i1:i2,i2])
          Ds[i2,i2:]   = np.max( Ds[(i1,i2),i2:],axis=0)
      elif method=='average':
          Ds[:i1,i2]   = ( Ds[:i1,i1]*s1 + Ds[:i1,i2]*s2 ) / S
          Ds[i1:i2,i2] = ( Ds[i1,i1:i2]*s1 + Ds[i1:i2,i2]*s2 ) / S
          Ds[i2,i2:]   = ( Ds[i1,i2:]*s1 + Ds[i2,i2:]*s2 ) / S
      elif method=='weighted':
          if i1>0:
              Ds[:i1,i2]   = np.mean( Ds[:i1,(i1,i2)],axis=1)
          Ds[i1:i2,i2] = ( Ds[i1,i1:i2] + Ds[i1:i2,i2] )*.5
          Ds[i2,i2:]   = np.mean( Ds[(i1,i2),i2:],axis=0)
      elif method=='ward':
          Ds[:i1,i2]   = np.sqrt((np.square(Ds[:i1,i1])*(s1+size[:i1])
                         -gmin*gmin*size[:i1]
                         +np.square(Ds[:i1,i2])*(s2+size[:i1]))/(S+size[:i1]))
          Ds[i1:i2,i2] = np.sqrt((np.square(Ds[i1,i1:i2])*(s1+size[i1:i2])
                         -gmin*gmin*size[i1:i2]
                         +np.square(Ds[i1:i2,i2])*(s2+size[i1:i2]))
                                 /(S+size[i1:i2]))
          Ds[i2,i2:]   = np.sqrt((np.square(Ds[i1,i2:])*(s1+size[i2:])
                         -gmin*gmin*size[i2:]
                         +np.square(Ds[i2,i2:])*(s2+size[i2:]))/(S+size[i2:]))
      elif method=='centroid':
          Ds[:i1,i2]   = np.sqrt((np.square(Ds[:i1,i1])*s1
                         +np.square(Ds[:i1,i2])*s2)*S-gmin*gmin*s1*s2) / S
          Ds[i1:i2,i2] = np.sqrt((np.square(Ds[i1,i1:i2])*s1
                         +np.square(Ds[i1:i2,i2])*s2)*S-gmin*gmin*s1*s2) / S
          Ds[i2,i2:]   = np.sqrt((np.square(Ds[i1,i2:])*s1
                         +np.square(Ds[i2,i2:])*s2)*S-gmin*gmin*s1*s2) / S
      elif method=='median':
          Ds[:i1,i2]   = np.sqrt((np.square(Ds[:i1,i1])
                                  +np.square(Ds[:i1,i2]))*2-gmin*gmin)*.5
          Ds[i1:i2,i2] = np.sqrt((np.square(Ds[i1,i1:i2])
                                  +np.square(Ds[i1:i2,i2]))*2-gmin*gmin)*.5
          Ds[i2,i2:]   = np.sqrt((np.square(Ds[i1,i2:])
                                  +np.square(Ds[i2,i2:]))*2-gmin*gmin)*.5
      else:
          raise ValueError('Unknown method.')

      Ds[i1, i1:n] = np.inf
      Ds[:i1, i1] = np.inf
      row_repr[n+i] = i2
      size[i2] = S
    print('OK.')

def test(repeats):
    if repeats:
        iterator = range(repeats)
    else:
        import itertools
        iterator = itertools.repeat(None)
        print('''
If everything is OK, the test program will run forever, without an error
message.
''')
    for _ in iterator:
        dim = np.random.randint(2, 13)
        n = np.random.randint(max(2*dim,5),200)

        print('Dimension: {0}'.format(dim))
        print('Number of points: {0}'.format(n))

        test_all(n,dim)

if __name__ == "__main__":
    test(None)
