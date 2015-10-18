#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test whether the fastcluster package correctly recognizes NaN values
and raises a FloatingPointError.'''
import sys
if sys.hexversion < 0x03000000: # uniform unicode handling for both Python 2.x and 3.x
    def u(x):
        return x.decode('utf-8')
else:
    def u(x):
        return x
print(u('''Test program for the 'fastcluster' package.

Copyright (c) 2011 Daniel Müllner, <http://danifold.net>
''').encode('utf-8'))
import numpy as np
import fastcluster

version = '1.1.20'
if fastcluster.__version__ != version:
    raise ValueError('Wrong module version: {} instead of {}.'.format(fastcluster.__version__, version))

import atexit
def print_seed():
  print("Seed: {0}".format(seed))
atexit.register(print_seed)

seed = np.random.randint(0,1e9)

np.random.seed(seed)

def test():
    n = np.random.random_integers(2,100)

    # Part 1: distance matrix input

    N = n*(n-1)//2
    D = np.random.rand(N)
    # Insert a single NaN value
    pos = np.random.randint(N)
    D[pos] = np.nan

    for method in ['single', 'complete', 'average', 'weighted', 'ward',
                   'centroid', 'median']:
        try:
            fastcluster.linkage(D, method=method)
            raise AssertionError('fastcluster did not detect a NaN value!')
        except FloatingPointError:
            pass

    # Next: the original array does not contain a NaN, but a NaN occurs
    # as an updated distance.
    for method in ['average', 'weighted', 'ward', 'centroid', 'median']:
        try:
            fastcluster.linkage([np.inf,-np.inf,-np.inf], method=method)
            raise AssertionError('fastcluster did not detect a NaN value!')
        except FloatingPointError:
            pass

    # Part 2: vector input

    dim = np.random.random_integers(2,12)
    X = np.random.rand(n,dim)
    pos = (np.random.randint(n), np.random.randint(dim))
    # Insert a single NaN coordinate
    X[pos] = np.nan

    for method in ['single', 'ward', 'centroid', 'median']:
        try:
            fastcluster.linkage_vector(X, method=method)
            raise AssertionError('fastcluster did not detect a NaN value!')
        except FloatingPointError:
            pass

    return True

if __name__ == "__main__":
    test()
    print('OK.')
