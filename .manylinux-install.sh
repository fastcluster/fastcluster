#!/usr/bin/env bash

set -e -x

ls -l /opt/python

cd /io

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    if [[ "${PYBIN}" != *"cp3"* ]] || \
       [[ "${PYBIN}" == *"cp35"* ]] || \
       [[ "${PYBIN}" == *"cp36"* ]] || \
       [[ "${PYBIN}" == *"cp37"* ]]; then
        continue
    fi

    "${PYBIN}/pip" install --upgrade pip wheel
    "${PYBIN}/pip" install oldest_supported_numpy
    "${PYBIN}/pip" install .[test]
    "${PYBIN}/python" setup.py bdist_wheel
    "${PYBIN}/python" setup.py test
done

# Bundle external shared libraries into the wheels
for whl in dist/fastcluster*.whl; do
    auditwheel repair "$whl"
done

rm dist/*
mv wheelhouse/*.whl dist/
