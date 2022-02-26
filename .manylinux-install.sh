#!/usr/bin/env bash

set -e -x

cd /io

ls -l /opt/python

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    if [[ "${PYBIN}" != *"cp38"* ]] && \
       [[ "${PYBIN}" != *"cp39"* ]] && \
       [[ "${PYBIN}" != *"cp310"* ]]; then
        continue
    fi

    "${PYBIN}/python" -m venv fastcluster-env
    fastcluster-env/bin/activate
    pip install --upgrade pip wheel
    pip install oldest_supported_numpy
    pip install .[test]
    python setup.py bdist_wheel
    python setup.py test
done

# Bundle external shared libraries into the wheels
for whl in dist/fastcluster*.whl; do
    auditwheel repair "$whl"
done

rm dist/*
mv wheelhouse/*.whl dist/
