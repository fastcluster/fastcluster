#!/usr/bin/env bash

set -e -x

arch=$(uname -m)

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    if [[ "${PYBIN}" == *"cp36"* && $arch != "aarch64" ]] || \
       [[ "${PYBIN}" == *"cp37"* ]] || \
       [[ "${PYBIN}" == *"cp38"* ]] || \
       [[ "${PYBIN}" == *"cp39"* ]]; then
        cd /io
        "${PYBIN}/pip" install --upgrade pip wheel
        "${PYBIN}/pip" install .[test]
        "${PYBIN}/python" setup.py bdist_wheel
        "${PYBIN}/python" setup.py test

        #"${PYBIN}/pip" install numpy scipy
        #"${PYBIN}/pip" install -e /io/
        #(cd /io && "${PYBIN}/python" setup.py test)
        #"${PYBIN}/pip" wheel /io/ -w dist/
        #rm -rf /io/build /io/*.egg-info
    fi
done

pwd
ls -la dist
ls -la /io
ls -la /io/dist

# Bundle external shared libraries into the wheels
for whl in dist/fastcluster*.whl; do
    auditwheel repair "$whl" -w /io/dist/
done
