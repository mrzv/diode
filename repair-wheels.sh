#!/bin/bash

DIODE_VERSION=$(python3 -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")

pushd dist
for fn in diode-$DIODE_VERSION-*-linux*.whl; do
  uvx auditwheel repair --plat manylinux_2_39_x86_64 $fn
done
popd
