#!/bin/bash
set -euo pipefail

uv lock --check

MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-11.0}"
MACOS_CMAKE_ARGS="${CMAKE_ARGS:-} -DCMAKE_OSX_DEPLOYMENT_TARGET=$MACOSX_DEPLOYMENT_TARGET"
export MACOSX_DEPLOYMENT_TARGET

rm -rf .venv
for PYTHON_VERSION in 3.9 3.10 3.11 3.12 3.13 3.14; do
  MACOSX_DEPLOYMENT_TARGET="$MACOSX_DEPLOYMENT_TARGET" \
    CMAKE_ARGS="$MACOS_CMAKE_ARGS" \
    uv run --python "$PYTHON_VERSION" --with build python -m build
done
./repair-macos-wheels.sh

rm -rf .venv
for PYTHON_VERSION in 3.9 3.10 3.11 3.12 3.13 3.14; do
  orb -m ubuntu-plucky uv run --python "$PYTHON_VERSION" --with build python -m build
done

orb -m ubuntu-plucky ./repair-wheels.sh
