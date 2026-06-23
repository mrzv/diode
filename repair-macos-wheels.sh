#!/bin/bash
set -euo pipefail
shopt -s nullglob

DIODE_VERSION=$(python3 -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")
MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-15.0}"
export MACOSX_DEPLOYMENT_TARGET

REPAIRED_DIR=$(mktemp -d "${TMPDIR:-/tmp}/diode-delocate.XXXXXX")
trap 'rm -rf "$REPAIRED_DIR"' EXIT
UNREPAIRED_DIR=dist/unrepaired-macos
mkdir -p "$UNREPAIRED_DIR"

WHEELS=(dist/diode-$DIODE_VERSION-*-macosx*.whl)
RAW_WHEELS=()
if (( ${#WHEELS[@]} > 0 )); then
  for fn in "${WHEELS[@]}"; do
    raw_fn="$UNREPAIRED_DIR/$(basename "$fn")"
    mv "$fn" "$raw_fn"
    RAW_WHEELS+=("$raw_fn")
  done
else
  RAW_WHEELS=("$UNREPAIRED_DIR"/diode-$DIODE_VERSION-*-macosx*.whl)
fi

if (( ${#RAW_WHEELS[@]} == 0 )); then
  echo "No macOS wheels found for diode $DIODE_VERSION" >&2
  exit 1
fi

for fn in "${RAW_WHEELS[@]}"; do
  uvx --from delocate delocate-wheel \
    --wheel-dir "$REPAIRED_DIR" \
    "$fn"
done

REPAIRED_WHEELS=("$REPAIRED_DIR"/*.whl)
if (( ${#REPAIRED_WHEELS[@]} == 0 )); then
  echo "delocate did not produce any repaired wheels" >&2
  exit 1
fi

for fn in "${REPAIRED_WHEELS[@]}"; do
  mv "$fn" dist/
done
