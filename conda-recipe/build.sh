#!/bin/bash
set -euxo pipefail

echo "Sanity check submodules"
test -d external/ggcat
test -d external/sshash

echo "Build GGCAT"
(
  cd external/ggcat/crates/capi/ggcat-cpp-api
  make
)

LIB_PATH=$(find . -name "libggcat_cpp_bindings.a" | head -n 1)
if [ -z "${LIB_PATH}" ]; then
  echo "ERROR: libggcat_cpp_bindings.a not found after cargo build"
  exit 1
fi

cp "${LIB_PATH}" external/ggcat/crates/capi/ggcat-cpp-api/lib/

echo "Build modified-Fulgor"
mkdir -p build
cd build
cmake .. ${CMAKE_ARGS}
make -j${CPU_COUNT}

mkdir -p "${PREFIX}/bin"
cp fulgor "${PREFIX}/bin/"
