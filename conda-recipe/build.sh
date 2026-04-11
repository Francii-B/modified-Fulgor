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

echo "Build modified-Fulgor"
mkdir -p build
cd build
cmake .. ${CMAKE_ARGS}
make -j${CPU_COUNT}

mkdir -p "${PREFIX}/bin"
cp fulgor "${PREFIX}/bin/"
