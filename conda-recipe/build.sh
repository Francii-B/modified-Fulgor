#!/bin/bash
set -euxo pipefail

echo "Sanity check submodules"
test -d external/ggcat
test -d external/sshash

echo "Pin Rust time crate to avoid E0282 with newer rustc"
cargo update --manifest-path external/ggcat/Cargo.toml -p time --precise 0.3.37

echo "Build GGCAT"
pushd external/ggcat/crates/capi/ggcat-cpp-api
make
popd

LIB_PATH=$(find external/ggcat -name "libggcat_cpp_bindings.a" | head -n 1)
if [ -z "${LIB_PATH}" ]; then
echo "ERROR: libggcat_cpp_bindings.a not found after cargo build!"
find external/ggcat -maxdepth 6 \( -name target -o -name lib \) -print
exit 1
fi

echo "Found bindings at: ${LIB_PATH}"
mkdir -p external/ggcat/crates/capi/ggcat-cpp-api/lib
cp "${LIB_PATH}" external/ggcat/crates/capi/ggcat-cpp-api/lib/

test -f external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cpp_bindings.a
test -f external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_api.a

echo "Build modified-Fulgor"
mkdir -p "${SRC_DIR}/build"
pushd "${SRC_DIR}/build"
cmake .. ${CMAKE_ARGS}
make -j${CPU_COUNT}
mkdir -p "${PREFIX}/bin"
cp fulgor "${PREFIX}/bin/"
popd
