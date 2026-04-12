#!/bin/bash
set -euxo pipefail

echo "Sanity check submodules"
test -d external/ggcat
test -d external/sshash

echo "Pin Rust time crate to avoid E0282 with newer rustc"
pushd external/ggcat/libs-crates/dynamic-dispatch-rs

python - <<'PY'
from pathlib import Path
p = Path("Cargo.toml")
s = p.read_text()
needle = "[dependencies]\n"
entry = 'time = "=0.3.37"\n'
if 'time = ' not in s:
    s = s.replace(needle, needle + entry, 1)
p.write_text(s)
PY

cargo update --manifest-path Cargo.toml -p time --precise 0.3.37
popd

echo "Build GGCAT"
pushd external/ggcat/crates/capi/ggcat-cpp-api

make

# Fail early if the Makefile did not put the library where Fulgor expects it
mkdir -p lib
test -f lib/libggcat_cpp_bindings.a || cp ../../../target/release/libggcat_cpp_bindings.a lib/

test -f lib/libggcat_cpp_bindings.a
test -f lib/libggcat_api.a

popd

echo "Build modified-Fulgor"
mkdir -p "${SRC_DIR}/build"
pushd "${SRC_DIR}/build"

cmake .. ${CMAKE_ARGS}
make -j${CPU_COUNT}

mkdir -p "${PREFIX}/bin"
cp fulgor "${PREFIX}/bin/"

popd
