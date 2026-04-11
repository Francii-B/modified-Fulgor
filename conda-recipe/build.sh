#!/bin/bash
set -euxo pipefail

echo "Sanity check submodules"
test -d external/ggcat
test -d external/sshash


echo "Pin Rust time crate to avoid E0282 with newer rustc"
cd external/ggcat/libs-crates/dynamic-dispatch-rs

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

cd ../../crates/capi/ggcat-cpp-api
make

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
