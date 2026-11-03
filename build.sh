#!/bin/bash

set -e

echo "==> Creating build directory ..."
mkdir -p build
cd build

echo "==> Setting the project with CMake ..."
cmake -DCMAKE_BUILD_TYPE=Release ..

echo "==> Compiling ..."
cmake --build . -- -j$(nproc)

echo "==> Build finished"
echo "==> The executable is in build/sos"
