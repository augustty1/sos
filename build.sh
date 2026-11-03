#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build}"
BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"

if command -v nproc >/dev/null 2>&1; then
  DEFAULT_JOBS="$(nproc)"
else
  DEFAULT_JOBS=2
fi
JOBS="${JOBS:-${DEFAULT_JOBS}}"

case "${1:-}" in
  --clean)
    echo "==> Removing build directory. ${BUILD_DIR}"
    rm -rf -- "${BUILD_DIR}"
    ;;
  -h|--help)
    sed -n '2,8p' "${BASH_SOURCE[0]}"
    exit 0
    ;;
  "") ;;
  *)
    echo "Uso: $0 [--clean|--help]" >&2
    exit 2
    ;;
esac

if [[ ! -f "${ROOT_DIR}/CMakeLists.txt" ]]; then
  echo "Error. CMakeLists.txt not found. ${ROOT_DIR}" >&2
  exit 1
fi

if ! command -v cmake >/dev/null 2>&1; then
  echo "Error. CMake not found." >&2
  exit 1
fi

CMAKE_ARGS=(-S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}")

if [[ -n "${SCIP_DIR:-}" ]]; then
  CMAKE_ARGS+=(-DSCIP_DIR="${SCIP_DIR}")
fi
if [[ -n "${SCIPPP_SOURCE_DIR:-}" ]]; then
  CMAKE_ARGS+=(-DSCIPPP_SOURCE_DIR="${SCIPPP_SOURCE_DIR}")
fi

echo "==> Setting project."
echo "    Source : ${ROOT_DIR}"
echo "    Build  : ${BUILD_DIR}"
echo "    Tipo   : ${BUILD_TYPE}"
cmake "${CMAKE_ARGS[@]}"

echo "==> Compiling with ${JOBS} job(s)"
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"

if [[ -x "${BUILD_DIR}/sos" ]]; then
  echo "==> Compilation finished. ${BUILD_DIR}/sos"
else
  echo "Warning. Compilation finished but file not found. ${BUILD_DIR}/sos." >&2
fi
