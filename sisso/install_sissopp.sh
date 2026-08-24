#!/usr/bin/env bash
set -euo pipefail

commit="43b99110118a51b9f4983b02c8d781ae6f25456c"
repository="https://gitlab.com/sissopp_developers/sissopp.git"
source_dir="${1:-${PWD}/sissopp-src}"
build_dir="${2:-${source_dir}/build-public}"
jobs="${JOBS:-4}"

command -v git >/dev/null || { echo "git is required" >&2; exit 1; }
command -v cmake >/dev/null || { echo "CMake 3.20.2 or newer is required" >&2; exit 1; }

if [[ ! -d "${source_dir}/.git" ]]; then
    git clone --recurse-submodules "${repository}" "${source_dir}"
fi

git -C "${source_dir}" fetch origin "${commit}"
git -C "${source_dir}" checkout --detach "${commit}"
git -C "${source_dir}" submodule update --init --recursive

cmake -S "${source_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DSISSO_BUILD_PYTHON=OFF \
    -DSISSO_BUILD_PARAMS=OFF \
    -DSISSO_BUILD_TESTS=OFF \
    -DSISSO_BUILD_EXE=ON
cmake --build "${build_dir}" --parallel "${jobs}"

binary="${build_dir}/bin/sisso++"
if [[ ! -x "${binary}" ]]; then
    echo "Build finished without the expected executable: ${binary}" >&2
    exit 1
fi
printf 'SISSO++ executable: %s\n' "${binary}"

