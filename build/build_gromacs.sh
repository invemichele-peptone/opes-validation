#!/bin/bash
# Build and install GROMACS patched with PLUMED in *runtime* mode.
#
# usage: build/build_gromacs.sh <version> <plumed-label>
#   version       GROMACS version, must have a patch in the PLUMED install (e.g. 2024.3)
#   plumed-label  PLUMED install (opt/plumed/<label>) used to apply the patch
#
# With the runtime patch the PLUMED kernel is loaded at run time from the
# PLUMED_KERNEL environment variable, so a single GROMACS build can be used with
# any PLUMED version. Use the *oldest* PLUMED you want to test to patch GROMACS:
# newer kernels stay backward compatible with older patches, not vice versa.
#
# env: JOBS (default: nproc), OPT_DIR (default: <repo-root>/opt),
#      MPI: MPI install prefix (e.g. opt/openmpi); builds gmx_mpi into opt/gromacs/<version>-mpi.
#           Patch it with an MPI build of PLUMED, and run it with MPI builds of PLUMED.

set -euo pipefail

if [ $# -ne 2 ]; then
  sed -n '2,13p' "$0"
  exit 1
fi
version=$1
plumed_label=$2

root=$(cd "$(dirname "$0")/.." && pwd)
opt=${OPT_DIR:-$root/opt}
jobs=${JOBS:-$(nproc)}
plumed_prefix=$opt/plumed/$plumed_label
src=$opt/src/gromacs-$version
prefix=$opt/gromacs/$version
[ -n "${MPI:-}" ] && MPI=$(cd "$MPI" && pwd) && src=$src-mpi && prefix=$prefix-mpi
tarball=$opt/src/gromacs-$version.tar.gz

[ -x "$plumed_prefix/bin/plumed" ] || { echo "PLUMED '$plumed_label' not found in $plumed_prefix"; exit 1; }

echo "### building GROMACS $version, runtime-patched with PLUMED '$plumed_label'"
mkdir -p "$opt/src"
[ -f "$tarball" ] || curl -fsSL -o "$tarball" "https://ftp.gromacs.org/gromacs/gromacs-$version.tar.gz"
rm -rf "$src" "$opt/src/gromacs-$version"
tar -xzf "$tarball" -C "$opt/src"
[ "$src" != "$opt/src/gromacs-$version" ] && mv "$opt/src/gromacs-$version" "$src"
cd "$src"
PATH=$plumed_prefix/bin:$PATH LD_LIBRARY_PATH=$plumed_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} plumed patch -p -e "gromacs-$version" --runtime > patch.log 2>&1 || { cat patch.log; exit 1; }

# FFTs are only needed for PME, which none of these systems use
cmake_opts=(
  -DCMAKE_INSTALL_PREFIX="$prefix"
  -DCMAKE_BUILD_TYPE=Release
  -DGMX_MPI=$([ -n "${MPI:-}" ] && echo ON || echo OFF)
  -DGMX_GPU=OFF
  -DGMX_FFT_LIBRARY=fftpack
  -DGMX_HWLOC=OFF
  -DGMX_INSTALL_LEGACY_API=OFF
  -DREGRESSIONTEST_DOWNLOAD=OFF
)
if [ -n "${MPI:-}" ]; then
  export PATH=$MPI/bin:$PATH CC=$MPI/bin/mpicc CXX=$MPI/bin/mpicxx
  cmake_opts+=(-DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX")
fi
mkdir build
cd build
cmake .. "${cmake_opts[@]}" > cmake.log 2>&1 || { tail -30 cmake.log; exit 1; }
make -j "$jobs" > make.log 2>&1 || { tail -30 make.log; exit 1; }
rm -rf "$prefix"
make install > install.log 2>&1 || { tail -30 install.log; exit 1; }

cat > "$prefix/BUILD_INFO" << EOF
version=$version
plumed_patch_from=$plumed_label ($(grep '^sha=' "$plumed_prefix/BUILD_INFO" | cut -d= -f2))
plumed_patch_mode=runtime
cmake=${cmake_opts[*]}
date=$(date -Iseconds)
host=$(uname -nm)
compiler=$(${CXX:-c++} --version | head -1)
EOF
echo "### installed in $prefix"
"$prefix/bin/"gmx* --version 2>/dev/null | grep -E 'GROMACS version|SIMD instructions|Precision'
