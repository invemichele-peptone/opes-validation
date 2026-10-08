#!/bin/bash
# Build and install PLUMED (opes+ves modules, no MPI) from a given git ref.
#
# usage: build/build_plumed.sh <label> <ref> [repo]
#   label  name of the install, e.g. "v2.9-base" -> opt/plumed/<label>
#   ref    branch, tag, commit SHA, or "pr/<N>" for a GitHub pull request head
#   repo   git URL or local path (default: https://github.com/plumed/plumed2.git)
#
# env: JOBS (default: nproc), OPT_DIR (default: <repo-root>/opt)

set -euo pipefail

if [ $# -lt 2 ]; then
  sed -n '2,10p' "$0"
  exit 1
fi
label=$1
ref=$2
repo=${3:-https://github.com/plumed/plumed2.git}

root=$(cd "$(dirname "$0")/.." && pwd)
opt=${OPT_DIR:-$root/opt}
jobs=${JOBS:-$(nproc)}
src=$opt/src/plumed-$label
prefix=$opt/plumed/$label

case $ref in
  pr/*) fetch_ref="pull/${ref#pr/}/head" ;;
  *) fetch_ref=$ref ;;
esac
depth="--depth 1"
[ -d "$repo" ] && depth=""  # local repos may refuse shallow fetches of arbitrary SHAs

echo "### building PLUMED '$label' from $repo @ $ref"
rm -rf "$src"
mkdir -p "$src"
cd "$src"
git init -q
git remote add origin "$repo"
git fetch -q $depth origin "$fetch_ref"
git checkout -q FETCH_HEAD
sha=$(git rev-parse HEAD)
echo "### commit $sha"

configure_opts=(--prefix="$prefix" --enable-modules=+opes+ves --disable-mpi --disable-python --disable-doc)
./configure "${configure_opts[@]}" CXXFLAGS="-O3" > configure.log 2>&1 || { tail -30 configure.log; exit 1; }
make -j "$jobs" > make.log 2>&1 || { tail -30 make.log; exit 1; }
rm -rf "$prefix"
make install > install.log 2>&1 || { tail -30 install.log; exit 1; }

cat > "$prefix/BUILD_INFO" << EOF
label=$label
repo=$repo
ref=$ref
sha=$sha
configure=${configure_opts[*]}
date=$(date -Iseconds)
host=$(uname -nm)
compiler=$(${CXX:-c++} --version | head -1)
EOF
echo "### installed in $prefix"
# the installed binaries have no rpath: libplumedKernel.so is found via LD_LIBRARY_PATH
LD_LIBRARY_PATH=$prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} "$prefix/bin/plumed" info --long-version
