#!/bin/bash
# Run one CI-like job on a PLUMED source tarball, inside the plumed-ci:22.04 image (build/ci/Dockerfile).
#
# usage (on the host): build/ci/run.sh <branch-or-commit> <job> [<plumed2 repo>]
# inside the container: job.sh <job> <cppcheck-version> <doxygen-version>
# jobs, following .github/workflows of PLUMED (v2.9: ci.yml; master: linuxWF.yml, codecheckWF.yml):
#   doc-mpi    MPI build with all modules and --enable-pdfdoc, make nmcheck, all regtests, make doc
#   debug      --enable-debug --enable-debug-glibcxx, all modules (without boost serialization), make nmcheck, regtests
#   debug-mpi  the same with MPI
#   codecheck  configure --enable-modules=all, make -C src/lib dirslinks, make codecheck
# env: REGTESTS  regtest subdirectories to run (default: all); JOBS (default: 16)
set -euo pipefail
job=$1
cppcheck_version=$2
doxygen_version=$3
jobs=${JOBS:-16}
export PATH=/opt/autoconf/bin:/opt/cppcheck-$cppcheck_version/bin:/opt/doxygen-$doxygen_version/bin:/opt/venv/bin:$PATH
export HOME=/tmp/home
mkdir -p $HOME /tmp/build
cd /tmp/build
tar xzf /src/plumed2.tgz
cd plumed2
prefix=/tmp/opt
step() { echo; echo "######## $(date +%T) $*"; }
regtests() {
  step "regtests (PLUMED_NUM_THREADS=$PLUMED_NUM_THREADS): ${REGTESTS:-all}"
  if [ -n "${REGTESTS:-}" ]; then
    for d in $REGTESTS; do make --no-print-directory -C regtest/$d > /dev/null 2>&1 || true; done
  else
    make --no-print-directory -j "$jobs" -C regtest testclean > /dev/null 2>&1 || true
  fi
  # as in the CI: checkfail lists the failed tests and fails if there is any
  make -C regtest checkfail
}
case $job in
  doc-mpi)
    step configure
    ./configure CC=mpicc CXX=mpic++ --enable-boost_serialization --disable-dependency-tracking --enable-modules=all \
      --enable-pdfdoc --prefix=$prefix > configure.log 2>&1 || { tail -40 configure.log; exit 1; }
    grep -E 'make_doc|make_pdfdoc' Makefile.conf
    step build; make -j "$jobs" > make.log 2>&1 || { grep -B5 -A5 'error' make.log | head -80; exit 1; }
    make install > install.log 2>&1
    step nmcheck; make nmcheck
    export PATH=$prefix/bin:$PATH LD_LIBRARY_PATH=$prefix/lib:${LD_LIBRARY_PATH:-}
    regtests
    step doc; export PLUMED_MAXBACKUP=0
    make doc > doc.log 2>&1 && echo "make doc: OK" || { echo "make doc: FAILED"; tail -60 doc.log; exit 1; }
    ;;
  debug | debug-mpi)
    mpi=(); [ $job = debug-mpi ] && mpi=(CC=mpicc CXX=mpic++)
    step configure
    ./configure "${mpi[@]}" --disable-dependency-tracking --enable-modules=all --enable-debug --enable-debug-glibcxx \
      --prefix=$prefix > configure.log 2>&1 || { tail -40 configure.log; exit 1; }
    step build; make -j "$jobs" > make.log 2>&1 || { grep -B5 -A5 'error' make.log | head -80; exit 1; }
    make install > install.log 2>&1
    step nmcheck; make nmcheck
    export PATH=$prefix/bin:$PATH LD_LIBRARY_PATH=$prefix/lib:${LD_LIBRARY_PATH:-}
    regtests
    ;;
  codecheck)
    # codecheck lists the files with git, as in the CI checkout
    git init -q && git add -A && git commit -q -m "source of the job"
    step "configure (autoconf $(autoconf --version | head -1), $(cppcheck --version))"
    ./configure --enable-modules=all > configure.log 2>&1 || { tail -40 configure.log; exit 1; }
    make -j "$jobs" -C astyle > /dev/null 2>&1
    make -C src/lib/ dirslinks > /dev/null
    step codecheck; make codecheck
    ;;
  *) echo "unknown job $job"; exit 1 ;;
esac
step "job $job: OK"
