# CI-like environment for PLUMED (as in .github/workflows of v2.9 and master): ubuntu 22.04, the
# packages installed by the CI, and the pinned tools of both branches (autoconf 2.69, cppcheck 2.10.3
# for v2.9 and 2.13.4 for master, doxygen 1.8.17 for v2.9 and 1.13.2 for master).
# build: docker build -t plumed-ci:22.04 build/ci
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -yq update && apt-get -yq install --no-install-recommends \
      build-essential gfortran git wget ca-certificates cmake flex bison m4 gawk xxd vim-common \
      python3 python3-pip python3-venv python3-dev libpcre3-dev zlib1g-dev \
      libatlas-base-dev libfftw3-dev gsl-bin libgsl-dev libboost-serialization-dev \
      libopenmpi-dev openmpi-bin graphviz doxygen-latex ghostscript \
 && apt-get -yq remove doxygen \
 && rm -rf /var/lib/apt/lists/*
ENV CC=gcc CXX=g++
# autoconf 2.69, as .ci/install.autoconf
RUN cd /tmp && wget -q https://ftpmirror.gnu.org/autoconf/autoconf-2.69.tar.gz \
 && tar xzf autoconf-2.69.tar.gz && cd autoconf-2.69 && ./configure --prefix=/opt/autoconf && make && make install \
 && rm -rf /tmp/autoconf-2.69*
# cppcheck, as .ci/install.cppcheck (MATCHCOMPILER=yes HAVE_RULES=yes), one prefix per version
RUN for v in 2.10.3 2.13.4; do \
      cd /tmp && git clone -q https://github.com/danmar/cppcheck.git --branch $v --depth 1 cppcheck-$v && cd cppcheck-$v \
      && make -j 16 install MATCHCOMPILER=yes FILESDIR=/opt/cppcheck-$v/share/cppcheck HAVE_RULES=yes \
         CXXFLAGS="-DNDEBUG -O2 -Wall -Wno-sign-compare -Wno-unused-function" PREFIX=/opt/cppcheck-$v > /dev/null \
      && cd /tmp && rm -rf cppcheck-$v && /opt/cppcheck-$v/bin/cppcheck --version; done
# doxygen, as .ci/install.doxygen, one prefix per version
RUN for v in 1_8_17 1_13_2; do \
      cd /tmp && git clone -q https://github.com/doxygen/doxygen.git --branch Release_$v --depth 1 doxygen-$v \
      && mkdir doxygen-$v/build && cd doxygen-$v/build \
      && cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=/opt/doxygen-$v .. > /dev/null && make -j 16 > /dev/null && make install > /dev/null \
      && cd /tmp && rm -rf doxygen-$v && /opt/doxygen-$v/bin/doxygen --version; done
# python: the manual of master (new-manual/requirements.txt) and the codecheck of master (ruff, Cython)
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install -q --upgrade pip \
 && /opt/venv/bin/pip install -q numpy requests PlumedToHTML networkx mkdocs-material python-markdown-math Cython==3.2.4 ruff==0.14.14 pytest
ENV OMPI_MCA_btl_base_warn_component_unused=0 OMPI_MCA_btl_base_verbose=0 OMPI_MCA_plm=isolated \
    OMPI_MCA_btl_vader_single_copy_mechanism=none OMPI_MCA_rmaps_base_oversubscribe=yes OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 PLUMED_STACK_TRACE=yes PLUMED_NUM_THREADS=2
WORKDIR /work
# curl: used by the manual of master to fetch citations; git identity for the codecheck job
RUN apt-get -yq update && apt-get -yq install --no-install-recommends curl && rm -rf /var/lib/apt/lists/* \
 && git config --system user.name ci && git config --system user.email ci@localhost \
 && git config --system --add safe.directory '*'
