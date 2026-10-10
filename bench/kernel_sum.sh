#!/bin/bash
# Run the kernel-sum micro-benchmark (bench/kernel_sum.cpp) on the ala4 kernels of a quick-tier run.
# The kernels are compressed by replaying the trajectory with one kernel per frame (as in
# build/replay_nlist.sh); the bias is then evaluated along the same trajectory.
#
# usage: bench/kernel_sum.sh <plumed-label> <run-set-with-trajectories> [outdir]
# env: CORE, the CPU core to pin the benchmark to (default 4)

set -euo pipefail
label=$1
runs=$2
root=$(cd "$(dirname "$0")/.." && pwd)
out=${3:-$root/runs/kernel-sum}
prefix=$root/opt/plumed/$label
mkdir -p "$out"
g++ -O3 -o "$out/kernel_sum" "$root/bench/kernel_sum.cpp"
for method in opes explore; do
  action=$([ $method = explore ] && echo OPES_METAD_EXPLORE || echo OPES_METAD)
  d=$out/ala4-$method
  mkdir -p "$d"
  cat > "$d/plumed.dat" << EOF
psi1: TORSION ATOMS=7,9,15,17
psi2: TORSION ATOMS=17,19,25,27
psi3: TORSION ATOMS=27,29,35,37
opes: $action ARG=psi1,psi2,psi3 PACE=1 BARRIER=50 TEMP=300 CUTOFF_SCHEME=PR FILE=/dev/null STATE_WFILE=STATE STATE_WSTRIDE=20000 FMT=%.17g
PRINT ARG=psi1,psi2,psi3 FILE=COLVAR FMT=%.17g
EOF
  (cd "$d" && rm -f STATE COLVAR && LD_LIBRARY_PATH=$prefix/lib "$prefix/bin/plumed" driver --plumed plumed.dat \
    --mf_xtc "$root/runs/$runs/quick/ala4/$method/000/md.xtc" > log.txt 2>&1)
  echo "== ala4 $action"
  taskset -c "${CORE:-4}" "$out/kernel_sum" "$d/STATE" "$d/COLVAR" 3
done
