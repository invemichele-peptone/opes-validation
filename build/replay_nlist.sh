#!/bin/bash
# Replay saved ala2 and ala4 trajectories with four OPES actions side by side: no neighbor list,
# the v2.9 one, NLIST_REBUILD=EXACT, and NLIST_REBUILD=EXACT NLIST_ZED=FULL (PR cutoff, PLAIN merge).
# Kernels are deposited at every frame, i.e. every 500 MD steps as in the runs.
#
# usage: build/replay_nlist.sh <plumed-label> <run-set-with-trajectories> [outdir]
# Writes <outdir>/<system>-<method>-<replica>/{bias.out,der.out}

set -euo pipefail
label=$1
runs=$2
root=$(cd "$(dirname "$0")/.." && pwd)
out=${3:-$root/runs/nlist-replay}
prefix=$root/opt/plumed/$label

for system in ala2 ala4; do
  if [ $system = ala2 ]; then
    cvs=$'phi: TORSION ATOMS=5,7,9,15\npsi: TORSION ATOMS=7,9,15,17'; args=phi,psi; barrier=30
  else
    cvs=$'psi1: TORSION ATOMS=7,9,15,17\npsi2: TORSION ATOMS=17,19,25,27\npsi3: TORSION ATOMS=27,29,35,37'
    args=psi1,psi2,psi3; barrier=50
  fi
  for method in opes explore; do
    action=$([ $method = explore ] && echo OPES_METAD_EXPLORE || echo OPES_METAD)
    common="ARG=$args PACE=1 BARRIER=$barrier TEMP=300 CUTOFF_SCHEME=PR FILE=/dev/null"
    for rep in 000 001 002; do
      d=$out/$system-$method-$rep
      mkdir -p "$d"
      cat > "$d/plumed.dat" << EOF
$cvs
ref: $action $common
nlp: $action $common NLIST
nla: $action $common NLIST NLIST_REBUILD=EXACT
nlab: $action $common NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL
PRINT ARG=ref.bias,nlp.bias,nla.bias,nlab.bias,ref.zed,nlp.zed,nla.zed,nlab.zed,ref.nker,nlp.nlker,nlab.nlker,nlp.nlsteps,nlab.nlsteps FILE=bias.out FMT=%.14g
DUMPDERIVATIVES ARG=ref.bias,nlp.bias,nla.bias,nlab.bias FILE=der.out FMT=%.14g
EOF
      (cd "$d" && LD_LIBRARY_PATH=$prefix/lib "$prefix/bin/plumed" driver --plumed plumed.dat \
        --mf_xtc "$root/runs/$runs/quick/$system/$method/$rep/md.xtc" > log.txt 2>&1 || echo "failed $d") &
    done
  done
done
wait
echo "written $out"
