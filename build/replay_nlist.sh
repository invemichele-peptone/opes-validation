#!/bin/bash
# Replay saved ala2 and ala4 trajectories with four OPES actions side by side: no neighbor list,
# the v2.9 one, NLIST_REBUILD=EXACT, and NLIST_REBUILD=EXACT NLIST_ZED=FULL (PR cutoff, PLAIN merge).
# Kernels are deposited at every frame, i.e. every 500 MD steps as in the runs.
#
# usage: build/replay_nlist.sh <plumed-label> <run-set-with-trajectories> [outdir]
# Writes <outdir>/<system>-<method>-<replica>/{bias.out,der.out}, or with WALKERS
# <outdir>/<system>-<method>-mw/{bias,der}.<walker>.out
# env: THREADS  OpenMP threads used by PLUMED (default 1)
#      NP       MPI processes per replica or walker (default 1; needs a PLUMED built with MPI)
#      WALKERS  if set, replicas 0..WALKERS-1 are the walkers of one WALKERS_MPI run
#      MPIRUN   MPI launcher (default <repo-root>/opt/openmpi/bin/mpirun)

set -euo pipefail
label=$1
runs=$2
root=$(cd "$(dirname "$0")/.." && pwd)
out=${3:-$root/runs/nlist-replay}
prefix=$root/opt/plumed/$label
threads=${THREADS:-1}
np=${NP:-1}
walkers=${WALKERS:-}
mpirun=${MPIRUN:-$root/opt/openmpi/bin/mpirun}
export LD_LIBRARY_PATH=$prefix/lib:$root/opt/openmpi/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

nprocs=$((np * ${walkers:-1}))
launch=()
if [ "$nprocs" -gt 1 ]; then
  launch=("$mpirun" -np "$nprocs" --oversubscribe --bind-to none -x PLUMED_NUM_THREADS -x LD_LIBRARY_PATH)
fi
multi=()
[ -n "$walkers" ] && multi=(--multi "$walkers")
maxjobs=$(( $(nproc) / (nprocs * threads) ))
[ "$maxjobs" -lt 1 ] && maxjobs=1

replay() {  # dir trajectory
  (cd "$1" && PLUMED_NUM_THREADS=$threads "${launch[@]}" "$prefix/bin/plumed" driver "${multi[@]}" \
    --plumed plumed.dat --mf_xtc "$2" > log.txt 2>&1 || echo "failed $1") &
  while [ "$(jobs -rp | wc -l)" -ge "$maxjobs" ]; do wait -n || true; done
}

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
    [ -n "$walkers" ] && common+=" WALKERS_MPI"
    reps=(000 001 002)
    [ -n "$walkers" ] && reps=(mw)
    for rep in "${reps[@]}"; do
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
      if [ -n "$walkers" ]; then
        # driver --multi reads md.<walker>.xtc
        for ((w = 0; w < walkers; w++)); do
          ln -sf "$root/runs/$runs/quick/$system/$method/$(printf %03d $w)/md.xtc" "$d/md.$w.xtc"
        done
        replay "$d" md.xtc
      else
        replay "$d" "$root/runs/$runs/quick/$system/$method/$rep/md.xtc"
      fi
    done
  done
done
wait
echo "written $out"
