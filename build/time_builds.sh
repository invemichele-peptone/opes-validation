#!/bin/bash
# Time ala2 and ala4 MD runs with different PLUMED builds and neighbor list settings (PR cutoff,
# PLAIN merge). Jobs of all builds are shuffled together, and each runs pinned to one core.
#
# usage: build/time_builds.sh <outdir> <plumed-label> [<plumed-label> ...]
# Writes <outdir>/times.dat: label system method variant replica wall_time nker
# env: VARIANTS  neighbor list settings, among none, exact (default: "none exact")
#      CORES     cores to use, one job per core at a time (default: all)
#      NSTEPS    MD steps per run (default 5000000, 10 ns)
#      REPS      replicas per case (default 3)

set -euo pipefail
out=$1
shift
labels=("$@")
root=$(cd "$(dirname "$0")/.." && pwd)
gmx=$root/opt/gromacs/2024.3/bin/gmx
nsteps=${NSTEPS:-5000000}
reps=${REPS:-3}
read -r -a variants <<< "${VARIANTS:-none exact}"
read -r -a cores <<< "${CORES:-$(seq -s ' ' 0 $(($(nproc) - 1)))}"

declare -A extra=(
  [none]=""
  [exact]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL"
)
mkdir -p "$out"
jobsfile=$out/jobs.txt
: > "$jobsfile"
for label in "${labels[@]}"; do
  for system in ala2 ala4; do
    for method in opes explore; do
      for variant in "${variants[@]}"; do
        for ((rep = 0; rep < reps; rep++)); do
          d=$out/$label-$system-$method-$variant-$rep
          mkdir -p "$d"
          # the inputs use NLIST: drop it, then add the keywords of the variant to the OPES action
          grep -v -E '^\s*NLIST\s*$' "$root/systems/$system/plumed-$method.dat" \
            | sed -E "s/^(opes: OPES_METAD(_EXPLORE)? \.\.\.)$/\1\n  CUTOFF_SCHEME=PR ${extra[$variant]}/" > "$d/plumed.dat"
          echo "$label $system $method $variant $rep $d" >> "$jobsfile"
        done
      done
    done
  done
done

# deal the shuffled jobs to the cores, each core runs its share one after the other
shuf "$jobsfile" > "$out/order.txt"
for ((c = 0; c < ${#cores[@]}; c++)); do
  awk -v n=${#cores[@]} -v c=$c 'NR % n == c' "$out/order.txt" | while read -r label system method variant rep d; do
    kernel=$root/opt/plumed/$label/lib/libplumedKernel.so
    (cd "$d" && /usr/bin/time -f %e -o time.txt env PLUMED_KERNEL="$kernel" OMP_NUM_THREADS=1 taskset -c "${cores[$c]}" \
      "$gmx" mdrun -s "$root/systems/$system/tpr/input.$rep.tpr" -plumed plumed.dat -deffnm md -nsteps "$nsteps" \
      -ntmpi 1 -ntomp 1 -pin off -nb cpu > out.log 2>&1 || echo "failed $d")
  done &
done
wait

echo "# label system method variant replica wall_time nker" > "$out/times.dat"
while read -r label system method variant rep d; do
  nker=$(awk '/^#! FIELDS/ {for (i=3; i<=NF; i++) if ($i=="opes.nker") c=i-2} !/^#/ {n=$c} END {print n}' "$d/COLVAR")
  echo "$label $system $method $variant $rep $(cat "$d/time.txt") $nker" >> "$out/times.dat"
done < "$jobsfile"
echo "written $out/times.dat"
