#!/bin/bash
# Time ala2 and ala4 MD runs with different neighbor list settings (PR cutoff, PLAIN merge).
#
# usage: build/time_nlist.sh <plumed-label> <outdir> [jobs]
# Writes <outdir>/times.dat: system method variant replica wall_time nker mean_nlker mean_nlsteps

set -euo pipefail
label=$1
out=$2
jobs=${3:-18}
root=$(cd "$(dirname "$0")/.." && pwd)
gmx=$root/opt/gromacs/2024.3/bin/gmx
kernel=$root/opt/plumed/$label/lib/libplumedKernel.so
nsteps=5000000  # 10 ns

declare -A extra=(
  [none]=""
  [plain]="NLIST"
  [exact]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL"
  [exact_p1.5]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=1.5,0.3"
  [exact_p2]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=2,0.4"
  [exact_p5]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=5,0.5"
)
mkdir -p "$out"
jobsfile=$out/jobs.txt
: > "$jobsfile"
for system in ala2 ala4; do
  for method in opes explore; do
    for variant in "${!extra[@]}"; do
      for rep in 0 1 2; do
        d=$out/$system-$method-$variant-$rep
        mkdir -p "$d"
        # the inputs use NLIST: drop it, then add the keywords of the variant to the OPES action
        grep -v -E '^\s*NLIST\s*$' "$root/systems/$system/plumed-$method.dat" \
          | sed -E "s/^(opes: OPES_METAD(_EXPLORE)? \.\.\.)$/\1\n  CUTOFF_SCHEME=PR ${extra[$variant]}/" > "$d/plumed.dat"
        echo "$system $method $variant $rep $d" >> "$jobsfile"
      done
    done
  done
done
shuf "$jobsfile" | while read -r system method variant rep d; do
  echo "cd $d && /usr/bin/time -f %e -o time.txt env PLUMED_KERNEL=$kernel OMP_NUM_THREADS=1 $gmx mdrun -s $root/systems/$system/tpr/input.$rep.tpr -plumed plumed.dat -deffnm md -nsteps $nsteps -ntmpi 1 -ntomp 1 -pin off -nb cpu > out.log 2>&1"
done | xargs -P "$jobs" -I{} bash -c "{}"

echo "# system method variant replica wall_time nker mean_nlker mean_nlsteps" > "$out/times.dat"
while read -r system method variant rep d; do
  stats=$(awk '/^#! FIELDS/ {for (i=3; i<=NF; i++) col[$i]=i-2}
               !/^#/ {n++; nk=$col["opes.nker"]; if ("opes.nlker" in col) {nl+=$col["opes.nlker"]; st+=$col["opes.nlsteps"]}}
               END {printf "%d %.1f %.1f", nk, ("opes.nlker" in col) ? nl/n : nk, ("opes.nlker" in col) ? st/n : 0}' "$d/COLVAR")
  echo "$system $method $variant $rep $(cat "$d/time.txt") $stats" >> "$out/times.dat"
done < "$jobsfile"
echo "written $out/times.dat"
