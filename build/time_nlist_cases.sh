#!/bin/bash
# Cost of OPES with and without the exact neighbor list, for low-dimensional cases (1 and 2 CVs),
# measured with PLUMED's detailed timers. Jobs are shuffled and each runs pinned to one core.
#
# usage: build/time_nlist_cases.sh <outdir> <plumed-label> [case ...]
# Cases (default: all short ones):
#   mueller1d-opes mueller1d-explore mueller2d-opes mueller2d-explore wolfe1d-opes   (ves engine)
#   ala2-1d-opes ala2-1d-explore ala2-2d-opes ala2-2d-explore ala2-2d-jpcl          (GROMACS)
#   ala4-3d-opes ala4-3d-explore                                                   (GROMACS, not by default)
# Writes <outdir>/times.dat:
#   case variant replica nsteps wall_time opes_time nker mean_nlker mean_nlsteps update_time
# opes_time is the time of the OPES action: forward loop, backward loop and update (the update
# also includes the other actions, e.g. PRINT, which is negligible). update_time is the update alone.
# env: VARIANTS  among none, exact, exact_p2, exact_p5, exact_p8 (default: none exact exact_p2)
#      CORES     cores to use, one job per core at a time (default: all)
#      REPS      replicas per case (default 3)
#      SCALE     multiplies the number of steps of every case (default 1)
#      SUMMARY_ONLY  if set, only rewrite times.dat from runs already done

set -euo pipefail
out=$1
label=$2
shift 2
root=$(cd "$(dirname "$0")/.." && pwd)
prefix=$root/opt/plumed/$label
gmx=$root/opt/gromacs/2024.3/bin/gmx
reps=${REPS:-3}
scale=${SCALE:-1}
read -r -a variants <<< "${VARIANTS:-none exact exact_p2}"
read -r -a cores <<< "${CORES:-$(seq -s ' ' 0 $(($(nproc) - 1)))}"
cases=("$@")
[ ${#cases[@]} -eq 0 ] && cases=(mueller1d-opes mueller1d-explore mueller2d-opes mueller2d-explore wolfe1d-opes
  ala2-1d-opes ala2-1d-explore ala2-2d-opes ala2-2d-explore ala2-2d-jpcl)

declare -A extra=(
  [none]=""
  [exact]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL"
  [exact_p2]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=2,0.4"
  [exact_p5]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=5,0.5"
  [exact_p8]="NLIST NLIST_REBUILD=EXACT NLIST_ZED=FULL NLIST_PARAMETERS=8,0.5"
)
# case -> system, input, sed expression on the ARG line, steps
declare -A system=( [mueller1d-opes]=mueller [mueller1d-explore]=mueller [mueller2d-opes]=mueller
  [mueller2d-explore]=mueller [wolfe1d-opes]=wolfe-quapp [ala2-1d-opes]=ala2 [ala2-1d-explore]=ala2
  [ala2-2d-opes]=ala2 [ala2-2d-explore]=ala2 [ala2-2d-jpcl]=ala2 [ala4-3d-opes]=ala4 [ala4-3d-explore]=ala4 )
declare -A input=( [mueller1d-opes]=opes [mueller1d-explore]=explore [mueller2d-opes]=opes
  [mueller2d-explore]=explore [wolfe1d-opes]=opes [ala2-1d-opes]=opes [ala2-1d-explore]=explore
  [ala2-2d-opes]=opes [ala2-2d-explore]=explore [ala2-2d-jpcl]=opes_jpcl [ala4-3d-opes]=opes [ala4-3d-explore]=explore )
declare -A argsed=( [mueller2d-opes]='s/ARG=p.x$/ARG=p.x,p.y/' [mueller2d-explore]='s/ARG=p.x$/ARG=p.x,p.y/'
  [ala2-1d-opes]='s/ARG=phi,psi$/ARG=phi/' [ala2-1d-explore]='s/ARG=phi,psi$/ARG=phi/' )
declare -A steps=( [mueller1d-opes]=20000000 [mueller1d-explore]=20000000 [mueller2d-opes]=20000000
  [mueller2d-explore]=20000000 [wolfe1d-opes]=4000000 [ala2-1d-opes]=5000000 [ala2-1d-explore]=5000000
  [ala2-2d-opes]=5000000 [ala2-2d-explore]=5000000 [ala2-2d-jpcl]=5000000 [ala4-3d-opes]=5000000
  [ala4-3d-explore]=5000000 )

mkdir -p "$out"
jobsfile=$out/jobs.txt
: > "$jobsfile"
for c in "${cases[@]}"; do
  s=${system[$c]}
  nsteps=$((steps[$c] * scale))
  for variant in "${variants[@]}"; do
    for ((rep = 0; rep < reps; rep++)); do
      d=$out/$c-$variant-$rep
      mkdir -p "$d"
      # drop any NLIST line, add the timers, the ARG change of the case and the keywords of the variant
      grep -v -E '^\s*NLIST\s*$' "$root/systems/$s/plumed-${input[$c]}.dat" \
        | sed -E "${argsed[$c]:-}" \
        | sed -E "s/^(opes: OPES_METAD(_EXPLORE)? \.\.\.)$/\1\n  CUTOFF_SCHEME=PR ${extra[$variant]}/" \
        | sed -E 's/^(opes: OPES_METAD(_EXPLORE)? \.\.\.)$/DEBUG DETAILED_TIMERS\n\1/' > "$d/plumed.dat"
      if [ $s = ala2 ] || [ $s = ala4 ]; then
        echo "$c $variant $rep $nsteps $d $gmx mdrun -s $root/systems/$s/tpr/input.$rep.tpr -plumed plumed.dat -deffnm md -nsteps $nsteps -ntmpi 1 -ntomp 1 -pin off -nb cpu" >> "$jobsfile"
      else
        read -r _ seed x0 y0 <<< "$(grep -v '^#' "$root/systems/$s/replicas.dat" | sed -n "$((rep + 1))p")"
        sed "s/@NSTEP@/$nsteps/; s/@SEED@/$seed/; s/@X0@/$x0/; s/@Y0@/$y0/" "$root/systems/$s/md_input.dat" > "$d/md_input.dat"
        cp "$root/systems/$s/potential.dat" "$d/"
        echo "$c $variant $rep $nsteps $d $prefix/bin/plumed ves_md_linearexpansion md_input.dat" >> "$jobsfile"
      fi
    done
  done
done

[ -n "${SUMMARY_ONLY:-}" ] || {
shuf "$jobsfile" > "$out/order.txt"
for ((k = 0; k < ${#cores[@]}; k++)); do
  awk -v n=${#cores[@]} -v k=$k 'NR % n == k' "$out/order.txt" | while read -r c variant rep nsteps d cmd; do
    (cd "$d" && /usr/bin/time -f %e -o time.txt env PLUMED_KERNEL="$prefix/lib/libplumedKernel.so" \
      LD_LIBRARY_PATH="$prefix/lib" OMP_NUM_THREADS=1 taskset -c "${cores[$k]}" $cmd > out.log 2>&1 || echo "failed $d")
  done &
done
wait
}

echo "# case variant replica nsteps wall_time opes_time nker mean_nlker mean_nlsteps update_time" > "$out/times.dat"
while read -r c variant rep nsteps d cmd; do
  log=$d/md.log
  [ -f "$log" ] || log=$d/plumed.log
  opes_time=$(awk '/^PLUMED: 4A +[0-9]+ opes / || /^PLUMED: 5A +[0-9]+ opes / {t+=$6} /^PLUMED: 6 Update / {t+=$5} END {print t}' "$log")
  update_time=$(awk '/^PLUMED: 6 Update / {print $5}' "$log")
  stats=$(awk '/^#! FIELDS/ {for (i=3; i<=NF; i++) col[$i]=i-2}
               !/^#/ {n++; nk=$col["opes.nker"]; if ("opes.nlker" in col) {nl+=$col["opes.nlker"]; if ($col["opes.nlsteps"]>0) {st+=$col["opes.nlsteps"]; ns++}}}
               END {printf "%d %.1f %.1f", nk, (("opes.nlker" in col) ? nl/n : nk), ((ns>0) ? st/ns : 0)}' "$d/COLVAR")
  echo "$c $variant $rep $nsteps $(cat "$d/time.txt") $opes_time $stats $update_time" >> "$out/times.dat"
done < "$jobsfile"
echo "written $out/times.dat"
