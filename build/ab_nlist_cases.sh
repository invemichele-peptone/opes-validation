#!/bin/bash
# Concurrent A/B timing of two PLUMED builds with build/time_nlist_cases.sh. Both builds run at the
# same time, each on one half of the cores, then the halves are swapped, so that they share the
# machine load. Timings of different builds also carry code-placement noise (up to ~10% here):
# build both with CXXFLAGS="-O3 -falign-functions=64" to remove it.
#
# usage: build/ab_nlist_cases.sh <outdir> <label-a> <label-b> [case ...]
# Writes <outdir>/<label>-<round>/times.dat, round 1 and 2. To pool the rounds:
#   uv run python build/summarize_nlist_cases.py --pool <outdir>/<label>-1 <outdir>/<label>-2
# env: HALF1, HALF2  the two halves of the cores (default: 5-9 and 15-19, the Cortex-X925 cores here)
#      VARIANTS (default: none exact), REPS, SCALE  as in build/time_nlist_cases.sh

set -euo pipefail
out=$1
label_a=$2
label_b=$3
shift 3
root=$(cd "$(dirname "$0")/.." && pwd)
half1=${HALF1:-5 6 7 8 9}
half2=${HALF2:-15 16 17 18 19}
export VARIANTS=${VARIANTS:-none exact}

mkdir -p "$out"
for round in 1 2; do
  if [ $round = 1 ]; then cores_a=$half1; cores_b=$half2; else cores_a=$half2; cores_b=$half1; fi
  CORES="$cores_a" "$root/build/time_nlist_cases.sh" "$out/$label_a-$round" "$label_a" "$@" > "$out/$label_a-$round.log" 2>&1 &
  CORES="$cores_b" "$root/build/time_nlist_cases.sh" "$out/$label_b-$round" "$label_b" "$@" > "$out/$label_b-$round.log" 2>&1 &
  wait
done
grep -h failed "$out"/*.log || true
echo "written $out"
