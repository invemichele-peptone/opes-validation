#!/bin/bash
# Generate the GROMACS run inputs of a system: one .tpr per replica.
#
# usage: build/make_tprs.sh <system> [gromacs-version] [replicas]
#   system           directory in systems/ with a gromacs/ subdirectory containing
#                    md.mdp (with @SEED@ placeholder), topol.top and conf.<i>.gro
#   gromacs-version  install in opt/gromacs/<version> (default: 2024.3)
#   replicas         number of replicas (default: number of conf.<i>.gro)
#
# Replica i starts from conf.<i % nconf>.gro with gen_seed = ld_seed = SEED_BASE + i,
# so that the inputs are deterministic and identical for every PLUMED version.

set -euo pipefail

if [ $# -lt 1 ]; then
  sed -n '2,11p' "$0"
  exit 1
fi
system=$1
version=${2:-2024.3}
SEED_BASE=${SEED_BASE:-1000}

root=$(cd "$(dirname "$0")/.." && pwd)
gmx=${OPT_DIR:-$root/opt}/gromacs/$version/bin/gmx
sysdir=$root/systems/$system
[ -x "$gmx" ] || { echo "GROMACS not found: $gmx"; exit 1; }
[ -d "$sysdir/gromacs" ] || { echo "no gromacs inputs in $sysdir"; exit 1; }

work=$(mktemp -d)
trap 'rm -rf -- "$work"' EXIT
mkdir -p "$sysdir/tpr"
nconf=$(ls "$sysdir"/gromacs/conf.*.gro | wc -l)
nrep=${3:-$nconf}
for ((i=0; i<nrep; i++)); do
  sed "s/@SEED@/$((SEED_BASE+i))/" "$sysdir/gromacs/md.mdp" > "$work/md.$i.mdp"
  GMX_MAXBACKUP=-1 "$gmx" grompp -f "$work/md.$i.mdp" -p "$sysdir/gromacs/topol.top" -c "$sysdir/gromacs/conf.$((i%nconf)).gro" \
    -o "$sysdir/tpr/input.$i.tpr" -po "$work/mdout.mdp" > "$work/grompp.$i.log" 2>&1 \
    || { tail -30 "$work/grompp.$i.log"; exit 1; }
  # the same notes are repeated for every replica
  [ $i -eq 0 ] && { grep -h -A4 -E '^(NOTE|WARNING)' "$work/grompp.$i.log" || true; }
done
"$gmx" --version 2>/dev/null | grep 'GROMACS version' > "$sysdir/tpr/GROMACS_VERSION"
echo "### written $nrep tpr files in $sysdir/tpr"
