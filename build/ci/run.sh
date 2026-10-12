#!/bin/bash
# Run a CI-like job (build/ci/job.sh) on a branch or commit of PLUMED, in the plumed-ci:22.04 image.
# The source is taken with git archive, i.e. exactly the committed tree, as the CI sees it.
#
# usage: build/ci/run.sh <branch-or-commit> <job> [<plumed2 repo>]
# The tool versions follow the branch: cppcheck 2.10.3 and doxygen 1.8.17 for 2.9 and 2.10,
# cppcheck 2.13.4 and doxygen 1.13.2 for later versions. Logs go to opt/logs/ci-<label>-<job>.log.
# env: REGTESTS, JOBS (passed to job.sh), CPUS (docker --cpus)
set -euo pipefail
ref=$1
job=$2
repo=${3:-$(cd "$(dirname "$0")/../../../plumed2" && pwd)}
root=$(cd "$(dirname "$0")/../.." && pwd)
label=$(echo "$ref" | tr '/' '_')
work=$(mktemp -d "${TMPDIR:-/tmp}/plumed-ci.XXXXXX")
git -C "$repo" archive "$ref" -o "$work/plumed2.tgz" --prefix plumed2/
version=$(tar xzOf "$work/plumed2.tgz" plumed2/VERSION.txt | grep -v "^#" | head -1)
case $version in
  2.9.* | 2.10.*) cppcheck=2.10.3 doxygen=1_8_17 ;;
  *) cppcheck=2.13.4 doxygen=1_13_2 ;;
esac
log=$root/opt/logs/ci-$label-$job.log
echo "### $ref ($(git -C "$repo" rev-parse --short "$ref"), PLUMED $version): job $job, cppcheck $cppcheck, doxygen $doxygen -> $log"
docker run --rm --user "$(id -u):$(id -g)" ${CPUS:+--cpus "$CPUS"} \
  -e REGTESTS="${REGTESTS:-}" -e JOBS="${JOBS:-16}" \
  -v "$work:/src:ro" -v "$root/build/ci/job.sh:/job.sh:ro" \
  plumed-ci:22.04 bash /job.sh "$job" "$cppcheck" "$doxygen" > "$log" 2>&1 && status=0 || status=$?
rm -rf "$work"
echo "### $ref job $job: $([ $status = 0 ] && echo OK || echo "FAILED (exit $status)")"
exit $status
