#!/usr/bin/env bash

set -u

MPIEXEC=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
PROBE=${PROBE:-./ads_lifecycle_probe}
output_file=$(mktemp)

cleanup() {
    rm -f -- "$output_file"
}
trap cleanup EXIT

"$MPIEXEC" -n 2 "$PROBE" >"$output_file" 2>&1
status=$?

if [[ $status -ne 0 ]] && \
   grep -Fq -- 'Number of degrees of freedom smaller than number of processors' \
       "$output_file"; then
    printf 'PASS initialize_setup rejects fewer DOFs than MPI ranks\n'
    printf 'OK (1 ADS lifecycle CLI check)\n'
    exit 0
fi

printf 'FAIL initialize_setup rejects fewer DOFs than MPI ranks (status %d)\n' \
    "$status"
sed 's/^/  /' "$output_file"
exit 1
