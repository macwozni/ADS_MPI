#!/usr/bin/env bash
set -eu

mpiexec=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
probe=${1:-./legacy_abort_probe}
output=$(mktemp)
trap 'rm -f "$output"' EXIT

set +e
timeout 20s "$mpiexec" -n 2 "$probe" >"$output" 2>&1
status=$?
set -e

if [ "$status" -eq 0 ]; then
    printf 'FAIL: legacy MUMPS failure returned success\n' >&2
    sed -n '1,80p' "$output" >&2
    exit 1
fi
if [ "$status" -eq 124 ]; then
    printf 'FAIL: legacy MUMPS failure left another MPI rank hanging\n' >&2
    sed -n '1,80p' "$output" >&2
    exit 1
fi
if ! grep -F 'MUMPS failed during analysis' "$output" >/dev/null; then
    printf 'FAIL: expected MUMPS analysis diagnostic is missing\n' >&2
    sed -n '1,80p' "$output" >&2
    exit 1
fi

printf 'OK (legacy five-argument MUMPS failure aborts the complete MPI job)\n'
