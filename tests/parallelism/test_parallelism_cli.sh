#!/usr/bin/env bash

set -u

MPIEXEC=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
PROBE=${PROBE:-./parallelism_probe}

checks=0
failures=0

expect_failure() {
    local label=$1
    local expected_message=$2
    shift 2

    local output_file
    local status
    output_file=$(mktemp)

    "$MPIEXEC" -n 1 "$PROBE" "$@" >"$output_file" 2>&1
    status=$?
    checks=$((checks + 1))

    if [[ $status -ne 0 ]] && grep -Fq "$expected_message" "$output_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (status %d)\n' "$label" "$status"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    fi

    rm -f "$output_file"
}

expect_failure 'zero X process-grid dimension' \
    'Process-grid dimensions must be positive:' 0 1 1
expect_failure 'negative Y process-grid dimension' \
    'Process-grid dimensions must be positive:' 1 -1 1
expect_failure 'zero Z process-grid dimension' \
    'Process-grid dimensions must be positive:' 1 1 0
expect_failure 'process-grid and world-size mismatch' \
    'Process-grid size mismatch:' 2 1 1

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d parallelism CLI checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d parallelism CLI checks)\n' "$failures" "$checks"
exit 1
