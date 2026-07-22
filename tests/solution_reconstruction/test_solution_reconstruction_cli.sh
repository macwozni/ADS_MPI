#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./solution_reconstruction_probe}
stdout_file=$(mktemp)
stderr_file=$(mktemp)
checks=0
failures=0

cleanup() {
    rm -f -- "$stdout_file" "$stderr_file"
}
trap cleanup EXIT

expect_failure() {
    local label=$1
    local mode=$2
    local status

    : >"$stdout_file"
    : >"$stderr_file"
    if "$probe" "$mode" >"$stdout_file" 2>"$stderr_file"; then
        status=0
    else
        status=$?
    fi
    checks=$((checks + 1))

    if [[ $status -ne 1 ]]; then
        printf 'FAIL %s (status %d, expected 1)\n' "$label" "$status"
        sed 's/^/  /' "$stderr_file"
        failures=$((failures + 1))
    elif [[ -s "$stdout_file" ]]; then
        printf 'FAIL %s (unexpected stdout)\n' "$label"
        sed 's/^/  /' "$stdout_file"
        failures=$((failures + 1))
    elif grep -Fq -- 'wrong substep' "$stderr_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (missing diagnostic)\n' "$label"
        sed 's/^/  /' "$stderr_file"
        failures=$((failures + 1))
    fi
}

expect_failure 'selector below supported range' zero
expect_failure 'selector above supported range' four

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d solution-reconstruction CLI checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d solution-reconstruction CLI checks)\n' "$failures" "$checks"
exit 1
