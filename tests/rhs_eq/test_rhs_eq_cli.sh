#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./rhs_eq_probe}
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
    local expected_message=$3
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
    elif grep -Fq -- "$expected_message" "$stderr_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (missing diagnostic)\n' "$label"
        sed 's/^/  /' "$stderr_file"
        failures=$((failures + 1))
    fi
}

expect_failure 'substep below supported range' substep-zero 'wrong substep: 0'
expect_failure 'substep above supported range' substep-four 'wrong substep: 4'
expect_failure 'negative derivative-state selector' state-negative \
    'wrong RHS derivative state: -1, term: 1, substep: 1'
expect_failure 'derivative-state selector above supported range' state-four \
    'wrong RHS derivative state: 4, term: 6, substep: 1'

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d RHS_eq CLI checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d RHS_eq CLI checks)\n' "$failures" "$checks"
exit 1
