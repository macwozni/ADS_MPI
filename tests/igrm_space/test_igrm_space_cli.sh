#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./igrm_space_probe}
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

    if [[ $status -eq 0 ]]; then
        printf 'FAIL %s (unexpected success)\n' "$label"
        failures=$((failures + 1))
    elif [[ $status -ne 1 ]]; then
        printf 'FAIL %s (unexpected status %d)\n' "$label" "$status"
        sed 's/^/  /' "$stderr_file"
        failures=$((failures + 1))
    elif [[ -s "$stdout_file" ]]; then
        printf 'FAIL %s (unexpected stdout)\n' "$label"
        sed 's/^/  /' "$stdout_file"
        failures=$((failures + 1))
    elif grep -Fq -- "$expected_message" "$stderr_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (status %d, missing diagnostic)\n' "$label" "$status"
        sed 's/^/  /' "$stderr_file"
        failures=$((failures + 1))
    fi
}

expect_failure 'equal degree, continuity-only enrichment' equal-degree \
    'iGRM requires test degree greater than trial degree:'
expect_failure 'test degree lower than trial degree' lower-degree \
    'iGRM requires test degree greater than trial degree:'
expect_failure 'knot displacement above tolerance' shifted-location \
    'iGRM requires identical test/trial knot locations'
expect_failure 'additional test knot location' extra-test-location \
    'iGRM requires identical test/trial knot locations'
expect_failure 'additional trial knot location' extra-trial-location \
    'iGRM requires identical test/trial knot locations'
expect_failure 'element-count metadata mismatch' metadata-mismatch \
    'iGRM mesh metadata mismatch:'

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d iGRM-space CLI checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d iGRM-space CLI checks)\n' "$failures" "$checks"
exit 1
