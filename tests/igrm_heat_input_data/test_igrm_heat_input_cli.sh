#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./igrm_heat_input_probe}
case_count=0

fail() {
    printf 'FAIL: %s\n' "$1" >&2
    exit 1
}

assert_output() {
    local expected=$1
    shift
    local actual

    actual=$("$@") || fail "command failed unexpectedly: $*"
    [[ "$actual" == "$expected" ]] || \
        fail "expected [$expected], got [$actual]"
    case_count=$((case_count + 1))
}

assert_failure() {
    local expected_message=$1
    shift
    local actual status

    set +e
    actual=$("$@" 2>&1)
    status=$?
    set -e

    [[ $status -eq 5 ]] || fail "expected exit status 5, got $status: $*"
    [[ "$actual" == *"$expected_message"* ]] || \
        fail "missing diagnostic [$expected_message] in [$actual]"
    case_count=$((case_count + 1))
}

valid=(4 5 6 3 4 5 2 2 3 1 2 3 0 0.125)
expected_prefix="4 5 6 3 4 5 2 2 3 1 2 3 0 125000"

assert_output "$expected_prefix 1 dg" "$probe" "${valid[@]}"
assert_output "2 3 4 1 2 3 0 1 2 3 2 1 7 1500000 1 dg" \
    "$probe" 2 3 4 1 2 3 0 1 2 3 2 1 7 1.5 DG
assert_output "$expected_prefix 1 dg" \
    "$probe" "${valid[@]}" Douglas_Gunn
assert_output "$expected_prefix 2 pr" \
    "$probe" "${valid[@]}" PEACEMAN-RACHFORD
assert_output "$expected_prefix 3 be" \
    "$probe" "${valid[@]}" BackwardEuler

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${valid[@]}" dg extra
assert_failure "invalid integer argument" \
    "$probe" 4 5 invalid 3 4 5 2 2 3 1 2 3 0 0.125
assert_failure "invalid integer argument" \
    "$probe" 4 5 "6 junk" 3 4 5 2 2 3 1 2 3 0 0.125
assert_failure "invalid real argument" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 3 0 invalid
assert_failure "invalid real argument" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 3 0 "0.125 junk"

assert_failure "number of elements must be positive" \
    "$probe" 0 5 6 3 4 5 2 2 3 1 2 3 0 0.125
assert_failure "number of elements must be positive" \
    "$probe" 4 -5 6 3 4 5 2 2 3 1 2 3 0 0.125
assert_failure "number of elements must be positive" \
    "$probe" 4 5 0 3 4 5 2 2 3 1 2 3 0 0.125
assert_failure "test polynomial degree must be non-negative" \
    "$probe" 4 5 6 -1 4 5 2 2 3 1 2 3 0 0.125
assert_failure "test polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 -1 5 2 2 3 1 2 3 0 0.125
assert_failure "test polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 4 -1 2 2 3 1 2 3 0 0.125
assert_failure "trial polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 4 5 -1 2 3 1 2 3 0 0.125
assert_failure "trial polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 4 5 2 -1 3 1 2 3 0 0.125
assert_failure "trial polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 4 5 2 2 -1 1 2 3 0 0.125

assert_failure "test polynomial degree must be greater than trial polynomial degree" \
    "$probe" 4 5 6 2 4 5 2 2 3 1 2 3 0 0.125
assert_failure "test polynomial degree must be greater than trial polynomial degree" \
    "$probe" 4 5 6 3 2 5 2 3 3 1 2 3 0 0.125
assert_failure "test polynomial degree must be greater than trial polynomial degree" \
    "$probe" 4 5 6 3 4 3 2 2 3 1 2 3 0 0.125
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 5 6 3 4 5 2 2 3 1 2 3 0 0.125

assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 0 2 3 0 0.125
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 -2 3 0 0.125
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 0 0 0.125
assert_failure "number of time steps must be non-negative" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 3 -1 0.125
assert_failure "time step must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 3 0 0
assert_failure "time step must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 3 0 -0.125
assert_failure "forward euler is not an iGRM time-scheme option" \
    "$probe" "${valid[@]}" forward_euler
assert_failure "unknown time scheme" "$probe" "${valid[@]}" crank-nicolson

printf 'OK (%d command-line cases)\n' "$case_count"
