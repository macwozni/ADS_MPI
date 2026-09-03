#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./pure_diffusion_igrm_input_probe}
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
    [[ "$actual" == "$expected" ]] || fail "expected [$expected], got [$actual]"
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

reference=(8 3 1 2 3 4 0.125)
assert_output "8 3 1 2 3 4 125000 dg" "$probe" "${reference[@]}"
assert_output "8 3 1 2 3 4 125000 pr" "$probe" "${reference[@]}" PR
assert_output "8 3 1 2 3 4 125000 be" \
    "$probe" "${reference[@]}" backward_euler
assert_output "1 0 2 1 4 0 1000000 dg" "$probe" 1 0 2 1 4 0 1.0 dg

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${reference[@]}" dg extra
assert_failure "invalid integer argument" "$probe" 8 bad 1 2 3 4 0.125
assert_failure "invalid real argument" "$probe" 8 3 1 2 3 4 invalid
assert_failure "number of elements must be positive" "$probe" 0 3 1 2 3 4 0.125
assert_failure "polynomial order must be non-negative" "$probe" 8 -1 1 2 3 4 0.125
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 0 1 2 3 4 0.125
assert_failure "number of time steps must be non-negative" \
    "$probe" 8 3 1 2 3 -1 0.125
assert_failure "time step must be positive" "$probe" 8 3 1 2 3 4 0
assert_failure "process-grid dimension must be positive" \
    "$probe" 8 3 0 2 3 4 0.125
assert_failure "process-grid dimension must be positive" \
    "$probe" 8 3 1 -2 3 4 0.125
assert_failure "process-grid dimension must be positive" \
    "$probe" 8 3 1 2 0 4 0.125
assert_failure "forward euler is not an iGRM time-scheme option" \
    "$probe" "${reference[@]}" fe
assert_failure "unknown time scheme" "$probe" "${reference[@]}" crank-nicolson

printf 'OK (%d command-line cases)\n' "$case_count"
