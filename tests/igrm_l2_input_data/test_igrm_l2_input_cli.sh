#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./igrm_l2_input_probe}
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

base=(4 5 6 3 4 5 2 2 3 1 2 3)
assert_output "4 5 6 3 4 5 2 2 3 1 2 3 1 dg 1.000000" \
    "$probe" "${base[@]}"
assert_output "4 5 6 3 4 5 2 2 3 1 2 3 2 pr 1.000000" \
    "$probe" "${base[@]}" PR
assert_output "4 5 6 3 4 5 2 2 3 1 2 3 3 be 1.000000" \
    "$probe" "${base[@]}" backward-euler
assert_output "4 5 6 3 4 5 2 2 3 1 2 3 1 dg .250000" \
    "$probe" "${base[@]}" 0.25 Douglas_Gunn

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${base[@]:0:11}"
assert_failure "proper usage with arguments" "$probe" "${base[@]}" 1.0 dg extra
assert_failure "invalid integer argument" \
    "$probe" 4 5 6 3 invalid 5 2 2 3 1 2 3
assert_failure "number of elements must be positive" \
    "$probe" 0 5 6 3 4 5 2 2 3 1 2 3
assert_failure "number of elements must be positive" \
    "$probe" 4 -5 6 3 4 5 2 2 3 1 2 3
assert_failure "number of elements must be positive" \
    "$probe" 4 5 0 3 4 5 2 2 3 1 2 3
assert_failure "test polynomial degree must be non-negative" \
    "$probe" 4 5 6 -1 4 5 2 2 3 1 2 3
assert_failure "trial polynomial degree must be non-negative" \
    "$probe" 4 5 6 3 4 5 2 -1 3 1 2 3
assert_failure "test polynomial degree must be greater than trial polynomial degree" \
    "$probe" 4 5 6 2 4 5 2 2 3 1 2 3
assert_failure "test polynomial degree must be greater than trial polynomial degree" \
    "$probe" 4 5 6 3 2 5 2 3 3 1 2 3
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 5 6 3 4 5 2 2 3 1 2 3
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 0 2 3
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 -2 3
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 5 6 3 4 5 2 2 3 1 2 0
assert_failure "invalid real argument" \
    "$probe" "${base[@]}" invalid dg
assert_failure "time-scheme tau must be positive" \
    "$probe" "${base[@]}" 0 dg
assert_failure "unknown time scheme" "$probe" "${base[@]}" crank-nicolson
assert_failure "forward euler is not an iGRM time-scheme option" \
    "$probe" "${base[@]}" 1.0 forward-euler

printf 'OK (%d command-line cases)\n' "$case_count"
