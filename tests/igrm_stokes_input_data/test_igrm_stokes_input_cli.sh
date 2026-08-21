#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./igrm_stokes_input_probe}
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

reference=(4 4 4 2 2 2 2 2 2 1 1 1)
assert_output "4 4 4 2 2 2 2 2 2 1 1 1" "$probe" "${reference[@]}"
assert_output "4 5 6 4 3 2 3 2 1 1 2 3" \
    "$probe" 4 5 6 4 3 2 3 2 1 1 2 3

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${reference[@]}" 99
assert_failure "invalid integer argument" \
    "$probe" 4 4 4 2 invalid 2 2 2 2 1 1 1
assert_failure "number of elements must be positive" \
    "$probe" 0 4 4 2 2 2 2 2 2 1 1 1
assert_failure "number of elements must be positive" \
    "$probe" 4 -1 4 2 2 2 2 2 2 1 1 1
assert_failure "test polynomial degree must be positive" \
    "$probe" 4 4 4 0 2 2 2 2 2 1 1 1
assert_failure "trial polynomial degree must be positive" \
    "$probe" 4 4 4 2 2 2 2 0 2 1 1 1
assert_failure "unsupported polynomial degree in iGRM test space" \
    "$probe" 4 4 4 2 10 2 2 2 2 1 1 1
assert_failure "unsupported polynomial degree in iGRM trial space" \
    "$probe" 4 4 4 9 9 9 2 2 10 1 1 1
assert_failure \
    "test polynomial degree must be greater than or equal to trial polynomial degree" \
    "$probe" 4 4 4 1 2 2 2 2 2 1 1 1
assert_failure \
    "test polynomial degree must be greater than or equal to trial polynomial degree" \
    "$probe" 4 4 4 2 2 2 2 3 2 1 1 1
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 4 4 2 2 2 2 2 2 1 1 1
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 4 4 2147483647 2 2 2 2 2 2 1 1 1
assert_failure \
    "discontinuous test-space dimensions exceed supported integer range" \
    "$probe" 214748365 4 4 9 2 2 1 1 1 1 1 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 4 4 2 2 2 2 2 2 0 1 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 4 4 2 2 2 2 2 2 1 -2 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 4 4 2 2 2 2 2 2 1 1 0

printf 'OK (%d command-line cases)\n' "$case_count"
