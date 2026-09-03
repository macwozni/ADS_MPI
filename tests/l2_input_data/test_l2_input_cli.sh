#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./l2_input_probe}
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

valid=(4 5 6 2 1 2 3)
assert_output "4 5 6 2 1 2 3" "$probe" "${valid[@]}"
assert_output "1 2 3 0 4 5 6" "$probe" 1 2 3 0 4 5 6

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${valid[@]}" 8
assert_failure "invalid integer argument" "$probe" 4 5 invalid 2 1 2 3
assert_failure "number of elements must be positive" "$probe" 0 5 6 2 1 2 3
assert_failure "number of elements must be positive" "$probe" 4 -5 6 2 1 2 3
assert_failure "number of elements must be positive" "$probe" 4 5 0 2 1 2 3
assert_failure "polynomial order must be non-negative" "$probe" 4 5 6 -1 1 2 3
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 5 6 0 1 2 3
assert_failure "process-grid dimension must be positive" "$probe" 4 5 6 2 0 2 3
assert_failure "process-grid dimension must be positive" "$probe" 4 5 6 2 1 -2 3
assert_failure "process-grid dimension must be positive" "$probe" 4 5 6 2 1 2 0

printf 'OK (%d command-line cases)\n' "$case_count"
