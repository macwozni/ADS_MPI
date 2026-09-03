#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./igrm_pollution_input_probe}
case_count=0

fail() {
    printf 'FAIL: %s\n' "$1" >&2
    exit 1
}

assert_output() {
    local expected=$1
    shift
    local actual

    actual=$(env -u ADS_POLLUTION_OUTPUT_RESOLUTION "$@") || \
        fail "command failed unexpectedly: $*"
    [[ "$actual" == "$expected" ]] || \
        fail "expected [$expected], got [$actual]"
    case_count=$((case_count + 1))
}

assert_env_output() {
    local resolution=$1
    local expected=$2
    shift 2
    local actual

    actual=$(env ADS_POLLUTION_OUTPUT_RESOLUTION="$resolution" "$@") || \
        fail "command failed unexpectedly with resolution $resolution: $*"
    [[ "$actual" == "$expected" ]] || \
        fail "expected [$expected], got [$actual]"
    case_count=$((case_count + 1))
}

assert_failure() {
    local expected_message=$1
    shift
    local actual status

    set +e
    actual=$(env -u ADS_POLLUTION_OUTPUT_RESOLUTION "$@" 2>&1)
    status=$?
    set -e

    [[ $status -eq 5 ]] || fail "expected exit status 5, got $status: $*"
    [[ "$actual" == *"$expected_message"* ]] || \
        fail "missing diagnostic [$expected_message] in [$actual]"
    case_count=$((case_count + 1))
}

assert_env_failure() {
    local resolution=$1
    local expected_message=$2
    shift 2
    local actual status

    set +e
    actual=$(env ADS_POLLUTION_OUTPUT_RESOLUTION="$resolution" "$@" 2>&1)
    status=$?
    set -e

    [[ $status -eq 5 ]] || \
        fail "expected exit status 5 for resolution $resolution, got $status"
    [[ "$actual" == *"$expected_message"* ]] || \
        fail "missing diagnostic [$expected_message] in [$actual]"
    case_count=$((case_count + 1))
}

reference=(4 0 2 1 3 0 2 1 1 1)
assert_output "4 0 2 1 3 0 2 1 1 1 100" "$probe" "${reference[@]}"
assert_output "4 1 3 2 3 1 7 1 2 3 100" \
    "$probe" 4 1 3 2 3 1 7 1 2 3
assert_output "5 0 3 1 3 1 1 2 1 1 100" \
    "$probe" 5 0 3 1 3 1 1 2 1 1
assert_env_output 4 "4 0 2 1 3 0 2 1 1 1 4" \
    "$probe" "${reference[@]}"
assert_env_output +7 "4 0 2 1 3 0 2 1 1 1 7" \
    "$probe" "${reference[@]}"

assert_failure "proper usage with arguments" "$probe"
assert_failure "proper usage with arguments" "$probe" "${reference[@]}" 99
assert_failure "invalid integer argument" \
    "$probe" 4 bad 2 1 3 0 2 1 1 1
assert_failure "number of elements must be positive" \
    "$probe" 0 0 2 1 3 0 2 1 1 1
assert_failure "adapt flag must be zero or one" \
    "$probe" 4 -1 2 1 3 0 2 1 1 1
assert_failure "adapt flag must be zero or one" \
    "$probe" 4 2 2 1 3 0 2 1 1 1
assert_failure "trial polynomial degree must be positive" \
    "$probe" 4 0 0 0 3 0 2 1 1 1
assert_failure "test polynomial degree must be positive" \
    "$probe" 4 0 2 1 0 0 2 1 1 1
assert_failure "unsupported polynomial degree in iGRM trial space" \
    "$probe" 4 0 10 1 3 0 2 1 1 1
assert_failure "unsupported polynomial degree in iGRM test space" \
    "$probe" 4 0 2 1 10 0 2 1 1 1
assert_failure \
    "trial continuity must be between zero and polynomial degree minus one" \
    "$probe" 4 0 2 -1 3 0 2 1 1 1
assert_failure \
    "trial continuity must be between zero and polynomial degree minus one" \
    "$probe" 4 0 2 2 3 0 2 1 1 1
assert_failure \
    "test continuity must be between zero and polynomial degree minus one" \
    "$probe" 4 0 2 1 3 -1 2 1 1 1
assert_failure \
    "test continuity must be between zero and polynomial degree minus one" \
    "$probe" 4 0 2 1 3 3 2 1 1 1
assert_failure "number of time steps must be positive" \
    "$probe" 4 0 2 1 3 0 0 1 1 1
assert_failure "number of time steps must be positive" \
    "$probe" 4 0 2 1 3 0 -3 1 1 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 0 2 1 3 0 2 0 1 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 0 2 1 3 0 2 1 -1 1
assert_failure "process-grid dimension must be positive" \
    "$probe" 4 0 2 1 3 0 2 1 1 0
assert_failure \
    "test-space dimension must be greater than or equal to trial-space dimension" \
    "$probe" 4 0 3 0 3 2 2 1 1 1
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" 2147483647 0 9 0 9 0 2 1 1 1

assert_env_failure 0 \
    "ADS_POLLUTION_OUTPUT_RESOLUTION must be a positive integer" \
    "$probe" "${reference[@]}"
assert_env_failure -2 \
    "ADS_POLLUTION_OUTPUT_RESOLUTION must be a positive integer" \
    "$probe" "${reference[@]}"
assert_env_failure invalid \
    "ADS_POLLUTION_OUTPUT_RESOLUTION must be a positive integer" \
    "$probe" "${reference[@]}"
assert_env_failure 2147483648 \
    "ADS_POLLUTION_OUTPUT_RESOLUTION must be a positive integer" \
    "$probe" "${reference[@]}"

printf 'OK (%d command-line cases)\n' "$case_count"
