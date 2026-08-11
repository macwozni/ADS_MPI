#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./int2str_probe}
stdout_file=$(mktemp)
stderr_file=$(mktemp)
case_count=0

cleanup() {
    rm -f -- "$stdout_file" "$stderr_file"
}
trap cleanup EXIT

fail() {
    printf 'FAIL: %s\n' "$1" >&2
    exit 1
}

run_probe() {
    : >"$stdout_file"
    : >"$stderr_file"
    "$probe" "$1" "$2" >"$stdout_file" 2>"$stderr_file" || \
        fail "probe failed for value $1 and length $2"
}

assert_no_warning() {
    local number=$1
    local length=$2
    local expected_output=$3
    local actual_output

    run_probe "$number" "$length"
    actual_output=$(<"$stdout_file")

    [[ "$actual_output" == "$expected_output" ]] || \
        fail "expected stdout [$expected_output], got [$actual_output]"
    [[ ! -s "$stderr_file" ]] || \
        fail "unexpected warning: $(<"$stderr_file")"
    case_count=$((case_count + 1))
}

assert_warning() {
    local number=$1
    local length=$2
    local expected_output=$3
    local expected_warning=$4
    local actual_output actual_warning

    run_probe "$number" "$length"
    actual_output=$(<"$stdout_file")
    actual_warning=$(<"$stderr_file")

    [[ "$actual_output" == "$expected_output" ]] || \
        fail "expected stdout [$expected_output], got [$actual_output]"
    [[ "$actual_warning" == "$expected_warning" ]] || \
        fail "expected warning [$expected_warning], got [$actual_warning]"
    case_count=$((case_count + 1))
}

assert_no_warning 0 4 "[0   ]"
assert_no_warning 42 6 "[42    ]"
assert_no_warning 123 3 "[123]"
assert_no_warning -123 4 "[-123]"

assert_warning 123 2 "[12]" \
    "int2str: WARNING: can't fit 123 into a   2-character variable"
assert_warning -123 3 "[-12]" \
    "int2str: WARNING: can't fit -123 into a   3-character variable"

printf 'OK (%d int2str CLI cases)\n' "$case_count"
