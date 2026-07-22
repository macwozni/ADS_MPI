#!/usr/bin/env bash

set -euo pipefail

probe=${1:-./argument_parser_probe}
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

assert_output "raw:Alpha-Beta_42" "$probe" raw "Alpha-Beta_42"
assert_output "raw:  Mixed Value" "$probe" raw "  Mixed Value  "
assert_output "integer:12345" "$probe" integer "12345"
assert_output "integer:12345" "$probe" integer "+12345"
assert_output "integer:-6789" "$probe" integer "-6789"
assert_output "integer:2147483647" "$probe" integer "2147483647"
assert_output "integer:-2147483648" "$probe" integer "-2147483648"
assert_output "real:ok" "$probe" real "-12.5" "negative-decimal"
assert_output "real:ok" "$probe" real "6.25d-3" "d-exponent"
assert_output "real:ok" "$probe" real "+6.25D-3" "d-exponent"
assert_output "real:ok" "$probe" real ".5" "leading-dot"
assert_output "real:ok" "$probe" real "5." "trailing-dot"
assert_output "real:ok" "$probe" real "1" "integer-form"
assert_output "real:ok" "$probe" real "1e+3" "e-exponent"
assert_output "positive-integer:ok" "$probe" positive-integer "1"
assert_output "nonnegative-integer:ok" "$probe" nonnegative-integer "0"
assert_output "positive-real:ok" "$probe" positive-real "0.5"
assert_output "safe-spline:ok" "$probe" safe-spline 1 0
assert_output "safe-spline:ok" "$probe" safe-spline 2147483646 0
assert_output "safe-spline:ok" "$probe" safe-spline 1 1073741822
assert_output "safe-spline:ok" "$probe" safe-spline 2147483644 0 1
assert_output "string:mixed_value" "$probe" string "  MiXeD_Value  "
assert_output "string:" "$probe" string ""
assert_output "scheme:2:pr" "$probe" scheme "  PeAcEmAn_RaChFoRd  "
assert_output "scheme:3:be" "$probe" scheme "BackwardEuler"

exact_length_argument="abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuv"
[[ ${#exact_length_argument} -eq 100 ]] || fail "boundary fixture is not exactly 100 characters"
assert_output "raw:$exact_length_argument" "$probe" raw "$exact_length_argument"

assert_failure "invalid command argument" "$probe" raw
assert_failure "invalid integer argument" "$probe" integer "12x"
assert_failure "invalid integer argument" "$probe" integer "12 junk"
assert_failure "invalid integer argument" "$probe" integer "12 34"
assert_failure "invalid integer argument" "$probe" integer "12,34"
assert_failure "invalid integer argument" "$probe" integer "12,"
assert_failure "invalid integer argument" "$probe" integer "12 /"
assert_failure "invalid integer argument" "$probe" integer "2147483648"
assert_failure "invalid command argument" "$probe" integer
assert_failure "invalid real argument" "$probe" real "not-a-real"
assert_failure "invalid real argument" "$probe" real "1.0 junk" unused
assert_failure "invalid real argument" "$probe" real "1.0 2.0" unused
assert_failure "invalid real argument" "$probe" real "1.0,2.0" unused
assert_failure "invalid real argument" "$probe" real "1.0," unused
assert_failure "invalid real argument" "$probe" real "1.0 /" unused
assert_failure "invalid real argument" "$probe" real "NaN" unused
assert_failure "invalid real argument" "$probe" real "Inf" unused
assert_failure "invalid real argument" "$probe" real "Infinity" unused
assert_failure "probe integer must be positive" "$probe" positive-integer "0"
assert_failure "probe integer must be positive" "$probe" positive-integer "-1"
assert_failure "probe integer must be non-negative" "$probe" nonnegative-integer "-1"
assert_failure "probe real must be positive" "$probe" positive-real "0"
assert_failure "probe real must be positive" "$probe" positive-real "-0.5"
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" safe-spline 2147483647 0
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" safe-spline 1 1073741823
assert_failure "spline dimensions exceed supported integer range" \
    "$probe" safe-spline 2147483645 0 1
for forward_euler_alias in fe forward-euler Forward_Euler ForwardEuler; do
    assert_failure "forward euler is not an iGRM time-scheme option" \
        "$probe" scheme "$forward_euler_alias"
done
assert_failure "unknown time scheme" "$probe" scheme "crank-nicolson"
assert_failure "unknown time scheme" "$probe" scheme ""
assert_failure "invalid command argument" "$probe" scheme

long_argument="${exact_length_argument}x"
[[ ${#long_argument} -gt 100 ]] || fail "long-argument fixture is not longer than 100 characters"
assert_failure "invalid command argument" "$probe" raw "$long_argument"

printf 'OK (%d command-line cases)\n' "$case_count"
