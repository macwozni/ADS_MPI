#!/usr/bin/env bash
set -eu

probe=${1:-./time_scheme_validation_probe}
case_count=0

expect_failure() {
    label=$1
    expected=$2
    pattern=$3
    shift 3

    stdout_file=$(mktemp)
    stderr_file=$(mktemp)
    set +e
    "$probe" "$@" >"$stdout_file" 2>"$stderr_file"
    status=$?
    set -e

    if [ "$status" -ne "$expected" ]; then
        printf 'FAIL: %s returned %d, expected %d\n' "$label" "$status" "$expected" >&2
        sed -n '1,20p' "$stdout_file" >&2
        sed -n '1,20p' "$stderr_file" >&2
        rm -f "$stdout_file" "$stderr_file"
        exit 1
    fi
    if ! grep -F "$pattern" "$stderr_file" "$stdout_file" >/dev/null; then
        printf 'FAIL: %s did not contain: %s\n' "$label" "$pattern" >&2
        sed -n '1,20p' "$stdout_file" >&2
        sed -n '1,20p' "$stderr_file" >&2
        rm -f "$stdout_file" "$stderr_file"
        exit 1
    fi

    rm -f "$stdout_file" "$stderr_file"
    case_count=$((case_count + 1))
}

for axis_spec in '1 x' '2 y' '3 z'; do
    set -- $axis_spec
    axis=$1
    axis_name=$2
    expect_failure "standard ng=0 axis $axis_name" 5 \
        "unsupported Gauss quadrature size in ADS space along $axis_name" \
        standard "$axis" ng 0
    expect_failure "standard ng=11 axis $axis_name" 5 \
        "unsupported Gauss quadrature size in ADS space along $axis_name" \
        standard "$axis" ng 11
done

expect_failure 'negative standard degree' 5 \
    'unsupported polynomial degree in ADS space along x' standard 1 p -1
expect_failure 'too large standard degree' 5 \
    'unsupported polynomial degree in ADS space along z' standard 3 p 10
expect_failure 'invalid iGRM test quadrature' 5 \
    'unsupported Gauss quadrature size in iGRM test space along y' test 2 ng 0
expect_failure 'invalid iGRM trial quadrature' 5 \
    'unsupported Gauss quadrature size in iGRM trial space along z' trial 3 ng 11
expect_failure 'invalid iGRM test degree' 5 \
    'unsupported polynomial degree in iGRM test space along x' test 1 p 10
expect_failure 'invalid iGRM trial degree' 5 \
    'unsupported polynomial degree in iGRM trial space along y' trial 2 p -1

printf 'OK (%d time-scheme validation failures)\n' "$case_count"
