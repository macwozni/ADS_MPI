#!/usr/bin/env bash
set -u

probe=${1:?missing solution-output probe}
contract_target=${2:?missing error-contract target}
limit=${FAULT_TIMEOUT:-10s}
overall=0

timeout "$limit" "$probe" mpi-failure
result=$?
if [ "$result" -ne 0 ]; then
    if [ "$result" -eq 124 ]; then
        echo "FAIL: solution-output fault probe exceeded timeout $limit"
    fi
    overall=1
fi

build_log=$(mktemp)
trap 'rm -f "$build_log"' EXIT
if make --no-print-directory "$contract_target" >"$build_log" 2>&1; then
    timeout "$limit" "./$contract_target"
    result=$?
    if [ "$result" -ne 0 ]; then
        if [ "$result" -eq 124 ]; then
            echo "FAIL: solution-output error-contract probe exceeded timeout $limit"
        fi
        overall=1
    fi
else
    echo 'FAIL: PrintSolution has no caller-visible MPI error contract'
    tail -n 8 "$build_log" | sed 's/^/  /'
    overall=1
fi

exit "$overall"
