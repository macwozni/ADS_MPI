#!/usr/bin/env bash
set -u

probe=${1:?missing fault-probe executable}
contract_probe=${2:?missing error-contract probe target}
limit=${FAULT_TIMEOUT:-5s}
overall=0

run_cases() {
    local executable=$1
    local kind=$2
    local mode result

    for mode in irecv isend waitall; do
        if timeout "$limit" "$executable" "$mode"; then
            :
        else
            result=$?
            if [ "$result" -eq 124 ]; then
                echo "FAIL: $kind $mode probe exceeded timeout $limit"
            fi
            overall=1
        fi
    done
}

run_cases "$probe" behavior

build_log=$(mktemp)
trap 'rm -f "$build_log"' EXIT
if make --no-print-directory "$contract_probe" >"$build_log" 2>&1; then
    run_cases "./$contract_probe" error-contract
else
    echo 'FAIL: DistributeSpline has no caller-visible MPI error contract'
    tail -n 8 "$build_log" | sed 's/^/  /'
    overall=1
fi

exit "$overall"
