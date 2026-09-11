#!/usr/bin/env bash
set -u

probe=${1:?missing synchronization-error probe}
limit=${SYNC_CASE_TIMEOUT:-10s}
overall=0

for test_case in initial-allreduce solver-allreduce solver-bcast; do
    if timeout "$limit" "$probe" "$test_case"; then
        :
    else
        result=$?
        if [ "$result" -eq 124 ]; then
            echo "FAIL: $test_case synchronization probe exceeded timeout $limit"
        fi
        overall=1
    fi
done

exit "$overall"
