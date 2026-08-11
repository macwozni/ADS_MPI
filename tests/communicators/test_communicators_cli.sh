#!/usr/bin/env bash

set -u

PROBE=${PROBE:-./communicators_error_probe}
checks=0
failures=0

expect_failure() {
    local label=$1
    local expected_message=$2
    local mode=$3
    local output_file status

    output_file=$(mktemp)
    "$PROBE" "$mode" >"$output_file" 2>&1
    status=$?
    checks=$((checks + 1))

    if [[ $status -ne 0 ]] && grep -Fq -- "$expected_message" "$output_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (status %d)\n' "$label" "$status"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    fi
    rm -f -- "$output_file"
}

expect_success() {
    local label=$1
    local expected_message=$2
    local mode=$3
    local output_file status

    output_file=$(mktemp)
    "$PROBE" "$mode" >"$output_file" 2>&1
    status=$?
    checks=$((checks + 1))

    if [[ $status -eq 0 ]] && grep -Fq -- "$expected_message" "$output_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (status %d)\n' "$label" "$status"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    fi
    rm -f -- "$output_file"
}

expect_failure 'mpi_comm_group error is fatal' \
    'error calling mpi_comm_group' comm-group
expect_failure 'Z mpi_group_incl error is fatal' \
    'error calling mpi_group_incl for Z' group-z
expect_failure 'Y mpi_group_incl error is fatal' \
    'error calling mpi_group_incl for Y' group-y
expect_failure 'X mpi_group_incl error is fatal' \
    'error calling mpi_group_incl for X' group-x
expect_failure 'world-group release error is fatal' \
    'error freeing MPI_COMM_WORLD group' world-group-free
expect_failure 'Z mpi_comm_create error is fatal' \
    'error calling mpi_com_create for Z' comm-z
expect_failure 'Y mpi_comm_create error is fatal' \
    'error calling mpi_com_create for Y' comm-y
expect_failure 'X mpi_comm_create error is fatal' \
    'error calling mpi_com_create for X' comm-x
expect_success 'cleanup retains first error and attempts every release' \
    'SUCCESS cleanup retained first error and continued' cleanup

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d communicator stub checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d communicator stub checks)\n' "$failures" "$checks"
exit 1
