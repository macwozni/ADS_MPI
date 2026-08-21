#!/usr/bin/env bash

set -u

MPIEXEC=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
EXEC_DIR=${EXEC_DIR:-../../mymake/EXEC}
DRIVER_SMOKE_TIMEOUT=${DRIVER_SMOKE_TIMEOUT:-60s}
ASAN_OPTIONS="${ASAN_OPTIONS:+${ASAN_OPTIONS}:}detect_leaks=0"
export ASAN_OPTIONS

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/mpi-ads-driver-smoke.XXXXXX")
checks=0
failures=0

cleanup() {
    case "$work_dir" in
        "${TMPDIR:-/tmp}"/mpi-ads-driver-smoke.*)
            rm -rf -- "$work_dir"
            ;;
    esac
}
trap cleanup EXIT

run_smoke() {
    local label=$1
    local driver=$2
    shift 2

    local output_file="$work_dir/${label}.log"
    local status

    (
        cd "$work_dir" || exit 1
        OMP_DYNAMIC=FALSE OMP_NUM_THREADS=2 \
            timeout "$DRIVER_SMOKE_TIMEOUT" "$MPIEXEC" -n 1 "$EXEC_DIR/$driver" "$@"
    ) >"$output_file" 2>&1
    status=$?
    checks=$((checks + 1))

    if [[ $status -eq 0 ]]; then
        printf 'PASS %s positive smoke\n' "$label"
    elif [[ $status -eq 124 ]]; then
        printf 'FAIL %s positive smoke (timeout)\n' "$label"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    else
        printf 'FAIL %s positive smoke (status %d)\n' "$label" "$status"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    fi
}

run_smoke L2 l2 2 2 2 1 1 1 1
run_smoke heat heat 2 1 1 0.01 1 1 1
run_smoke Eriksson eriksson 2 1 1 0.01 1 1 1
run_smoke pure-diffusion-iGRM pure_diffusion_igrm 3 1 1 1 1 1 0.1 dg
run_smoke oil oil \
    2 1 1 1 1 1 0.1 \
    1 0.5 0.5 0.5 \
    1 0.25 0.25 0.25
run_smoke iGRM-L2 igrm_l2 \
    2 2 2 3 3 3 1 1 1 \
    1 1 1 pr
run_smoke iGRM-heat igrm_heat \
    2 2 2 3 3 3 2 2 2 \
    1 1 1 1 0.001 dg
run_smoke iGRM-Eriksson-MUMPS igrm_eirksson \
    2 2 2 2 2 2 1 1 1 \
    1 1 1
run_smoke iGRM-Stokes igrm_stokes \
    2 2 2 2 2 2 2 2 2 \
    1 1 1

if [[ $failures -ne 0 ]]; then
    printf 'FAILED (%d/%d positive driver smoke checks)\n' "$failures" "$checks"
    exit 1
fi

printf 'OK (%d positive driver smoke checks)\n' "$checks"
