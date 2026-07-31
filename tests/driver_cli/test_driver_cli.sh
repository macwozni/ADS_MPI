#!/usr/bin/env bash

set -u

MPIEXEC=${MPIEXEC:-/opt/lib/mpich-5.0.0/bin/mpiexec}
EXEC_DIR=${EXEC_DIR:-../../mymake/EXEC}
ASAN_OPTIONS="${ASAN_OPTIONS:+${ASAN_OPTIONS}:}detect_leaks=0"
export ASAN_OPTIONS

checks=0
failures=0
temp_files=()

cleanup() {
    rm -f -- "${temp_files[@]}"
}
trap cleanup EXIT

expect_failure() {
    local label=$1
    local ranks=$2
    local driver=$3
    local expected_message=$4
    shift 4

    local output_file status
    output_file=$(mktemp)
    temp_files+=("$output_file")

    if [[ $ranks -eq 1 ]]; then
        timeout 20s "$EXEC_DIR/$driver" "$@" >"$output_file" 2>&1
    else
        timeout 20s "$MPIEXEC" -n "$ranks" "$EXEC_DIR/$driver" "$@" \
            >"$output_file" 2>&1
    fi
    status=$?
    checks=$((checks + 1))

    if [[ $status -eq 124 ]]; then
        printf 'FAIL %s (timeout)\n' "$label"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    elif [[ $status -eq 0 ]]; then
        printf 'FAIL %s (unexpected success)\n' "$label"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    elif grep -Fq -- "$expected_message" "$output_file"; then
        printf 'PASS %s\n' "$label"
    else
        printf 'FAIL %s (status %d, missing diagnostic)\n' "$label" "$status"
        sed 's/^/  /' "$output_file"
        failures=$((failures + 1))
    fi

    rm -f -- "$output_file"
}

# Each real main program must delegate an invalid argument count to its
# production input_data module before entering the numerical workflow.
expect_failure 'L2 argument count' 1 l2 'proper usage with arguments:'
expect_failure 'heat argument count' 1 heat 'proper usage with arguments:'
expect_failure 'Eriksson argument count' 1 eriksson 'proper usage with arguments:'
expect_failure 'pure diffusion argument count' 1 pure_diffusion_igrm 'proper usage with arguments:'
expect_failure 'oil argument count' 1 oil 'proper usage with arguments:'
expect_failure 'iGRM L2 argument count' 1 igrm_l2 'proper usage with arguments:'
expect_failure 'iGRM heat argument count' 1 igrm_heat 'proper usage with arguments:'

# A command-line argument must contain exactly one integer token.
expect_failure 'L2 partial integer' 1 l2 'invalid integer argument:' \
    '2 junk' 2 2 1 1 1 1
expect_failure 'heat partial integer' 1 heat 'invalid integer argument:' \
    '2 junk' 1 0 0.1 1 1 1
expect_failure 'Eriksson partial integer' 1 eriksson 'invalid integer argument:' \
    '2 junk' 1 0 0.1 1 1 1
expect_failure 'pure diffusion partial integer' 1 pure_diffusion_igrm 'invalid integer argument:' \
    '3 junk' 1 1 1 1 0 0.1
expect_failure 'oil partial integer' 1 oil 'invalid integer argument:' \
    '2 junk' 1 1 1 1 0 0.1 0 0
expect_failure 'iGRM L2 partial integer' 1 igrm_l2 'invalid integer argument:' \
    '3 junk' 3 3 2 2 2 1 1 1 1 1 1
expect_failure 'iGRM heat partial integer' 1 igrm_heat 'invalid integer argument:' \
    '3 junk' 3 3 2 2 2 1 1 1 1 1 1 0 0.1

# The same exact-token rule applies to all real-valued driver arguments.
expect_failure 'heat partial real' 1 heat 'invalid real argument:' \
    2 1 0 '0.1 junk' 1 1 1
expect_failure 'Eriksson partial real' 1 eriksson 'invalid real argument:' \
    2 1 0 '0.1 junk' 1 1 1
expect_failure 'pure diffusion partial real' 1 pure_diffusion_igrm 'invalid real argument:' \
    3 1 1 1 1 0 '0.1 junk'
expect_failure 'oil partial real' 1 oil 'invalid real argument:' \
    2 1 1 1 1 0 '0.1 junk' 0 0
expect_failure 'iGRM L2 partial tau' 1 igrm_l2 'invalid real argument:' \
    3 3 3 2 2 2 1 1 1 1 1 1 '0.1 junk' dg
expect_failure 'iGRM heat partial real' 1 igrm_heat 'invalid real argument:' \
    3 3 3 2 2 2 1 1 1 1 1 1 0 '0.1 junk'

# Transient drivers reject a negative iteration count before MPI setup.
expect_failure 'heat negative steps' 1 heat 'number of time steps must be non-negative' \
    2 1 -1 0.1 1 1 1
expect_failure 'Eriksson negative steps' 1 eriksson 'number of time steps must be non-negative' \
    2 1 -1 0.1 1 1 1
expect_failure 'pure diffusion negative steps' 1 pure_diffusion_igrm \
    'number of time steps must be non-negative' 3 1 1 1 1 -1 0.1
expect_failure 'oil negative steps' 1 oil 'number of time steps must be non-negative' \
    2 1 1 1 1 -1 0.1 0 0
expect_failure 'iGRM heat negative steps' 1 igrm_heat \
    'number of time steps must be non-negative' \
    3 3 3 2 2 2 1 1 1 1 1 1 -1 0.1

# Basic discretization and process-grid constraints are checked before MPI.
expect_failure 'L2 zero elements' 1 l2 'number of elements must be positive' \
    0 2 2 1 1 1 1
expect_failure 'L2 negative order' 1 l2 'polynomial order must be non-negative' \
    2 2 2 -1 1 1 1
expect_failure 'L2 zero process-grid dimension' 1 l2 \
    'process-grid dimension must be positive' 2 2 2 1 0 1 1

expect_failure 'heat zero elements' 1 heat 'number of elements must be positive' \
    0 1 0 0.1 1 1 1
expect_failure 'heat negative order' 1 heat 'polynomial order must be non-negative' \
    2 -1 0 0.1 1 1 1
expect_failure 'heat zero time step' 1 heat 'time step must be positive' \
    2 1 0 0 1 1 1
expect_failure 'heat zero process-grid dimension' 1 heat \
    'process-grid dimension must be positive' 2 1 0 0.1 1 0 1
expect_failure 'Eriksson negative time step' 1 eriksson 'time step must be positive' \
    2 1 0 -0.1 1 1 1

expect_failure 'pure diffusion zero elements' 1 pure_diffusion_igrm \
    'number of elements must be positive' 0 1 1 1 1 0 0.1
expect_failure 'pure diffusion negative order' 1 pure_diffusion_igrm \
    'polynomial order must be non-negative' 3 -1 1 1 1 0 0.1
expect_failure 'pure diffusion zero time step' 1 pure_diffusion_igrm \
    'time step must be positive' 3 1 1 1 1 0 0
expect_failure 'pure diffusion zero process-grid dimension' 1 pure_diffusion_igrm \
    'process-grid dimension must be positive' 3 1 1 0 1 0 0.1

expect_failure 'oil zero elements' 1 oil 'number of elements must be positive' \
    0 1 1 1 1 0 0.1 0 0
expect_failure 'oil negative order' 1 oil 'polynomial order must be non-negative' \
    2 -1 1 1 1 0 0.1 0 0
expect_failure 'oil zero time step' 1 oil 'time step must be positive' \
    2 1 1 1 1 0 0 0 0
expect_failure 'oil zero process-grid dimension' 1 oil \
    'process-grid dimension must be positive' 2 1 1 1 0 0 0.1 0 0

# Derived spline extents must be checked before default-integer arithmetic.
expect_failure 'L2 overflowing spline extent' 1 l2 \
    'spline dimensions exceed supported integer range' \
    2147483647 2 2 1 1 1 1
expect_failure 'heat overflowing spline extent' 1 heat \
    'spline dimensions exceed supported integer range' \
    2147483647 1 0 0.1 1 1 1
expect_failure 'Eriksson overflowing spline extent' 1 eriksson \
    'spline dimensions exceed supported integer range' \
    2147483647 1 0 0.1 1 1 1
expect_failure 'pure diffusion enriched degree overflow' 1 pure_diffusion_igrm \
    'spline dimensions exceed supported integer range' \
    2147483645 0 1 1 1 0 0.1
expect_failure 'oil overflowing spline extent' 1 oil \
    'spline dimensions exceed supported integer range' \
    2147483647 1 1 1 1 0 0.1 0 0

# iGRM keeps the deliberate strict contract p_test > p_trial in each axis.
expect_failure 'iGRM L2 zero elements' 1 igrm_l2 'number of elements must be positive' \
    0 3 3 2 2 2 1 1 1 1 1 1
expect_failure 'iGRM L2 negative test degree' 1 igrm_l2 \
    'test polynomial degree must be non-negative' \
    3 3 3 -1 2 2 1 1 1 1 1 1
expect_failure 'iGRM L2 negative trial degree' 1 igrm_l2 \
    'trial polynomial degree must be non-negative' \
    3 3 3 2 2 2 -1 1 1 1 1 1
expect_failure 'iGRM L2 unenriched test degree' 1 igrm_l2 \
    'test polynomial degree must be greater than trial polynomial degree' \
    3 3 3 1 2 2 1 1 1 1 1 1
expect_failure 'iGRM L2 zero process-grid dimension' 1 igrm_l2 \
    'process-grid dimension must be positive' \
    3 3 3 2 2 2 1 1 1 1 0 1
expect_failure 'iGRM L2 zero tau' 1 igrm_l2 'time-scheme tau must be positive' \
    3 3 3 2 2 2 1 1 1 1 1 1 0 dg

expect_failure 'iGRM heat zero elements' 1 igrm_heat 'number of elements must be positive' \
    0 3 3 2 2 2 1 1 1 1 1 1 0 0.1
expect_failure 'iGRM heat negative test degree' 1 igrm_heat \
    'test polynomial degree must be non-negative' \
    3 3 3 -1 2 2 1 1 1 1 1 1 0 0.1
expect_failure 'iGRM heat negative trial degree' 1 igrm_heat \
    'trial polynomial degree must be non-negative' \
    3 3 3 2 2 2 -1 1 1 1 1 1 0 0.1
expect_failure 'iGRM heat unenriched test degree' 1 igrm_heat \
    'test polynomial degree must be greater than trial polynomial degree' \
    3 3 3 1 2 2 1 1 1 1 1 1 0 0.1
expect_failure 'iGRM heat zero process-grid dimension' 1 igrm_heat \
    'process-grid dimension must be positive' \
    3 3 3 2 2 2 1 1 1 1 0 1 0 0.1
expect_failure 'iGRM heat zero time step' 1 igrm_heat 'time step must be positive' \
    3 3 3 2 2 2 1 1 1 1 1 1 0 0

expect_failure 'iGRM L2 overflowing spline extent' 1 igrm_l2 \
    'spline dimensions exceed supported integer range' \
    2147483647 3 3 2 2 2 1 1 1 1 1 1
expect_failure 'iGRM heat overflowing spline extent' 1 igrm_heat \
    'spline dimensions exceed supported integer range' \
    2147483647 3 3 2 2 2 1 1 1 1 1 1 0 0.1

# Space validation rejects degrees requiring unavailable Gauss rules.
expect_failure 'L2 unsupported degree' 1 l2 \
    'unsupported polynomial degree in ADS space' \
    1 1 1 10 1 1 1
expect_failure 'heat unsupported degree' 1 heat \
    'unsupported polynomial degree in ADS space' \
    1 10 0 0.1 1 1 1
expect_failure 'Eriksson unsupported degree' 1 eriksson \
    'unsupported polynomial degree in ADS space' \
    1 10 0 0.1 1 1 1
expect_failure 'oil unsupported degree' 1 oil \
    'unsupported polynomial degree in ADS space' \
    1 10 1 1 1 0 0.1 0 0
expect_failure 'pure diffusion unsupported enriched degree' 1 pure_diffusion_igrm \
    'unsupported polynomial degree in iGRM test space' \
    1 9 1 1 1 0 0.1 dg
expect_failure 'iGRM L2 unsupported test degree' 1 igrm_l2 \
    'unsupported polynomial degree in iGRM test space' \
    1 1 1 10 10 10 9 9 9 1 1 1 dg
expect_failure 'iGRM L2 unsupported y degree' 1 igrm_l2 \
    'iGRM test space along y: 10' \
    1 1 1 2 10 2 1 9 1 1 1 1 dg
expect_failure 'iGRM L2 unsupported z degree' 1 igrm_l2 \
    'iGRM test space along z: 10' \
    1 1 1 2 2 10 1 1 9 1 1 1 dg
expect_failure 'iGRM heat unsupported test degree' 1 igrm_heat \
    'unsupported polynomial degree in iGRM test space' \
    1 1 1 10 10 10 9 9 9 1 1 1 0 0.1 dg

# iGRM drivers accept DG/PR/BE only; Forward Euler has its own diagnostic.
expect_failure 'pure diffusion unknown scheme' 1 pure_diffusion_igrm 'unknown time scheme:' \
    3 1 1 1 1 0 0.1 cn
expect_failure 'pure diffusion Forward Euler' 1 pure_diffusion_igrm \
    'forward euler is not an iGRM time-scheme option' 3 1 1 1 1 0 0.1 fe
expect_failure 'iGRM L2 unknown scheme' 1 igrm_l2 'unknown time scheme:' \
    3 3 3 2 2 2 1 1 1 1 1 1 cn
expect_failure 'iGRM L2 Forward Euler' 1 igrm_l2 \
    'forward euler is not an iGRM time-scheme option' \
    3 3 3 2 2 2 1 1 1 1 1 1 fe
expect_failure 'iGRM heat unknown scheme' 1 igrm_heat 'unknown time scheme:' \
    3 3 3 2 2 2 1 1 1 1 1 1 0 0.1 cn
expect_failure 'iGRM heat Forward Euler' 1 igrm_heat \
    'forward euler is not an iGRM time-scheme option' \
    3 3 3 2 2 2 1 1 1 1 1 1 0 0.1 fe

# Oil has a variable-length suffix describing pumps and drains.
expect_failure 'oil missing pump/drain suffix' 1 oil 'proper usage with arguments:' \
    2 1 1 1 1 0 0.1
expect_failure 'oil negative pump count' 1 oil 'number of pumps must be non-negative' \
    2 1 1 1 1 0 0.1 -1 0
expect_failure 'oil missing pump coordinates' 1 oil 'proper usage with arguments:' \
    2 1 1 1 1 0 0.1 1 0
expect_failure 'oil invalid pump coordinate' 1 oil 'invalid real argument:' \
    2 1 1 1 1 0 0.1 1 bad 0 0 0
expect_failure 'oil negative drain count' 1 oil 'number of drains must be non-negative' \
    2 1 1 1 1 0 0.1 0 -1
expect_failure 'oil missing drain coordinates' 1 oil 'proper usage with arguments:' \
    2 1 1 1 1 0 0.1 0 1 0
expect_failure 'oil invalid drain coordinate' 1 oil 'invalid real argument:' \
    2 1 1 1 1 0 0.1 0 1 bad 0 0
expect_failure 'oil extra suffix argument' 1 oil 'proper usage with arguments:' \
    2 1 1 1 1 0 0.1 0 0 extra
expect_failure 'oil huge pump count is rejected before allocation' 1 oil \
    'proper usage with arguments:' 2 1 1 1 1 0 0.1 2147483647 0

# This failure occurs after MPI initialization and used to return success
# because initialize_setup used an unnumbered STOP.
if [[ ${SKIP_MPI_CASES:-0} != 1 ]]; then
    expect_failure 'elements smaller than process grid' 2 l2 \
        'Number of elements smaller than number of processors' \
        1 1 1 1 2 1 1
fi

if [[ $failures -eq 0 ]]; then
    printf 'OK (%d negative driver CLI checks)\n' "$checks"
    exit 0
fi

printf 'FAILED (%d of %d negative driver CLI checks)\n' "$failures" "$checks"
exit 1
