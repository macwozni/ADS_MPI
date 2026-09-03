EXEC_NAME := igrm_heat
SOURCES := input_data.F90 RHS_fun.F90 main.F90
PROBLEM_MODULES := input_data rhs_fun
DEFAULT_ARGS := 2 2 2 3 3 3 2 2 2 1 1 1 1 0.001 dg
USAGE := <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <steps> <dt> [dg|pr|be]

include ../problem.mk
