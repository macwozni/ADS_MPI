EXEC_NAME := igrm_stokes
SOURCES := input_data.F90 RHS_fun.F90 igrm_stokes_solver.F90 main.F90
PROBLEM_MODULES := input_data rhs_fun igrm_stokes_solver
DEFAULT_ARGS := 4 4 4 2 2 2 2 2 2 1 1 1
USAGE := <nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>

include ../problem.mk
