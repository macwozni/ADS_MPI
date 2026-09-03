EXEC_NAME := igrm_pollution
SOURCES := input_data.F90 RHS_fun.F90 pollution_dpg_solver.F90 main.F90
PROBLEM_MODULES := input_data rhs_fun pollution_dpg_solver
DEFAULT_ARGS := 4 0 1 0 2 1 1 1 1 1
USAGE := <N> <adapt:0|1> <p_trial> <C_trial> <p_test> <C_test> <steps> <procx> <procy> <procz>

include ../problem.mk
