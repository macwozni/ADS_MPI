EXEC_NAME := pure_diffusion_igrm
SOURCES := input_data.F90 RHS_fun.F90 main.F90
PROBLEM_MODULES := input_data rhs_fun
DEFAULT_ARGS := 3 1 1 1 1 1 0.1 dg
USAGE := <size> <order> <procx> <procy> <procz> <steps> <dt> [dg|pr|be]

include ../problem.mk
