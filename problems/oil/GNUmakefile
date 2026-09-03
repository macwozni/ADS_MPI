EXEC_NAME := oil
SOURCES := input_data.F90 RHS_fun.F90 main.F90
PROBLEM_MODULES := input_data rhs_fun
DEFAULT_ARGS := 2 1 1 1 1 1 0.1 1 0.5 0.5 0.5 1 0.25 0.25 0.25
USAGE := <size> <order> <procx> <procy> <procz> <steps> <dt> <npumps> (<pump_x> <pump_y> <pump_z>)... <ndrains> (<drain_x> <drain_y> <drain_z>)...
PROBLEM_RUN_ENV = $(if $(strip $(OIL_SEED)),ADS_OIL_RANDOM_SEED='$(OIL_SEED)')

include ../problem.mk
