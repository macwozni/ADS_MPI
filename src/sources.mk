# Ordered manifest of the ADS core library.
#
# Keep this list in Fortran module dependency order.  src/GNUmakefile turns the
# order into prerequisites, so this file is the single source of truth for the
# library build and for source/test coverage checks.
SRC_FILES := \
	Setup.F90 \
	Interfaces.F90 \
	argument_parser.F90 \
	knot_vector.F90 \
	parallelism.F90 \
	communicators.F90 \
	gauss.F90 \
	basis.F90 \
	math.F90 \
	reorderRHS.F90 \
	my_mpi.F90 \
	sparse.F90 \
	RHS_eq.F90 \
	igrm_space.F90 \
	operator_assembly.F90 \
	rhs_assembly.F90 \
	solution_reconstruction.F90 \
	projection_engine.F90 \
	ads_lifecycle.F90 \
	mumps_solver.F90 \
	ads_directional_solve.F90 \
	utils.F90 \
	plot.F90 \
	gnuplot.F90 \
	vtk.F90 \
	solution_output.F90 \
	ADS.F90 \
	time_scheme.F90

# Module files produced by SRC_FILES.  These are listed explicitly so `clean`
# removes only files owned by src and leaves problem-local modules untouched.
SRC_MODULES := setup interfaces argument_parser knot_vector parallelism \
	communicators gauss basis math reorderrhs my_mpi sparse rhs_eq igrm_space \
	operator_assembly rhs_assembly solution_reconstruction projection_engine \
	ads_lifecycle mumps_solver ads_directional_solve utils plot gnuplot vtk \
	solution_output adss time_scheme
