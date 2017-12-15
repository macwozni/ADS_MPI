cmake_minimum_required(VERSION 3.5)
project (ADS)
enable_language (Fortran)


# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()


# make sure that the default is a RELEASE
#if (NOT CMAKE_BUILD_TYPE)
#  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
#      "Choose the type of build, options are: None Debug Release."
#      FORCE)
#endif (NOT CMAKE_BUILD_TYPE)


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none  -fno-f2c")
    set(bounds "-fbounds-check")
    set(realeaseopts "-O3  -pthread -funroll-all-loops")
    set(debugopts "-O0 -g3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-std08  -module . -implicitnone")
    set(bounds " -check bounds")
    set(realeaseopts "-funroll-all-loops -O3")
    set(debugopts "-O0 -g -traceback")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debugopts}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_RELEASE_FLAGS} ${realeaseopts}")


include_directories("${PROJECT_SOURCE_DIR}/src")

set(SRC src)
set(ads_SRC
  ${SRC}/knot_vector.F90
  ${SRC}/parallelism.F90
  ${SRC}/communicators.F90
  ${SRC}/gauss.F90
  ${SRC}/basis.F90
  ${SRC}/math.F90
  ${SRC}/reorderRHS.F90
  ${SRC}/my_mpi.F90
  ${SRC}/projection_engine.F90
  ${SRC}/analysis.F90
  ${SRC}/int2str.F90
  ${SRC}/plot.F90
  ${SRC}/gnuplot.F90
  ${SRC}/vtk.F90
  ${SRC}/ADS.F90
  )

# --------------------------------------------------------------------
# Libraries
# --------------------------------------------------------------------
set(ads_LIBS)

# Numerical libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

list(APPEND ads_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
list(APPEND ads_LIBS ${MPI_Fortran_LIBRARIES})

# Target library definition
add_library(ads STATIC ${ads_SRC})
target_link_libraries(ads ${ads_LIBS})


# --------------------------------------------------------------------
# Problems
# --------------------------------------------------------------------

option(SKIP_PROBLEMS "skip compiling example problems" OFF)

function(define_problem name)
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} ads)
endfunction()

if (NOT SKIP_PROBLEMS)

  define_problem(oil
    problems/oil/input_data.F90
    problems/oil/RHS_eq.F90
    problems/oil/main.F90)

endif()
