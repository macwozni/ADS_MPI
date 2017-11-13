cmake_minimum_required(VERSION 3.1)
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

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-pthread -funroll-all-loops -fno-f2c -O3 -std=f2003")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -std=f03")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -std03")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -O0 -g -std03")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2  -std03 -std=f03")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -std03 -std=03")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

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
