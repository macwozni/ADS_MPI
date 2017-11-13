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
  set (CMAKE_Fortran_FLAGS_RELEASE "-pthread -funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested) 
  ######################3
  # NOT WORKING #
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -module . -std95")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

include_directories("${PROJECT_SOURCE_DIR}/src")

set(SRC src)
set(ads_SRC
  ${SRC}/knot_vector.F95
  ${SRC}/parallelism.F95
  ${SRC}/communicators.F95
  ${SRC}/gauss.F95
  ${SRC}/basis.F95
  ${SRC}/math.F95
  ${SRC}/reorderRHS.F95
  ${SRC}/my_mpi.F95
  ${SRC}/projection_engine.F95
  ${SRC}/analysis.F95
  ${SRC}/int2str.F95
  ${SRC}/plot.F95
  ${SRC}/gnuplot.F95
  ${SRC}/vtk.F95
  ${SRC}/ADS.F95
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
    problems/oil/input_data.F95
    problems/oil/RHS_eq.F95
    problems/oil/main.F95)

endif()
