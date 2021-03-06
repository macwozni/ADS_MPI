cmake_minimum_required(VERSION 3.5)
project (ADS)
enable_language (Fortran)


# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
#  set(CMAKE_BUILD_TYPE "Debug")
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
    set(dialect "-ffree-form -std=f2008 -fimplicit-none -cpp -Wall -pedantic")
    set(bounds "-fbounds-check")
    set(realeaseopts "-O3  -pthread -funroll-all-loops")
    set(debugopts "-O0 -g3 -fbacktrace -fcheck=all -fdump-fortran-optimized")
#    -DIINFO -DIDEBUG -DIPRINT")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-std08  -module . -implicitnone -fpp -Wall -pedantic")
    set(bounds " -check bounds")
    set(realeaseopts "-funroll-all-loops -O3")
    set(debugopts "-O0 -g -traceback -check all -debug extended -debug-parameters all")
#    -DIINFO -DIDEBUG -DIPRINT")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(bounds " -check bounds")
    set(realeaseopts "-Munroll -Mipa -O3 -DPGI=1 -ta=host")
    set(debugopts "-O0 -g -traceback -check all -debug extended -debug-parameters all -DPGI=1")
#    -DIINFO -DIDEBUG -DIPRINT")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debugopts}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_RELEASE_FLAGS} ${realeaseopts}")

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DPERFORMANCE=1")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -DPERFORMANCE=1")
#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories("${PROJECT_SOURCE_DIR}/src")

set(SRC src)
set(ads_SRC
  ${SRC}/Setup.F90
  ${SRC}/knot_vector.F90
  ${SRC}/parallelism.F90
  ${SRC}/communicators.F90
  ${SRC}/gauss.F90
  ${SRC}/basis.F90
  ${SRC}/math.F90
  ${SRC}/reorderRHS.F90
  ${SRC}/my_mpi.F90
  ${SRC}/sparse.F90
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

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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

  define_problem(L2
    problems/L2Projection/input_data.F90
    problems/L2Projection/RHS_eq.F90
    problems/L2Projection/main.F90)

  define_problem(heat
    problems/heat/input_data.F90
    problems/heat/RHS_eq.F90
    problems/heat/main.F90)
    
endif()
