#!/bin/bash

export MPIF90=/opt/mpich_gnu/bin/mpif90 
export F90=gfortran
export PFUNIT=/opt/pfunit/pfunit-parallel
# make clean
make all
