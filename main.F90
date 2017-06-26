

program main

use parallelism, ONLY : MYRANK
use debug, ONLY : iinfo,idebug,iprint
use stuff
use RHS_eq
use ADS
use input_data

implicit none

include "mpif.h"

! Iteration counter
integer :: iter = 0

integer(kind=4) :: ierr

t = 0

! -------------------------------------------------------------------
! Code
! -------------------------------------------------------------------

  call InitializeParameters

  ! prepare the problem dimensions
  px = ORDER ! order
  py = ORDER ! order
  pz = ORDER ! order
  nx = SIZE  ! intervals
  ny = SIZE  ! intervals
  nz = SIZE  ! intervals

  call Initialize

  allocate(Kqvals(px+1,py+1,pz+1,maxex-minex+1,maxey-miney+1,maxez-minez+1))
  call InitInputData
  call PrecomputeKq

  ! Iterations
  do iter = 0,steps

    l2norm=0
    call Step(iter,ComputePointForRHS)
    call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (MYRANK == 0) then
      write(*,*)iter, 'L2 norm:', fullnorm
    endif
    t = t + Dt
    
  enddo

  call ComputeResults
  call Cleanup

end

