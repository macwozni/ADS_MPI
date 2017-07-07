

program main

use parallelism, ONLY : MYRANK
use parallelism, ONLY : PRINTRANK,InitializeParallelism
use communicators, ONLY : CreateCommunicators
use RHS_eq
use ADSS
use input_data

implicit none

include "mpif.h"

! Iteration counter
integer :: iter = 0

integer(kind=4) :: ierr

type   (ADS_setup) :: ads

#ifdef DEBUG
   write (*,*) 'debug'
#endif
   
t = 0

! -------------------------------------------------------------------
! Code
! -------------------------------------------------------------------

  call InitializeParameters

  ! prepare the problem dimensions

  call InitializeParallelism
  call CreateCommunicators
  call Initialize(SIZE,SIZE,SIZE,ORDER,ORDER,ORDER,ads)

  allocate(Kqvals(ads%p(1)+1,ads%p(2)+1,ads%p(3)+1,ads%maxex-ads%minex+1,ads%maxey-ads%miney+1,ads%maxez-ads%minez+1))
  call InitInputData
  call PrecomputeKq(ads)

  ! Iterations
  do iter = 0,steps

    l2norm=0
    call Step(iter,ComputePointForRHS,ads)
    call MPI_Reduce(l2norm, fullnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (MYRANK == 0) then
      write(*,*)iter, 'L2 norm:', fullnorm
    endif
    t = t + Dt
    
  enddo

  call ComputeResults
  call Cleanup(ads)

end

