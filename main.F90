

program main

use parallelism, ONLY : NRPROCX,NRPROCY,NRPROCZ
use debug, ONLY : iinfo,idebug,iprint
use input_data
use stuff
use ADS

implicit none

! Iteration counter
integer :: iter = 0

t = 0

! -------------------------------------------------------------------
! Code
! -------------------------------------------------------------------

  call InitializeParameters
  call Initialize

  ! prepare the problem dimensions
  px = ORDER ! order
  py = ORDER ! order
  pz = ORDER ! order
  nx = SIZE  ! intervals
  ny = SIZE  ! intervals
  nz = SIZE  ! intervals

  if (iinfo == 1) then
    write(*,*)'px',px,'py',py,'pz',pz, &
    'nx',nx,'ny',ny,'nz',nz, &
    'size of Ux',nx+px+2,'size of Uy',ny+py+2,'size of Uz',nz+pz+2
  endif

  if (SIZE<NRPROCX .or. SIZE<NRPROCY .or. SIZE<NRPROCZ) then
    write(*,*)'Number of elements smaller than number of processors'
    stop
  endif

  KLx = px
  KUx = px
  KLy = py
  KUy = py
  KLz = pz
  KUz = pz

  call ComputeDecomposition

  if (idebug == 1) then
    call ValidateDimensions
  endif

  if (iprint == 1) then
    call PrintDecompositionInfo
  endif

  call AllocateArrays
  call PrepareKnot(Ux,nx,px)
  call PrepareKnot(Uy,ny,py)
  call PrepareKnot(Uz,nz,pz)
  allocate(Kqvals(px+1,py+1,pz+1,maxex-minex+1,maxey-miney+1,maxez-minez+1))
  call InitInputData
  call PrecomputeKq

  ! Iterations
  do iter = 0,steps

    ! if (MYRANK == 0) then
    !   write(*,*)'Iteration',iter,'/',steps
    !   write(*,*)'t = ',t
    ! endif

    call Step(iter,ComputePointForRHS)
    t = t + Dt

    ! if (mod(iter, 100) == 0) then
    !   call PrintSolution(iter, t)
    ! endif
  enddo

  call ComputeResults
  call Cleanup

end

