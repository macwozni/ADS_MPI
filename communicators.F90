module communicators

use parallelism

implicit none

! Total number of processors along X, Y and Z
integer, parameter :: NRPROCXMAX = 128
integer, parameter :: NRPROCYMAX = 128
integer, parameter :: NRPROCZMAX = 128

! Global ranks of processors in the cube
integer :: processors(NRPROCXMAX,NRPROCYMAX,NRPROCZMAX)

! Group involving processors along X, Y, Z
integer :: GROUPX(NRPROCYMAX,NRPROCZMAX)
integer :: GROUPY(NRPROCXMAX,NRPROCZMAX)
integer :: GROUPZ(NRPROCXMAX,NRPROCYMAX)

! Communicatorx along X, Y, Z
integer :: COMMXALL(NRPROCYMAX,NRPROCZMAX)
integer :: COMMYALL(NRPROCXMAX,NRPROCZMAX)
integer :: COMMZALL(NRPROCXMAX,NRPROCYMAX)

! Local communicators
integer :: COMMX,COMMY,COMMZ

! Corresponding contexts for SCALAPACK calls
integer :: CONTEXTX
integer :: CONTEXTY
integer :: CONTEXTZ

contains


! -------------------------------------------------------------------
! Creates process groups and communicators for each 'fibre' in each
! direction.
! -------------------------------------------------------------------
subroutine CreateCommunicators
include "mpif.h"
integer :: group_comm_world
integer :: comm_myrank_local
integer :: processors_X(NRPROCX)
integer :: processors_Y(NRPROCY)
integer :: processors_Z(NRPROCZ)
integer :: i,j,k
integer :: irank
integer :: iprint
integer :: ierr

  iprint = 0
  if (iprint == 1) then
    write(*,*)MYRANK,'NRPROC',NRPROC
  endif

  call mpi_comm_group(MPI_COMM_WORLD,group_comm_world,ierr)

  if (ierr /= 0) then
    write(*,*)MYRANK,': main: error calling mpi_comm_group!'
    call abort
  endif
  if (iprint == 1) then
    write(*,*)MYRANK,'got group',group_comm_world
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  do i = 1,NRPROCX
    do j = 1,NRPROCY
      do k = 1,NRPROCZ
        processors(i,j,k) = (i-1)+(j-1)*NRPROCX+(k-1)*NRPROCX*NRPROCY
      enddo
    enddo
  enddo

  do i = 1,NRPROCX
    do j = 1,NRPROCY
      processors_Z(1:NRPROCZ) = processors(i,j,1:NRPROCZ)
      call mpi_group_incl(group_comm_world,NRPROCZ,processors_Z,GROUPZ(i,j),ierr)
      if (ierr /=0 ) then
        write(*,*)MYRANK,': main: error calling mpi_group_incl for Z',i,j
        call abort
      endif
    enddo
  enddo
  do i = 1,NRPROCX
    do k = 1,NRPROCZ
      processors_Y(1:NRPROCY) = processors(i,1:NRPROCY,k)
      call mpi_group_incl(group_comm_world,NRPROCY,processors_Y,GROUPY(i,k),ierr)
      if (ierr /= 0) then
        write(*,*)MYRANK,': main: error calling mpi_group_incl for Y',i,k
        call abort
      endif
    enddo
  enddo
  do j = 1,NRPROCY
    do k = 1,NRPROCZ
      processors_X(1:NRPROCX) = processors(1:NRPROCX,j,k)
      call mpi_group_incl(group_comm_world,NRPROCX,processors_X,GROUPX(j,k),ierr)
      if (ierr /= 0) then
        write(*,*)MYRANK,': main: error calling mpi_group_incl for X',j,k
        call abort
      endif
    enddo
  enddo

  if (iprint == 1) then
    call PrintGroups
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  ! create the new communicators
  do i = 1,NRPROCX
    do j = 1,NRPROCY
      call mpi_comm_create(MPI_COMM_WORLD,GROUPZ(i,j),comm_myrank_local,ierr)
      COMMZALL(i,j) = comm_myrank_local
      if (ierr /= 0) then
        write(*,*)MYRANK,': main: error calling mpi_com_create for Z',i,j
        call abort
      endif
    enddo
  enddo
  do i = 1,NRPROCX
    do k = 1,NRPROCZ
      call mpi_comm_create(MPI_COMM_WORLD,GROUPY(i,k),comm_myrank_local,ierr)
      COMMYALL(i,k) = comm_myrank_local
      if (ierr /= 0) then
        write(*,*)MYRANK,': main: error calling mpi_com_create for Y',i,k
        call abort
      endif
    enddo
  enddo
  do j = 1,NRPROCY
    do k = 1,NRPROCZ
      call mpi_comm_create(MPI_COMM_WORLD,GROUPX(j,k),comm_myrank_local,ierr)
      COMMXALL(j,k) = comm_myrank_local
      if (ierr /= 0) then
        write(*,*)MYRANK,': main: error calling mpi_com_create for X',j,k
        call abort
      endif
    enddo
  enddo
  if (iprint == 1) then
    call PrintCommunicators
  endif      

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  ! extract local communicators
  COMMX = COMMXALL(myranky+1,myrankz+1)
  COMMY = COMMYALL(myrankx+1,myrankz+1)
  COMMZ = COMMZALL(myrankx+1,myranky+1)

  if (iprint == 1) then
    write(*,*)PRINTRANK,'COMMX(Y,Z)',COMMX,COMMY,COMMZ
  endif

end subroutine


! -------------------------------------------------------------------
! Prints groups used later to create communicators
! -------------------------------------------------------------------
subroutine PrintGroups
integer :: i, j, k

  do i = 1,NRPROCX
    do j = 1,NRPROCY
      write(*,*)PRINTRANK,'GROUPZ(',i,j,')',GROUPZ(i,j)
    enddo     
  enddo
  do i = 1,NRPROCX
    do k = 1,NRPROCZ
      write(*,*)PRINTRANK,'GROUPY(',i,k,')',GROUPY(i,k)
    enddo     
  enddo
  do j = 1,NRPROCY
    do k = 1,NRPROCZ
      write(*,*)PRINTRANK,'GROUPX(',j,k,')',GROUPX(j,k)
    enddo     
  enddo

end subroutine


! -------------------------------------------------------------------
! Prints communicators
! -------------------------------------------------------------------
subroutine PrintCommunicators
integer :: i, j, k

  do i = 1,NRPROCX
    do j = 1,NRPROCY
      write(*,*)PRINTRANK,'COMMZALL(',i,j,')',COMMZALL(i,j)
    enddo     
  enddo
  do i = 1,NRPROCX
    do k = 1,NRPROCZ
      write(*,*)PRINTRANK,'COMMYALL(',i,k,')',COMMYALL(i,k)
    enddo     
  enddo
  do j = 1,NRPROCY
    do k = 1,NRPROCZ
      write(*,*)PRINTRANK,'COMMXALL(',j,k,')',COMMXALL(j,k)
    enddo     
  enddo

end subroutine

end module
