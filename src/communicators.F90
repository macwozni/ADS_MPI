!------------------------------------------------------------------------------
! AGH University of Science and Technology, Krakow, Poland
!------------------------------------------------------------------------------
!
! MODULE: communicators
!
!> @author
!> Maciej Wozniak
!
! DESCRIPTION:
!> This module contains communicators logic.
!
! REVISION HISTORY:
! 21 11 2017 - Initial Version
! 
!------------------------------------------------------------------------------

module communicators

implicit none

!> Total number of processors along X
integer(kind = 4), parameter :: NRPROCXMAX = 128
!> Total number of processors along Y
integer(kind = 4), parameter :: NRPROCYMAX = 128
!> Total number of processors along Z
integer(kind = 4), parameter :: NRPROCZMAX = 128

!> Global ranks of processors in the cube
integer(kind = 4) :: processors(NRPROCXMAX, NRPROCYMAX, NRPROCZMAX)

!> Group involving processors along X
integer(kind = 4) :: GROUPX(NRPROCYMAX, NRPROCZMAX)
!> Group involving processors along Y
integer(kind = 4) :: GROUPY(NRPROCXMAX, NRPROCZMAX)
!> Group involving processors along Z
integer(kind = 4) :: GROUPZ(NRPROCXMAX, NRPROCYMAX)

!> Communicatorx along Z
integer(kind = 4) :: COMMXALL(NRPROCYMAX, NRPROCZMAX)
!> Communicatorx along Y
integer(kind = 4) :: COMMYALL(NRPROCXMAX, NRPROCZMAX)
!> Communicatorx along Z
integer(kind = 4) :: COMMZALL(NRPROCXMAX, NRPROCYMAX)

!> Local communicators X
integer(kind = 4) :: COMMX
!> Local communicators Y
integer(kind = 4) :: COMMY
!> Local communicators Z
integer(kind = 4) :: COMMZ

!> Corresponding contexts for SCALAPACK calls
integer(kind = 4) :: CONTEXTX
!> Corresponding contexts for SCALAPACK calls
integer(kind = 4) :: CONTEXTY
!> Corresponding contexts for SCALAPACK calls
integer(kind = 4) :: CONTEXTZ

PRIVATE :: GROUPX, GROUPY, GROUPZ
PRIVATE :: COMMXALL, COMMYALL, COMMZALL
PRIVATE :: CONTEXTX, CONTEXTY, CONTEXTZ
!PRIVATE :: PrintGroups, PrintCommunicators
PROTECTED :: processors, COMMX, COMMY, COMMZ

contains


!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Creates process groups and communicators for each 'fibre' in each direction.
!> A subroutine initializes some variables in module.
!
! Output:
! -------
!> @param[out] mierr - error code
! -------------------------------------------------------------------
subroutine CreateCommunicators(mierr)
   use parallelism, ONLY: MYRANK, NRPROC, MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
   NRPROCX, NRPROCY, NRPROCZ
   USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   use mpi
   implicit none
   integer(kind = 4), intent(out) :: mierr
   integer(kind = 4) :: group_comm_world
   integer(kind = 4) :: comm_myrank_local
   integer(kind = 4) :: processors_X(NRPROCX)
   integer(kind = 4) :: processors_Y(NRPROCY)
   integer(kind = 4) :: processors_Z(NRPROCZ)
   integer(kind = 4) :: i, j, k
   integer(kind = 4) :: irank
   integer(kind = 4) :: ierr

#ifdef IPRINT
   write(*, *) MYRANK, 'NRPROC', NRPROC
#endif

call mpi_comm_group(MPI_COMM_WORLD, group_comm_world, ierr)

   if (ierr /= 0) then
      write(*, *) MYRANK, ': main: error calling mpi_comm_group!'
      STOP 4
   endif
#ifdef IPRINT
   write(*, *) MYRANK, 'got group', group_comm_world
#endif

call mpi_barrier(MPI_COMM_WORLD, ierr)

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         do k = 1, NRPROCZ
            processors(i, j, k) = (i - 1)+(j - 1) * NRPROCX + (k - 1) * NRPROCX * NRPROCY
         enddo
      enddo
   enddo

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         processors_Z(1:NRPROCZ) = processors(i, j, 1:NRPROCZ)
         call mpi_group_incl(group_comm_world, NRPROCZ, processors_Z, GROUPZ(i, j), ierr)
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for Z', i, j
            STOP 4
         endif
      enddo
   enddo
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         processors_Y(1:NRPROCY) = processors(i, 1:NRPROCY, k)
         call mpi_group_incl(group_comm_world, NRPROCY, processors_Y, GROUPY(i, k), ierr)
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for Y', i, k
            STOP 4
         endif
      enddo
   enddo
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         processors_X(1:NRPROCX) = processors(1:NRPROCX, j, k)
         call mpi_group_incl(group_comm_world, NRPROCX, processors_X, GROUPX(j, k), ierr)
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for X', j, k
            STOP 4
         endif
      enddo
   enddo

#ifdef IPRINT
   call PrintGroups
#endif

call mpi_barrier(MPI_COMM_WORLD, ierr)

   ! create the new communicators
   do i = 1, NRPROCX
      do j = 1, NRPROCY
         call mpi_comm_create(MPI_COMM_WORLD, GROUPZ(i, j), comm_myrank_local, ierr)
         COMMZALL(i, j) = comm_myrank_local
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for Z', i, j  
            STOP 4
         endif
      enddo
   enddo
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         call mpi_comm_create(MPI_COMM_WORLD, GROUPY(i, k), comm_myrank_local, ierr)
         COMMYALL(i, k) = comm_myrank_local
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for Y', i, k
            STOP 4
         endif
      enddo
   enddo
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         call mpi_comm_create(MPI_COMM_WORLD, GROUPX(j, k), comm_myrank_local, ierr)
         COMMXALL(j, k) = comm_myrank_local
         if (ierr /= 0) then
            write(ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for X', j, k   
            STOP 4
         endif
      enddo
   enddo
#ifdef IPRINT
   call PrintCommunicators
#endif     

call mpi_barrier(MPI_COMM_WORLD, ierr)

   ! extract local communicators
   COMMX = COMMXALL(myranky + 1, myrankz + 1)
   COMMY = COMMYALL(myrankx + 1, myrankz + 1)
   COMMZ = COMMZALL(myrankx + 1, myranky + 1)

#ifdef IPRINT
   write(*, *) PRINTRANK, 'COMMX(Y,Z)', COMMX, COMMY, COMMZ
#endif

mierr = 0

end subroutine CreateCommunicators

!!!!! przeniesc do debug
!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Prints groups used later to create communicators.
!> For debug only
! -------------------------------------------------------------------
subroutine PrintGroups
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK
   implicit none
   integer(kind = 4) :: i, j, k

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         write(*, *) PRINTRANK, 'GROUPZ(', i, j, ')', GROUPZ(i, j)
      enddo
   enddo
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         write(*, *) PRINTRANK, 'GROUPY(', i, k, ')', GROUPY(i, k)
      enddo
   enddo
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         write(*, *) PRINTRANK, 'GROUPX(', j, k, ')', GROUPX(j, k)
      enddo
   enddo

end subroutine PrintGroups

!!! przeniesc do debug
!---------------------------------------------------------------------------  
!> @author Maciej Wozniak
!>
!> @brief
!> Prints communicators.
!> For debug only
! -------------------------------------------------------------------
subroutine PrintCommunicators
   use parallelism, ONLY: PRINTRANK, NRPROCX, NRPROCY, NRPROCZ
   implicit none
   integer(kind = 4) :: i, j, k

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         write(*, *) PRINTRANK, 'COMMZALL(', i, j, ')', COMMZALL(i, j)
      enddo
   enddo
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         write(*, *) PRINTRANK, 'COMMYALL(', i, k, ')', COMMYALL(i, k)
      enddo
   enddo
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         write(*, *) PRINTRANK, 'COMMXALL(', j, k, ')', COMMXALL(j, k)
      enddo
   enddo

end subroutine PrintCommunicators

!!!!! dodac czyszczenie komunikatorow


end module communicators
