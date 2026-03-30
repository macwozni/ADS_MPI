!------------------------------------------------------------------------------
!
! MODULE: communicators
!
! DESCRIPTION:
!> @file communicators.F90
!> @brief Module creating and storing MPI subcommunicators aligned with
!> the logical process grid.
!>
!> @details
!> This module provides the communicator-management layer used by the
!> project to build one-dimensional MPI fibres along the three logical
!> coordinate directions of the process cube.
!>
!> The stored data include:
!> - the global-rank layout of the logical process cube,
!> - MPI groups associated with fibres parallel to the coordinate axes,
!> - communicator tables for all such fibres,
!> - the three communicators local to the current process.
!>
!> These communicators are subsequently consumed by higher-level
!> routines performing directional gather, scatter, and solve stages,
!> in particular within the ADS workflow and associated MPI wrappers.
!
!------------------------------------------------------------------------------
module communicators

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Maximum supported number of processes in the first
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4), parameter :: NRPROCXMAX = 128

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Maximum supported number of processes in the second
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4), parameter :: NRPROCYMAX = 128

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Maximum supported number of processes in the third
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4), parameter :: NRPROCZMAX = 128

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Table of global MPI ranks arranged as a logical process cube.
!>
!> @details
!> The entry `processors(i,j,k)` stores the linear rank associated with
!> process-grid coordinates \f$(i-1,j-1,k-1)\f$.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: processors(NRPROCXMAX, NRPROCYMAX, NRPROCZMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief MPI groups corresponding to fibres parallel to the first
!> process-grid direction.
!>
!> @details
!> The entry `GROUPX(j,k)` stores the MPI group containing all processes
!> with fixed \f$(y,z)\f$ coordinates and varying \f$x\f$ coordinate.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: GROUPX(NRPROCYMAX, NRPROCZMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief MPI groups corresponding to fibres parallel to the second
!> process-grid direction.
!>
!> @details
!> The entry `GROUPY(i,k)` stores the MPI group containing all processes
!> with fixed \f$(x,z)\f$ coordinates and varying \f$y\f$ coordinate.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: GROUPY(NRPROCXMAX, NRPROCZMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief MPI groups corresponding to fibres parallel to the third
!> process-grid direction.
!>
!> @details
!> The entry `GROUPZ(i,j)` stores the MPI group containing all processes
!> with fixed \f$(x,y)\f$ coordinates and varying \f$z\f$ coordinate.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: GROUPZ(NRPROCXMAX, NRPROCYMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Table of communicators for fibres parallel to the first
!> process-grid direction.
!>
!> @details
!> The entry `COMMXALL(j,k)` stores the communicator spanning all
!> processes with fixed \f$(y,z)\f$ coordinates.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMXALL(NRPROCYMAX, NRPROCZMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Table of communicators for fibres parallel to the second
!> process-grid direction.
!>
!> @details
!> The entry `COMMYALL(i,k)` stores the communicator spanning all
!> processes with fixed \f$(x,z)\f$ coordinates.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMYALL(NRPROCXMAX, NRPROCZMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Table of communicators for fibres parallel to the third
!> process-grid direction.
!>
!> @details
!> The entry `COMMZALL(i,j)` stores the communicator spanning all
!> processes with fixed \f$(x,y)\f$ coordinates.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMZALL(NRPROCXMAX, NRPROCYMAX)

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Communicator local to the current process for the fibre
!> parallel to the first process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMX

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Communicator local to the current process for the fibre
!> parallel to the second process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMY

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Communicator local to the current process for the fibre
!> parallel to the third process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: COMMZ

   PRIVATE :: GROUPX, GROUPY, GROUPZ
   PRIVATE :: COMMXALL, COMMYALL, COMMZALL
   PROTECTED :: processors, COMMX, COMMY, COMMZ

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Creates MPI groups and communicators for all one-dimensional
!> fibres of the logical process grid.
!>
!> @details
!> This routine constructs the full communication infrastructure used by
!> directional phases of the project. Starting from `MPI_COMM_WORLD`, it:
!> - obtains the world communicator group,
!> - fills the logical process-cube table \ref processors,
!> - creates MPI groups for all fibres parallel to the \f$x\f$,
!>   \f$y\f$, and \f$z\f$ directions,
!> - creates communicators associated with those groups,
!> - extracts the three communicators local to the current process and
!>   stores them in \ref COMMX, \ref COMMY, and \ref COMMZ.
!>
!> The communicator layout is consistent with the rank decomposition
!> supplied by module \ref parallelism.
!
! Output:
! -------
!> @param[out] mierr
!> Returned status code.
!
! Notes:
! ------
!> @note
!> The communicator tables are created for all fibres of the logical
!> process cube, not only for the fibre containing the current process.
!
!> @warning
!> The implementation assumes that the active process-grid dimensions do
!> not exceed `NRPROCXMAX`, `NRPROCYMAX`, and `NRPROCZMAX`.
!
!> @warning
!> MPI groups and communicators created here are not released within
!> this module in the current implementation.
!
!---------------------------------------------------------------------------
subroutine CreateCommunicators(mierr)
   use parallelism, ONLY: MYRANK, MYRANKX, MYRANKY, MYRANKZ, &
                           NRPROCX, NRPROCY, NRPROCZ
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   use mpi
   implicit none
!> @brief Returned status code.
   integer(kind=4), intent(out) :: mierr
!> @brief MPI group associated with `MPI_COMM_WORLD`.
   integer(kind=4) :: group_comm_world
!> @brief Temporary communicator handle returned by `mpi_comm_create`.
   integer(kind=4) :: comm_myrank_local
!> @brief Temporary list of global ranks for an \f$x\f$-aligned fibre.
   integer(kind=4) :: processors_X(NRPROCX)
!> @brief Temporary list of global ranks for a \f$y\f$-aligned fibre.
   integer(kind=4) :: processors_Y(NRPROCY)
!> @brief Temporary list of global ranks for a \f$z\f$-aligned fibre.
   integer(kind=4) :: processors_Z(NRPROCZ)
!> @brief Loop counters over logical process-grid coordinates.
   integer(kind=4) :: i, j, k
!> @brief MPI return code.
   integer(kind=4) :: ierr

#ifdef IPRINT
   write (*, *) MYRANK, 'NRPROC', NRPROC
#endif

   call mpi_comm_group(MPI_COMM_WORLD, group_comm_world, ierr)

   if (ierr /= 0) then
      write (*, *) MYRANK, ': main: error calling mpi_comm_group!'
      STOP 4
   end if
#ifdef IPRINT
   write (*, *) MYRANK, 'got group', group_comm_world
#endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         do k = 1, NRPROCZ
            processors(i, j, k) = (i - 1) + (j - 1)*NRPROCX + (k - 1)*NRPROCX*NRPROCY
         end do
      end do
   end do

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         processors_Z(1:NRPROCZ) = processors(i, j, 1:NRPROCZ)
         call mpi_group_incl(group_comm_world, NRPROCZ, processors_Z, GROUPZ(i, j), ierr)
         if (ierr /= 0) then
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for Z', i, j
            STOP 4
         end if
      end do
   end do
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         processors_Y(1:NRPROCY) = processors(i, 1:NRPROCY, k)
         call mpi_group_incl(group_comm_world, NRPROCY, processors_Y, GROUPY(i, k), ierr)
         if (ierr /= 0) then
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for Y', i, k
            STOP 4
         end if
      end do
   end do
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         processors_X(1:NRPROCX) = processors(1:NRPROCX, j, k)
         call mpi_group_incl(group_comm_world, NRPROCX, processors_X, GROUPX(j, k), ierr)
         if (ierr /= 0) then
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_group_incl for X', j, k
            STOP 4
         end if
      end do
   end do

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
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for Z', i, j
            STOP 4
         end if
      end do
   end do
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         call mpi_comm_create(MPI_COMM_WORLD, GROUPY(i, k), comm_myrank_local, ierr)
         COMMYALL(i, k) = comm_myrank_local
         if (ierr /= 0) then
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for Y', i, k
            STOP 4
         end if
      end do
   end do
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         call mpi_comm_create(MPI_COMM_WORLD, GROUPX(j, k), comm_myrank_local, ierr)
         COMMXALL(j, k) = comm_myrank_local
         if (ierr /= 0) then
            write (ERROR_UNIT, *) MYRANK, ': main: error calling mpi_com_create for X', j, k
            STOP 4
         end if
      end do
   end do
#ifdef IPRINT
   call PrintCommunicators
#endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)

   ! extract local communicators
   COMMX = COMMXALL(myranky + 1, myrankz + 1)
   COMMY = COMMYALL(myrankx + 1, myrankz + 1)
   COMMZ = COMMZALL(myrankx + 1, myranky + 1)

#ifdef IPRINT
   write (*, *) PRINTRANK, 'COMMX(Y,Z)', COMMX, COMMY, COMMZ
#endif

   mierr = 0

end subroutine CreateCommunicators

!!!!! dodac czyszczenie komunikatorow

end module communicators