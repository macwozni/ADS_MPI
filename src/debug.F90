!------------------------------------------------------------------------------
!
! MODULE: DEBUGG
!
! DESCRIPTION:
!> @file debug.F90
!> @brief Module providing diagnostic and validation utilities for
!> decomposition, redistribution, and communicator inspection.
!>
!> @details
!> This module groups auxiliary debugging routines used across the
!> project to:
!> - inspect data redistributed to neighbouring processes,
!> - validate dimension vectors used in MPI gather/scatter operations,
!> - print domain-decomposition metadata,
!> - test whether tensor-product indices belong to a local ownership
!>   range,
!> - display MPI groups and communicators used by the parallel runtime.
!>
!> The routines are intended for development-time verification of the
!> infrastructure shared by:
!> - \ref parallelism,
!> - \ref communicators,
!> - \ref Setup,
!> - the ADS workflow and redistribution layers.
!
!------------------------------------------------------------------------------
module DEBUGG

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Prints neighbouring redistribution buffers stored in the local
!> coefficient-exchange array.
!>
!> @details
!> This routine displays the contents of the neighbour-data buffer
!> \p ads_data%R for all process neighbours surrounding the current
!> process in the logical process cube.
!>
!> For each neighbouring process block, the routine:
!> - computes the owned ranges in all three directions using
!>   \ref ComputeEndpoints,
!> - reconstructs the logical shape of the redistributed block,
!> - reshapes the corresponding one-dimensional storage slice of
!>   \p ads_data%R into a three-dimensional array,
!> - writes the result to standard output.
!>
!> The printed data are useful when verifying redistribution or halo-like
!> exchange stages preceding local spline reconstruction.
!
! Input:
! ------
!> @param[inout] ads
!> ADS setup structure providing decomposition metadata.
!>
!> @param[in] ads_data
!> Runtime-data structure containing redistributed neighbour buffers.
!
! Notes:
! ------
!> @note
!> The routine prints only the neighbourhood within one process step in
!> each direction around the current process.
!
!---------------------------------------------------------------------------
subroutine PrintDistributedData(ads, ads_data)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
                           NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
   implicit none
!> @brief ADS setup structure containing decomposition metadata.
   type(ADS_setup), intent(inout) :: ads
!> @brief Runtime-data structure holding redistributed neighbour blocks.
   type(ADS_compute_data), intent(in) :: ads_data
!> @brief Loop counters over neighbouring process-grid coordinates.
   integer(kind=4) :: i, j, k
!> @brief Owned basis-function bounds of a selected neighbour block.
   integer(kind=4) :: obegx, oendx, obegy, oendy, obegz, oendz
!> @brief Auxiliary element-range outputs from `ComputeEndpoints`.
   integer(kind=4) :: mine, maxe

   write (*, *) PRINTRANK, 'R:'

   do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
      do j = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
         do k = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
            write (*, *) '(i,j,k)=', i + 1, j + 1, k + 1

            call ComputeEndpoints(i - 1, NRPROCX, ads%n(1), ads%p(1), ads%nrcpp(1), obegx, oendx, mine, maxe)
            call ComputeEndpoints(j - 1, NRPROCY, ads%n(2), ads%p(2), ads%nrcpp(2), obegy, oendy, mine, maxe)
            call ComputeEndpoints(k - 1, NRPROCZ, ads%n(3), ads%p(3), ads%nrcpp(3), obegz, oendz, mine, maxe)

            write (*, *) reshape( &
               ads_data%R(:, i - MYRANKX + 1, j - MYRANKY + 1, k - MYRANKZ + 1), &
               (/oendz - obegz + 1, oendx - obegx + 1, oendy - obegy + 1/))
         end do
      end do
   end do

   write (*, *) '----'

end subroutine PrintDistributedData

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Performs a consistency check of dimension vectors used by
!> directional MPI redistribution.
!>
!> @details
!> This routine verifies that the dimension vectors associated with the
!> three coordinate directions sum to the expected global slice sizes.
!>
!> The expected totals are:
!> - for \p dimensionsX: \f$(n_x+1)\,s_y\,s_z\f$,
!> - for \p dimensionsY: \f$(n_y+1)\,s_x\,s_z\f$,
!> - for \p dimensionsZ: \f$(n_z+1)\,s_x\,s_y\f$.
!>
!> If any inconsistency is detected, diagnostic information is printed
!> and program execution is stopped.
!
! Input:
! ------
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] s
!> Local owned sizes in the three directions.
!>
!> @param[in] nrcpp
!> Nominal numbers of columns per processor in the three directions.
!>
!> @param[in] dimensionsX
!> Dimension vector for redistribution along the first direction.
!>
!> @param[in] dimensionsY
!> Dimension vector for redistribution along the second direction.
!>
!> @param[in] dimensionsZ
!> Dimension vector for redistribution along the third direction.
!
! Notes:
! ------
!> @warning
!> On failure, the routine terminates execution with `stop`.
!
!---------------------------------------------------------------------------
subroutine ValidateDimensions(n, s, nrcpp, &
                              dimensionsX, dimensionsY, dimensionsZ)
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK
   use mpi
   implicit none
!> @brief Numbers of basis functions minus one in the three directions.
   integer(kind=4), intent(in), dimension(3) :: n
!> @brief Local owned sizes in the three directions.
   integer(kind=4), intent(in), dimension(3) :: s
!> @brief Nominal numbers of columns per processor.
   integer(kind=4), intent(in), dimension(3) :: nrcpp
!> @brief Dimension vector for redistribution along the first direction.
   integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsX
!> @brief Dimension vector for redistribution along the second direction.
   integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsY
!> @brief Dimension vector for redistribution along the third direction.
   integer(kind=4), intent(in), allocatable, dimension(:) :: dimensionsZ

   integer(kind=4) :: i, k
   !> @brief Loop counter and accumulated dimension total.

   k = 0
   do i = 1, NRPROCX
      k = k + dimensionsX(i)
   end do
   if (k /= (n(1) + 1)*s(2)*s(3)) then
      write (*, *) PRINTRANK, 'problem with dimensionsX', dimensionsX
      write (*, *) PRINTRANK, 'nx+1', n(1) + 1
      write (*, *) PRINTRANK, 'sy', s(2)
      write (*, *) PRINTRANK, 'sz', s(3)
      write (*, *) PRINTRANK, 'nrcppx', nrcpp(1)
      stop
   end if

   k = 0
   do i = 1, NRPROCY
      k = k + dimensionsY(i)
   end do
   if (k /= (n(2) + 1)*s(1)*s(3)) then
      write (*, *) PRINTRANK, 'problem with dimensionsY', dimensionsY
      write (*, *) PRINTRANK, 'n+1', n(2) + 1
      write (*, *) PRINTRANK, 'sx', s(1)
      write (*, *) PRINTRANK, 'sz', s(3)
      stop
   end if

   k = 0
   do i = 1, NRPROCZ
      k = k + dimensionsZ(i)
   end do
   if (k /= (n(3) + 1)*s(1)*s(2)) then
      write (*, *) PRINTRANK, 'problem with dimensionsZ', dimensionsZ
      write (*, *) PRINTRANK, 'n+1', n(3) + 1
      write (*, *) PRINTRANK, 'sx', s(1)
      write (*, *) PRINTRANK, 'sy', s(2)
      stop
   end if

end subroutine ValidateDimensions

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Prints a summary of the computed local domain decomposition.
!>
!> @details
!> This routine writes the current process-grid coordinates, global
!> process-grid dimensions, global basis sizes, nominal columns per
!> processor, and the local owned ranges in all three coordinate
!> directions.
!>
!> It is intended as a lightweight inspection tool for verifying the
!> output of decomposition routines from module \ref parallelism.
!
! Input:
! ------
!> @param[in] n
!> Numbers of basis functions minus one in the three directions.
!>
!> @param[in] nrcpp
!> Nominal numbers of columns per processor in the three directions.
!>
!> @param[in] ibeg
!> First owned basis-function indices in the three directions.
!>
!> @param[in] iend
!> Last owned basis-function indices in the three directions.
!
!---------------------------------------------------------------------------
subroutine PrintDecompositionInfo(n, nrcpp, ibeg, iend)
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK, &
                           MYRANKX, MYRANKY, MYRANKZ
   implicit none
!> @brief Numbers of basis functions minus one in the three directions.
   integer(kind=4), intent(in), dimension(3) :: n
!> @brief Nominal numbers of columns per processor in the three directions.
   integer(kind=4), intent(in), dimension(3) :: nrcpp
!> @brief First owned basis-function indices.
   integer(kind=4), intent(in), dimension(3) :: ibeg
!> @brief Last owned basis-function indices.
   integer(kind=4), intent(in), dimension(3) :: iend

   write (*, *) PRINTRANK, 'MYRANKX,MYRANKY,MYRANKZ', MYRANKX, MYRANKY, MYRANKZ
   write (*, *) PRINTRANK, 'NRPROCX,NRPROCY,NRPROCZ', NRPROCX, NRPROCY, NRPROCZ
   write (*, *) PRINTRANK, 'nx+1', n(1) + 1
   write (*, *) PRINTRANK, 'ny+1', n(2) + 1
   write (*, *) PRINTRANK, 'nz+1', n(3) + 1
   write (*, *) PRINTRANK, 'nrcppx,nrcppy,nrcppz', nrcpp(1), nrcpp(2), nrcpp(3)
   write (*, *) PRINTRANK, 'ibegx,iendx', ibeg(1), iend(1)
   write (*, *) PRINTRANK, 'ibegy,iendy', ibeg(2), iend(2)
   write (*, *) PRINTRANK, 'ibegz,iendz', ibeg(3), iend(3)
end subroutine PrintDecompositionInfo

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Tests whether a three-dimensional index belongs to a local
!> ownership range.
!>
!> @details
!> This logical function checks whether the coordinate triple \p ind lies
!> within the half-open ownership convention used by the project, namely
!> between `ibeg-1` and `iend-1` inclusively in each direction.
!>
!> The function returns `.TRUE.` only if the index lies inside the local
!> range in all three coordinate directions.
!
! Input:
! ------
!> @param[in] ind
!> Three-dimensional index to be tested.
!>
!> @param[in] ibeg
!> First owned indices in the three directions.
!>
!> @param[in] iend
!> Last owned indices in the three directions.
!
! Output:
! -------
!> @return IndexInRange
!> Logical flag indicating whether the index lies within the supplied
!> range.
!
!---------------------------------------------------------------------------
logical function IndexInRange(ind, ibeg, iend)
   implicit none
!> @brief Three-dimensional index to be tested.
   integer(kind=4), dimension(3), intent(in) :: ind
!> @brief Lower and upper ownership bounds.
   integer(kind=4), dimension(3), intent(in) :: ibeg, iend

   IndexInRange = .true.
   if (ind(1) < ibeg(1) - 1 .or. ind(1) > iend(1) - 1) IndexInRange = .false.
   if (ind(2) < ibeg(2) - 1 .or. ind(2) > iend(2) - 1) IndexInRange = .false.
   if (ind(3) < ibeg(3) - 1 .or. ind(3) > iend(3) - 1) IndexInRange = .false.
end function IndexInRange

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Prints all communicator handles created for directional fibres
!> of the logical process grid.
!>
!> @details
!> This routine writes the communicator tables associated with the
!> \f$x\f$-, \f$y\f$-, and \f$z\f$-aligned fibres of the process cube.
!>
!> It is intended strictly for debugging communicator creation and for
!> inspecting the global communicator layout after initialization of
!> module \ref communicators.
!
! Notes:
! ------
!> @note
!> This routine is for diagnostic use only.
!
!---------------------------------------------------------------------------
subroutine PrintCommunicators
   use parallelism, ONLY: PRINTRANK, NRPROCX, NRPROCY, NRPROCZ
   implicit none
!> @brief Loop counters over process-grid coordinates.
   integer(kind=4) :: i, j, k

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         write (*, *) PRINTRANK, 'COMMZALL(', i, j, ')', COMMZALL(i, j)
      end do
   end do
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         write (*, *) PRINTRANK, 'COMMYALL(', i, k, ')', COMMYALL(i, k)
      end do
   end do
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         write (*, *) PRINTRANK, 'COMMXALL(', j, k, ')', COMMXALL(j, k)
      end do
   end do
end subroutine PrintCommunicators

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Prints MPI groups later used to create directional fibre
!> communicators.
!>
!> @details
!> This routine writes the MPI-group handles associated with the
!> \f$x\f$-, \f$y\f$-, and \f$z\f$-aligned fibres of the logical process
!> grid.
!>
!> It is intended strictly for debugging the communicator-construction
!> phase and verifying the group layout before or after communicator
!> creation.
!
! Notes:
! ------
!> @note
!> This routine is for diagnostic use only.
!
!---------------------------------------------------------------------------
subroutine PrintGroups
   use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK
   implicit none
!> @brief Loop counters over process-grid coordinates.
   integer(kind=4) :: i, j, k

   do i = 1, NRPROCX
      do j = 1, NRPROCY
         write (*, *) PRINTRANK, 'GROUPZ(', i, j, ')', GROUPZ(i, j)
      end do
   end do
   do i = 1, NRPROCX
      do k = 1, NRPROCZ
         write (*, *) PRINTRANK, 'GROUPY(', i, k, ')', GROUPY(i, k)
      end do
   end do
   do j = 1, NRPROCY
      do k = 1, NRPROCZ
         write (*, *) PRINTRANK, 'GROUPX(', j, k, ')', GROUPX(j, k)
      end do
   end do
end subroutine PrintGroups

end module DEBUGG