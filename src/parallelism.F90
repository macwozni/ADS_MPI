!------------------------------------------------------------------------------
!
! MODULE: parallelism
!
! DESCRIPTION
!> This module contains all functionality associated parallelism configuration.
!
!------------------------------------------------------------------------------

module parallelism

   implicit none

!> Rank of this processor
   integer(kind=4) :: MYRANK

!> Number of processors
   integer(kind=4) :: NRPROC

!> Integer coordinates of processor along X
   integer(kind=4) :: MYRANKX
!> Integer coordinates of processor along Y
   integer(kind=4) :: MYRANKY
!> Integer coordinates of processor along Z
   integer(kind=4) :: MYRANKZ

!> Total number of processors along X
   integer(kind=4) :: NRPROCX
!> Total number of processors along Y
   integer(kind=4) :: NRPROCY
!> Total number of processors along Z
   integer(kind=4) :: NRPROCZ

!> Rank of this processor converted to a string
   character(len=7) :: PRINTRANK

   PROTECTED :: MYRANK, NRPROC, NRPROCX, NRPROCY, NRPROCZ, PRINTRANK

contains

! -------------------------------------------------------------------
!> @brief
!> Initializes MPI communicators and global variables of this module.
!
! Input:
! ------
!> @param[in]  procx - number of processors along X axis
!> @param[in]  procy - number of processors along Y axis
!> @param[in]  procz - number of processors along Z axis
!
! Output:
! -------
!> @param[out]  ierr - error code
! -------------------------------------------------------------------
   subroutine InitializeParallelism(procx, procy, procz, ierr)
      use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT ! access computing environment
      use mpi
      implicit none
      integer(kind=4), intent(in) :: procx, procy, procz
      integer(kind=4), intent(out) :: ierr
      character(4) :: buffer
      integer(kind=4) :: i1, i2, i3

      NRPROCX = procx
      NRPROCY = procy
      NRPROCZ = procz

      ! Initialize MPI
      call mpi_init(i1)
      call mpi_comm_size(MPI_COMM_WORLD, NRPROC, i2)
      call mpi_comm_rank(MPI_COMM_WORLD, MYRANK, i3)

      if ((i1 + i2 + i3) /= 0) then
         write (ERROR_UNIT, *) MYRANK, ': main: error initializing MPI!'
         STOP 4
      end if

      call Decompose(MYRANK, MYRANKX, MYRANKY, MYRANKZ)

      ierr = 0

   end subroutine InitializeParallelism

! -------------------------------------------------------------------
!> @brief
!> Based on the linear index (process rank) computes its coordinates
!> in 3D cube \f$ (\text{NRPROCX} \times \text{NRPROCY} \times \text{NRPROCZ}) \f$. 
!
! Input:
! ------
!> @param[in]  rank    - linear rank of the process
!
! Output:
! -------
!> @param[out]  rankx   - x coordinate
!> @param[out]  ranky   - y coordinate
!> @param[out]  rankz   - z coordinate
!>
!> Order of components (from slowest changing): Z, Y, X.
!> An example for \f$ 1 \times 1 \times 1 \f$   
!> Rank  | Coords
!> ------------- | -------------
!> \f$ 0 \f$    | \f$ (0, 0, 0) \f$
!> \f$ 1 \f$    | \f$ (0, 0, 1) \f$
!> \f$ 2 \f$    | \f$ (0, 1, 0) \f$
!> \f$ 3 \f$    | \f$ (0, 1, 1) \f$
!> \f$ 4 \f$    | \f$ (1, 0, 0) \f$
!> etc.         | etc.
! -------------------------------------------------------------------
   subroutine Decompose(rank, rankx, ranky, rankz)
      implicit none
      integer(kind=4), intent(in) :: rank
      integer(kind=4), intent(out) :: rankx, ranky, rankz

      rankz = rank/(NRPROCX*NRPROCY)
      ranky = (rank - rankz*(NRPROCX*NRPROCY))/NRPROCX
      rankx = rank - rankz*(NRPROCX*NRPROCY)
      rankx = rankx - ranky*NRPROCX

   end subroutine Decompose

! -------------------------------------------------------------------
!> @brief
!> Computes global, linear index, based on coordinates in 3D cube.
!>
! Input:
! ------
!> @param[in]  rankx   - x coordinate
!> @param[in]  ranky   - y coordinate
!> @param[in]  rankz   - z coordinate
!
! Output:
! -------
!> @return rank - linear index based on (Z, Y, X) lexicographic order (Z changes slowest)
! -------------------------------------------------------------------
   function LinearIndex(rankx, ranky, rankz) result(rank)
      implicit none
      integer(kind=4), intent(in) :: rankx, ranky, rankz
      integer(kind=4) :: rank

      rank = rankz
      rank = rank*NRPROCY + ranky
      rank = rank*NRPROCX + rankx

   end function LinearIndex

! -------------------------------------------------------------------
!> @brief
!> Calculates the range of the direction that is assigned to processor
!> with specified rank.
!
! Input:
! ------
!> @param[in]  rank    - rank of the current process in this direction
!> @param[in]  nrproc  - number of processors for this direction
!> @param[in]  n       - size of the problem
!> @param[in]  p       - order of polynomial basis
!
! Output:
! -------
!> @param[out]  nrcpp   - number of columns per processor
!> @param[out]  ibeg    - index of first assigned slice
!> @param[out]  iend    - index of last assigned slice
!> @param[out]  mine    - index of first element corresponding to the assigned slice
!> @param[out]  maxe    - index of last element corresponding to the assigned slice
! -------------------------------------------------------------------
   subroutine ComputeEndpoints(rank, nrproc, n, p, nrcpp, ibeg, iend, mine, maxe)
      implicit none
      integer(kind=4), intent(in) :: rank, nrproc, n, p
      integer(kind=4), intent(out) :: nrcpp, ibeg, iend, mine, maxe
      integer(kind=4) :: elems

      elems = n + 1 - p
      nrcpp = (n + 1 + nrproc - 1)/nrproc
      ibeg = nrcpp*rank + 1

      if (rank == nrproc - 1) then
         iend = n + 1
      else
         iend = nrcpp*(rank + 1)
      end if

      mine = max(ibeg - p - 1, 1)
      maxe = min(iend, elems)

   end subroutine ComputeEndpoints



! -------------------------------------------------------------------
!> @brief
!> Calculates sizes and ranges of slices for each processor in
!> given direction.
!
!
! Input:
! ------
!> @param[in] dims     - sizes of 'slices' of dimension
!> @param[in] shifts   - offsets of slices
!> @param[in] nrcpp    - number of columns per processor
!> @param[in] stride   - size of full, 3d slice
!
!
! Output:
! -------
!> @param[out] n        - size of the problem
!> @param[out] nrproc   - number of processors for this direction
! -------------------------------------------------------------------
   subroutine FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)
      implicit none
      integer(kind=4), intent(in) :: nrcpp, stride, n, nrproc
      integer(kind=4), allocatable, dimension(:), intent(out) :: dims, shifts
      integer(kind=4) :: i

      allocate (dims(nrproc))
      allocate (shifts(nrproc))

      shifts = 0
      dims = 0

      do i = 1, nrproc - 1
         dims(i) = nrcpp*stride
         if (i > 1) shifts(i) = shifts(i - 1) + dims(i - 1)
      end do

      if (nrproc > 1) then
         dims(nrproc) = ((n + 1) - nrcpp*(nrproc - 1))*stride
         shifts(nrproc) = shifts(nrproc - 1) + dims(nrproc - 1)
      else
         dims(1) = (n + 1)*stride
         shifts(1) = 0
      end if

   end subroutine FillDimVector


! -------------------------------------------------------------------
!> @brief
!> Cleans all parallelism structures, and finilizes MPI
!
!
! Input:
! ------
!
! Output:
! -------
!> @param[out] ierr    - error code
! -------------------------------------------------------------------
   subroutine Cleanup_Parallelism(ierr)
      implicit none
      integer(kind=4), intent(out) :: ierr

      call mpi_finalize(ierr)
   end subroutine CleanParallelism

end module parallelism
