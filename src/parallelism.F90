!------------------------------------------------------------------------------
!
! MODULE: parallelism
!
! DESCRIPTION:
!> @file parallelism.F90
!> @brief Module providing MPI-related process-grid configuration and
!> decomposition utilities.
!>
!> @details
!> This module groups the infrastructure routines and global variables
!> responsible for:
!> - initialization and finalization of MPI,
!> - storage of global and directional process-rank metadata,
!> - conversion between linear and three-dimensional process indices,
!> - computation of owned index ranges and corresponding element ranges,
!> - construction of dimension and shift vectors used by gather/scatter
!>   communication layers.
!>
!> Within the project, this module provides the low-level parallel
!> metadata consumed by higher-level components such as:
!> - \ref ADSS for workflow orchestration,
!> - \ref Setup for decomposition-aware setup containers,
!> - MPI communication wrappers and redistribution utilities.
!
!------------------------------------------------------------------------------
module parallelism

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Global MPI rank of the current process in `MPI_COMM_WORLD`.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: MYRANK

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Total number of MPI processes in `MPI_COMM_WORLD`.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: NRPROC

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Cartesian coordinate of the current process in the first
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: MYRANKX

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Cartesian coordinate of the current process in the second
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: MYRANKY

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Cartesian coordinate of the current process in the third
!> process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: MYRANKZ

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Total number of processes in the first process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: NRPROCX

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Total number of processes in the second process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: NRPROCY

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Total number of processes in the third process-grid direction.
!
!---------------------------------------------------------------------------
   integer(kind=4) :: NRPROCZ

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Character representation of the current process rank used in
!> diagnostic output.
!
!---------------------------------------------------------------------------
   character(len=7) :: PRINTRANK

   PROTECTED :: MYRANK, NRPROC, NRPROCX, NRPROCY, NRPROCZ, PRINTRANK

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Initializes MPI and the global process-grid metadata of the
!> module.
!>
!> @details
!> This routine stores the requested logical process-grid dimensions,
!> initializes MPI, queries the global rank and process count in
!> `MPI_COMM_WORLD`, and computes the three-dimensional process-grid
!> coordinates of the current rank through \ref Decompose.
!>
!> The process-grid dimensions are stored in the module-global variables
!> \ref NRPROCX, \ref NRPROCY, and \ref NRPROCZ. The global rank and
!> communicator size are stored in \ref MYRANK and \ref NRPROC.
!>
!> If any of the MPI initialization calls fails, the routine writes an
!> error message to `ERROR_UNIT` and terminates execution.
!
! Input:
! ------
!> @param[in] procx
!> Number of processes in the first process-grid direction.
!>
!> @param[in] procy
!> Number of processes in the second process-grid direction.
!>
!> @param[in] procz
!> Number of processes in the third process-grid direction.
!
! Output:
! -------
!> @param[out] ierr
!> Returned status code.
!
! Notes:
! ------
!> @note
!> The actual MPI return codes are accumulated internally in the local
!> variables `i1`, `i2`, and `i3`.
!
!> @warning
!> The routine does not verify that
!> \f$\text{procx}\cdot\text{procy}\cdot\text{procz}=\text{NRPROC}\f$.
!
!---------------------------------------------------------------------------
subroutine InitializeParallelism(procx, procy, procz, ierr)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   use mpi
   implicit none
!> @brief Numbers of processes in the three process-grid directions.
   integer(kind=4), intent(in) :: procx, procy, procz
!> @brief Returned status code.
   integer(kind=4), intent(out) :: ierr
!> @brief Auxiliary character buffer reserved for diagnostic use.
   character(4) :: buffer
!> @brief MPI return codes from initialization and rank/size queries.
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a linear MPI rank into three-dimensional
!> process-grid coordinates.
!>
!> @details
!> This routine maps the linear rank \p rank to the process-grid
!> coordinates \p rankx, \p ranky, and \p rankz using the logical grid
!> dimensions stored in \ref NRPROCX, \ref NRPROCY, and \ref NRPROCZ.
!>
!> The ordering convention is tensor-product lexicographic ordering with:
!> - the first direction changing fastest,
!> - the third direction changing slowest.
!>
!> Equivalently, the linear rank is interpreted according to
!>
!> \f[
!> \text{rank} = (\text{rankz}\cdot \text{NRPROCY} + \text{ranky})\cdot
!> \text{NRPROCX} + \text{rankx}.
!> \f]
!
! Input:
! ------
!> @param[in] rank
!> Linear rank of the process in `MPI_COMM_WORLD`.
!
! Output:
! -------
!> @param[out] rankx
!> Coordinate of the process in the first direction.
!>
!> @param[out] ranky
!> Coordinate of the process in the second direction.
!>
!> @param[out] rankz
!> Coordinate of the process in the third direction.
!
! Notes:
! ------
!> @note
!> The component order is `(X,Y,Z)` in the interface, while the linear
!> storage order corresponds to `Z` slowest and `X` fastest.
!
!---------------------------------------------------------------------------
subroutine Decompose(rank, rankx, ranky, rankz)
   implicit none
!> @brief Linear rank to be decomposed.
   integer(kind=4), intent(in) :: rank
!> @brief Cartesian coordinates of the process in the logical grid.
   integer(kind=4), intent(out) :: rankx, ranky, rankz

   rankz = rank/(NRPROCX*NRPROCY)
   ranky = (rank - rankz*(NRPROCX*NRPROCY))/NRPROCX
   rankx = rank - rankz*(NRPROCX*NRPROCY)
   rankx = rankx - ranky*NRPROCX

end subroutine Decompose

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the linear MPI rank associated with
!> three-dimensional process-grid coordinates.
!>
!> @details
!> This function performs the inverse mapping of \ref Decompose.
!> Given the process-grid coordinates \p rankx, \p ranky, and \p rankz,
!> it returns the linear rank in the tensor-product ordering in which
!> the first direction changes fastest and the third direction changes
!> slowest.
!
! Input:
! ------
!> @param[in] rankx
!> Coordinate in the first process-grid direction.
!>
!> @param[in] ranky
!> Coordinate in the second process-grid direction.
!>
!> @param[in] rankz
!> Coordinate in the third process-grid direction.
!
! Output:
! -------
!> @return rank
!> Linear rank corresponding to the supplied process-grid coordinates.
!
!---------------------------------------------------------------------------
function LinearIndex(rankx, ranky, rankz) result(rank)
   implicit none
   !> @brief Cartesian process-grid coordinates.
   integer(kind=4), intent(in) :: rankx, ranky, rankz
   !> @brief Returned linear rank.
   integer(kind=4) :: rank

   rank = rankz
   rank = rank*NRPROCY + ranky
   rank = rank*NRPROCX + rankx

end function LinearIndex

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the owned basis-function range and corresponding
!> element range for one process in one coordinate direction.
!>
!> @details
!> This routine determines the local interval of globally indexed basis
!> functions assigned to the processor identified by \p rank in a
!> one-dimensional partition with \p nrproc processors.
!>
!> The routine computes:
!> - \p nrcpp, the nominal number of columns per processor,
!> - \p ibeg and \p iend, the owned basis-function range,
!> - \p mine and \p maxe, the element range overlapping that ownership.
!>
!> The element count is derived as
!> \f$\text{elems} = n + 1 - p\f$, consistent with the project
!> convention for open spline spaces.
!
! Input:
! ------
!> @param[in] rank
!> Rank of the current processor in the selected direction.
!>
!> @param[in] nrproc
!> Number of processors in the selected direction.
!>
!> @param[in] n
!> Number of basis functions minus one in the selected direction.
!>
!> @param[in] p
!> Polynomial degree in the selected direction.
!
! Output:
! -------
!> @param[out] nrcpp
!> Nominal number of columns per processor.
!>
!> @param[out] ibeg
!> First owned basis-function index, using one-based indexing.
!>
!> @param[out] iend
!> Last owned basis-function index, using one-based indexing.
!>
!> @param[out] mine
!> First local element index overlapping the owned basis range.
!>
!> @param[out] maxe
!> Last local element index overlapping the owned basis range.
!
! Notes:
! ------
!> @note
!> The last processor in the direction may own fewer columns than
!> \p nrcpp.
!
!---------------------------------------------------------------------------
subroutine ComputeEndpoints(rank, nrproc, n, p, nrcpp, ibeg, iend, mine, maxe)
   implicit none
!> @brief Directional processor rank, processor count, basis size, and degree.
   integer(kind=4), intent(in) :: rank, nrproc, n, p
!> @brief Computed ownership and overlapping element-range data.
   integer(kind=4), intent(out) :: nrcpp, ibeg, iend, mine, maxe
!> @brief Number of elements in the selected direction.
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds dimension and shift vectors for directional gather/scatter
!> communication.
!>
!> @details
!> This routine allocates and fills the arrays \p dims and \p shifts used
!> to describe the sizes and offsets of processor-local slices in one
!> selected direction.
!>
!> The local slice size is computed as:
!>
!> \f[
!> \text{slice\_size} = \text{owned\_columns} \times \text{stride},
!> \f]
!>
!> where \p stride represents the size of the full two-dimensional
!> transverse slice. The resulting vectors are suitable for MPI gather or
!> scatter calls that require counts and displacements.
!
! Input:
! ------
!> @param[in] nrcpp
!> Nominal number of columns per processor.
!>
!> @param[in] stride
!> Size of the full transverse slice associated with one column.
!>
!> @param[in] n
!> Number of basis functions minus one in the selected direction.
!>
!> @param[in] nrproc
!> Number of processors in the selected direction.
!
! Output:
! -------
!> @param[out] dims
!> Array of local slice sizes for all processors.
!>
!> @param[out] shifts
!> Array of offsets of local slices in gathered storage.
!
! Notes:
! ------
!> @note
!> The last processor may receive a slice size smaller than
!> \p nrcpp*stride.
!
!---------------------------------------------------------------------------
subroutine FillDimVector(dims, shifts, nrcpp, stride, n, nrproc)
   implicit none
!> @brief Partition and problem-size data used to build the vectors.
   integer(kind=4), intent(in) :: nrcpp, stride, n, nrproc
!> @brief Output arrays of slice sizes and offsets.
   integer(kind=4), allocatable, dimension(:), intent(out) :: dims, shifts
!> @brief Loop counter.
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

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Finalizes MPI and cleans the parallel runtime environment.
!>
!> @details
!> This routine calls `mpi_finalize` and returns its status code through
!> \p ierr.
!
! Output:
! -------
!> @param[out] ierr
!> Returned status code from MPI finalization.
!
!---------------------------------------------------------------------------
subroutine Cleanup_Parallelism(ierr)
   use mpi
   implicit none
!> @brief Returned status code from MPI finalization.
   integer(kind=4), intent(out) :: ierr

   call mpi_finalize(ierr)
end subroutine Cleanup_Parallelism

end module parallelism