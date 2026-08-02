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
   private :: BalancedBlock

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
!> communicator size are stored in \ref MYRANK and \ref NRPROC. A
!> zero-padded character representation of the rank is stored in
!> \ref PRINTRANK for diagnostic output.
!>
!> If any MPI initialization call fails, or if the requested process grid
!> is inconsistent with `MPI_COMM_WORLD`, the routine writes an error
!> message to `ERROR_UNIT` and terminates execution.
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
!> The routine terminates execution when
!> \f$\text{procx}\cdot\text{procy}\cdot\text{procz}\ne\text{NRPROC}\f$.
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
!> @brief Auxiliary character buffer used to format the process rank.
   character(len=7) :: buffer
!> @brief MPI return codes from initialization and rank/size queries.
   integer(kind=4) :: i1, i2, i3
!> @brief Remaining MPI ranks while checking the process-grid product safely.
   integer(kind=4) :: remaining_proc
!> @brief Whether the process-grid product equals the MPI world size.
   logical :: process_grid_matches

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

   if (procx <= 0 .or. procy <= 0 .or. procz <= 0) then
      if (MYRANK == 0) then
         write (ERROR_UNIT, *) 'Process-grid dimensions must be positive:', &
            procx, procy, procz
         flush (ERROR_UNIT)
      end if
      call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
      STOP 4
   end if

   process_grid_matches = .FALSE.
   remaining_proc = NRPROC
   if (mod(remaining_proc, procx) == 0) then
      remaining_proc = remaining_proc/procx
      if (mod(remaining_proc, procy) == 0) then
         remaining_proc = remaining_proc/procy
         process_grid_matches = remaining_proc == procz
      end if
   end if

   if (.NOT. process_grid_matches) then
      if (MYRANK == 0) then
         write (ERROR_UNIT, *) 'Process-grid size mismatch: dimensions =', &
            procx, procy, procz, 'MPI ranks =', NRPROC
         flush (ERROR_UNIT)
      end if
      call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
      STOP 4
   end if

   call Decompose(MYRANK, MYRANKX, MYRANKY, MYRANKZ)

   write (buffer, '(I0)') MYRANK
   PRINTRANK = repeat('0', max(0, 5 - len_trim(buffer)))//trim(buffer)

   ierr = 0

end subroutine InitializeParallelism

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Terminates the complete MPI job after a fatal distributed error.
!>
!> @details
!> Solver errors are synchronized before they reach problem drivers.  This
!> helper prevents a driver from using a partial result or overwriting the
!> original status during cleanup.  Rank zero reports the original code and
!> `MPI_Abort` terminates every rank, including ranks blocked in a collective.
!
!---------------------------------------------------------------------------
subroutine AbortOnError(status, operation)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   use mpi
   implicit none
   integer(kind=4), intent(in) :: status
   character(len=*), intent(in) :: operation
   integer(kind=4) :: abort_ierr

   if (status == 0) return

   if (MYRANK == 0) then
      write (ERROR_UNIT, '(A,A,I0)') trim(operation), &
         ' failed; aborting MPI job with original status ', status
      flush (ERROR_UNIT)
   end if
   call MPI_Abort(MPI_COMM_WORLD, status, abort_ierr)
   stop 5
end subroutine AbortOnError

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
!> @brief Computes one block of a balanced one-dimensional partition.
!>
!> @details
!> The \f$n+1\f$ globally indexed columns are split into contiguous blocks.
!> Every process receives either
!> \f$\lfloor(n+1)/\mathit{nrproc}\rfloor\f$ or one more column, with the
!> remainder assigned to the lowest process ranks.
!
!> @param[in] rank Zero-based rank in the selected direction.
!> @param[in] nrproc Number of processes in the selected direction.
!> @param[in] n Number of columns minus one.
!> @param[out] ibeg First owned one-based column index.
!> @param[out] iend Last owned one-based column index.
!
!---------------------------------------------------------------------------
pure subroutine BalancedBlock(rank, nrproc, n, ibeg, iend)
   implicit none
   integer(kind=4), intent(in) :: rank, nrproc, n
   integer(kind=4), intent(out) :: ibeg, iend
   integer(kind=4) :: base, owned, remainder

   base = (n + 1)/nrproc
   remainder = mod(n + 1, nrproc)

   owned = base
   if (rank < remainder) owned = owned + 1

   ibeg = rank*base + min(rank, remainder) + 1
   iend = ibeg + owned - 1

end subroutine BalancedBlock

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
!> - \p nrcpp, the maximum number of columns owned by one processor,
!> - \p ibeg and \p iend, the owned basis-function range,
!> - \p mine and \p maxe, the element range overlapping that ownership.
!>
!> The basis functions are partitioned into contiguous balanced blocks.
!> The first `mod(n + 1, nrproc)` processors receive one more basis
!> function than the remaining processors.
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
!> Maximum number of columns owned by one processor, equal to
!> `ceiling((n + 1)/nrproc)`.
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
!> Every processor owns either \p nrcpp or `nrcpp - 1` columns when
!> `n + 1 >= nrproc`; the larger blocks are assigned to lower ranks.
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
   call BalancedBlock(rank, nrproc, n, ibeg, iend)

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
!> Maximum number of columns per processor. Retained in the public
!> interface for compatibility; exact counts are derived from \p n and
!> \p nrproc.
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
!> The larger blocks are assigned to the lowest process ranks, consistently
!> with \ref ComputeEndpoints.
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
!> @brief Balanced block bounds for one process.
   integer(kind=4) :: block_begin, block_end

   allocate (dims(nrproc))
   allocate (shifts(nrproc))

   do i = 1, nrproc
      call BalancedBlock(i - 1, nrproc, n, block_begin, block_end)
      dims(i) = (block_end - block_begin + 1)*stride
      shifts(i) = (block_begin - 1)*stride
   end do

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
