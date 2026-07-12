!------------------------------------------------------------------------------
!
! MODULE: sparse
!
! DESCRIPTION:
!> @file sparse.F90
!> @brief Sparse-matrix assembly and conversion utilities for ADS.
!>
!> @details
!> The module stores sparse matrices row by row. Each row owns two compact
!> dynamic arrays:
!> - column indices,
!> - accumulated numerical values.
!>
!> Entries are not kept sorted during assembly. This keeps \ref add cheap:
!> the routine performs only a short row-local search and appends new
!> entries at the end of the row. Export routines sort each row on demand,
!> preserving the historical row-major, column-increasing traversal order
!> expected by the rest of the code and by MUMPS conversion tests.
!
!------------------------------------------------------------------------------
module sparse

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Row-local sparse storage.
!>
!> @details
!> The row stores unique column indices and accumulated values. The arrays
!> are intentionally unordered during assembly.
!
!---------------------------------------------------------------------------
   type sparse_row_storage
!> @brief Number of explicitly stored entries in the row.
      integer(kind=4) :: nnz = 0
!> @brief Allocated capacity of the row buffers.
      integer(kind=4) :: cap = 0
!> @brief Stored column indices.
      integer(kind=4), allocatable :: col(:)
!> @brief Stored numerical values.
      real(kind=8), allocatable :: val(:)
   end type sparse_row_storage

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Sparse matrix container.
!>
!> @details
!> The matrix uses zero-based row and column indices internally. The public
!> API keeps the historical pointer-based allocation contract so existing
!> callers do not need to change.
!
!---------------------------------------------------------------------------
   type sparse_matrix
!> @brief Number of rows.
      integer(kind=4) :: x = 0
!> @brief Number of columns.
      integer(kind=4) :: y = 0
!> @brief Number of explicitly stored unique sparse entries.
      integer(kind=8) :: total_entries = 0_8
!> @brief Row-wise sparse storage.
      type(sparse_row_storage), allocatable :: rows(:)
   end type sparse_matrix

   private :: sparse_row_storage
   private :: ensure_row_capacity, find_in_row, append_to_row
   private :: sorted_row_order, release_mumps_triplets
   private :: build_mumps_row_offsets, fill_mumps_triplets
   public  :: initialize_sparse, clear_matrix, add, add_existing
   public  :: to_dense_matrix, to_mumps_format, to_mumps_format_transposed
   public  :: sparse_matrix

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates and initializes an empty sparse matrix.
!
! Input:
! ------
!> @param[in] x
!> Number of rows.
!>
!> @param[in] y
!> Number of columns.
!
! Output:
! -------
!> @param[out] matrix
!> Newly allocated empty sparse matrix.
!
!---------------------------------------------------------------------------
subroutine initialize_sparse(x, y, matrix)
   implicit none
   integer(kind=4), intent(in) :: x
   integer(kind=4), intent(in) :: y
   type(sparse_matrix), pointer, intent(out) :: matrix

   allocate(matrix)
   matrix%x = x
   matrix%y = y
   matrix%total_entries = 0_8

   if (x > 0) then
      allocate(matrix%rows(0:x - 1))
   else
      allocate(matrix%rows(1:0))
   end if
end subroutine initialize_sparse

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Ensures that a row can store at least \p needed entries.
!
! Input:
! ------
!> @param[in] needed
!> Required minimum row capacity.
!
! Input/Output:
! -------------
!> @param[inout] row
!> Row whose storage may be reallocated.
!
!---------------------------------------------------------------------------
subroutine ensure_row_capacity(row, needed)
   implicit none
   type(sparse_row_storage), intent(inout) :: row
   integer(kind=4), intent(in) :: needed
   integer(kind=4) :: newcap
   integer(kind=4), allocatable :: newcol(:)
   real(kind=8), allocatable :: newval(:)

   if (row%cap >= needed) return

   if (row%cap == 0) then
      newcap = max(needed, 8)
   else
      newcap = row%cap
      do while (newcap < needed)
         newcap = 2 * newcap
      end do
   end if

   allocate(newcol(newcap))
   allocate(newval(newcap))

   if (row%nnz > 0) then
      newcol(1:row%nnz) = row%col(1:row%nnz)
      newval(1:row%nnz) = row%val(1:row%nnz)
   end if

   if (allocated(row%col)) deallocate(row%col)
   if (allocated(row%val)) deallocate(row%val)

   call move_alloc(newcol, row%col)
   call move_alloc(newval, row%val)
   row%cap = newcap
end subroutine ensure_row_capacity

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Finds a column index inside one row.
!
! Output:
! -------
!> @param[out] found
!> Flag indicating whether the column exists.
!>
!> @param[out] pos
!> Existing position when found, or append position otherwise.
!
!---------------------------------------------------------------------------
subroutine find_in_row(row, col, found, pos)
   implicit none
   type(sparse_row_storage), intent(in) :: row
   integer(kind=4), intent(in) :: col
   logical, intent(out) :: found
   integer(kind=4), intent(out) :: pos
   integer(kind=4) :: i

   found = .false.
   do i = 1, row%nnz
      if (row%col(i) == col) then
         found = .true.
         pos = i
         return
      end if
   end do

   pos = row%nnz + 1
end subroutine find_in_row

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Appends a new unique entry to a sparse row.
!
!---------------------------------------------------------------------------
subroutine append_to_row(row, col, val)
   implicit none
   type(sparse_row_storage), intent(inout) :: row
   integer(kind=4), intent(in) :: col
   real(kind=8), intent(in) :: val

   call ensure_row_capacity(row, row%nnz + 1)

   row%nnz = row%nnz + 1
   row%col(row%nnz) = col
   row%val(row%nnz) = val
end subroutine append_to_row

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Adds a contribution to a sparse matrix entry.
!>
!> @details
!> Repeated additions to the same \f$(x,y)\f$ entry are accumulated in
!> place and do not increase \ref sparse_matrix::total_entries.
!
!---------------------------------------------------------------------------
subroutine add(matrix, x, y, val)
   implicit none
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4), intent(in) :: x
   integer(kind=4), intent(in) :: y
   real(kind=8), intent(in) :: val
   logical :: found
   integer(kind=4) :: pos

   call find_in_row(matrix%rows(x), y, found, pos)

   if (found) then
      matrix%rows(x)%val(pos) = matrix%rows(x)%val(pos) + val
   else
      call append_to_row(matrix%rows(x), y, val)
      matrix%total_entries = matrix%total_entries + 1_8
   end if
end subroutine add

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Adds a contribution to an already allocated sparse entry.
!>
!> @details
!> This routine is intended for assembly kernels that build the sparse
!> pattern before the numerical phase. It never reallocates row storage and
!> never changes \ref sparse_matrix::total_entries.
!
!---------------------------------------------------------------------------
subroutine add_existing(matrix, x, y, val)
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   implicit none
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4), intent(in) :: x
   integer(kind=4), intent(in) :: y
   real(kind=8), intent(in) :: val
   logical :: found
   integer(kind=4) :: pos

   call find_in_row(matrix%rows(x), y, found, pos)

   if (.not. found) then
      write(ERROR_UNIT, *) 'add_existing called for a missing sparse entry:', x, y
      stop 1
   end if

   matrix%rows(x)%val(pos) = matrix%rows(x)%val(pos) + val
end subroutine add_existing

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Returns row-entry indices sorted by increasing column index.
!
!---------------------------------------------------------------------------
subroutine sorted_row_order(row, order)
   implicit none
   type(sparse_row_storage), intent(in) :: row
   integer(kind=4), allocatable, intent(out) :: order(:)
   integer(kind=4) :: i
   integer(kind=4) :: j
   integer(kind=4) :: key

   allocate(order(row%nnz))

   do i = 1, row%nnz
      order(i) = i
   end do

   do i = 2, row%nnz
      key = order(i)
      j = i - 1
      do while (j >= 1)
         if (row%col(order(j)) <= row%col(key)) exit
         order(j + 1) = order(j)
         j = j - 1
      end do
      order(j + 1) = key
   end do
end subroutine sorted_row_order

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Releases all storage owned by a sparse matrix.
!
!---------------------------------------------------------------------------
subroutine clear_matrix(matrix)
   implicit none
   type(sparse_matrix), pointer, intent(inout) :: matrix
   integer(kind=4) :: i

   if (.not. associated(matrix)) return

   if (allocated(matrix%rows)) then
      do i = lbound(matrix%rows, 1), ubound(matrix%rows, 1)
         if (allocated(matrix%rows(i)%col)) deallocate(matrix%rows(i)%col)
         if (allocated(matrix%rows(i)%val)) deallocate(matrix%rows(i)%val)
         matrix%rows(i)%nnz = 0
         matrix%rows(i)%cap = 0
      end do
      deallocate(matrix%rows)
   end if

   deallocate(matrix)
   nullify(matrix)
end subroutine clear_matrix

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to dense storage.
!
!---------------------------------------------------------------------------
subroutine to_dense_matrix(matrix, dmatrix)
   implicit none
   type(sparse_matrix), intent(in) :: matrix
   real(kind=8), dimension(0:matrix%x - 1, 0:matrix%y - 1), intent(out) :: dmatrix
   integer(kind=4) :: i
   integer(kind=4) :: k

   dmatrix = 0.0d0

   do i = 0, matrix%x - 1
      do k = 1, matrix%rows(i)%nnz
         dmatrix(i, matrix%rows(i)%col(k)) = matrix%rows(i)%val(k)
      end do
   end do
end subroutine to_dense_matrix

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Deallocates triplet arrays already attached to a MUMPS structure.
!
!---------------------------------------------------------------------------
subroutine release_mumps_triplets(mumps_par)
   implicit none
   include 'dmumps_struc.h'
   type(dmumps_struc), intent(inout) :: mumps_par

   if (associated(mumps_par%irn)) deallocate(mumps_par%irn)
   if (associated(mumps_par%jcn)) deallocate(mumps_par%jcn)
   if (associated(mumps_par%a))   deallocate(mumps_par%a)
end subroutine release_mumps_triplets

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds one-based row starts for deterministic MUMPS triplet export.
!
!---------------------------------------------------------------------------
subroutine build_mumps_row_offsets(matrix, row_offsets)
   implicit none
   type(sparse_matrix), pointer, intent(in) :: matrix
   integer(kind=8), allocatable, intent(out) :: row_offsets(:)
   integer(kind=4) :: i

   allocate(row_offsets(0:matrix%x))
   row_offsets(0) = 1_8

   do i = 0, matrix%x - 1
      row_offsets(i + 1) = row_offsets(i) + int(matrix%rows(i)%nnz, kind=8)
   end do
end subroutine build_mumps_row_offsets

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Fills MUMPS triplet arrays from independent row ranges.
!
!---------------------------------------------------------------------------
subroutine fill_mumps_triplets(matrix, row_offsets, transpose_triplets, mumps_par)
   implicit none
   include 'dmumps_struc.h'
   type(sparse_matrix), pointer, intent(in) :: matrix
   integer(kind=8), intent(in) :: row_offsets(0:matrix%x)
   logical, intent(in) :: transpose_triplets
   type(dmumps_struc), intent(inout) :: mumps_par
!> @brief Minimum number of triplets for OpenMP row conversion.
   integer(kind=8), parameter :: MUMPS_FORMAT_OMP_THRESHOLD = 4096_8
   integer(kind=4) :: i
   integer(kind=4) :: k
   integer(kind=4), allocatable :: order(:)
   integer(kind=8) :: pos

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k,pos,order) SCHEDULE(STATIC) &
!$OMP IF(matrix%total_entries > MUMPS_FORMAT_OMP_THRESHOLD)
   do i = 0, matrix%x - 1
      if (matrix%rows(i)%nnz <= 0) cycle

      pos = row_offsets(i)
      call sorted_row_order(matrix%rows(i), order)
      do k = 1, matrix%rows(i)%nnz
         if (transpose_triplets) then
            mumps_par%irn(pos + k - 1_8) = matrix%rows(i)%col(order(k)) + 1
            mumps_par%jcn(pos + k - 1_8) = i + 1
         else
            mumps_par%irn(pos + k - 1_8) = i + 1
            mumps_par%jcn(pos + k - 1_8) = matrix%rows(i)%col(order(k)) + 1
         end if
         mumps_par%a(pos + k - 1_8) = matrix%rows(i)%val(order(k))
      end do
      deallocate(order)
   end do
!$OMP END PARALLEL DO
end subroutine fill_mumps_triplets

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to the one-based MUMPS triplet format.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
   type(sparse_matrix), pointer, intent(in) :: matrix
   type(dmumps_struc), intent(inout) :: mumps_par
   integer(kind=8), allocatable :: row_offsets(:)

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if

   mumps_par%NZ = matrix%total_entries
   call release_mumps_triplets(mumps_par)

   allocate(mumps_par%irn(mumps_par%NZ))
   allocate(mumps_par%jcn(mumps_par%NZ))
   allocate(mumps_par%a(mumps_par%NZ))

   call build_mumps_row_offsets(matrix, row_offsets)
   call fill_mumps_triplets(matrix, row_offsets, .false., mumps_par)
   deallocate(row_offsets)
end subroutine to_mumps_format

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to a transposed MUMPS triplet format.
!
!> @details
!> The row and column indices are swapped and shifted to the one-based indexing
!> required by MUMPS.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format_transposed(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
   type(sparse_matrix), pointer, intent(in) :: matrix
   type(dmumps_struc), intent(inout) :: mumps_par
   integer(kind=8), allocatable :: row_offsets(:)

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if

   mumps_par%NZ = matrix%total_entries
   call release_mumps_triplets(mumps_par)

   allocate(mumps_par%irn(mumps_par%NZ))
   allocate(mumps_par%jcn(mumps_par%NZ))
   allocate(mumps_par%a(mumps_par%NZ))

   call build_mumps_row_offsets(matrix, row_offsets)
   call fill_mumps_triplets(matrix, row_offsets, .true., mumps_par)
   deallocate(row_offsets)
end subroutine to_mumps_format_transposed

end module sparse
