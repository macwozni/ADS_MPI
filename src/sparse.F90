!------------------------------------------------------------------------------
!
! MODULE: sparse
!
! DESCRIPTION:
!> @file sparse.F90
!> @brief Module providing sparse-matrix storage and helper conversion
!> routines.
!>
!> @details
!> This module defines a sparse-matrix infrastructure optimized for
!> additive assembly. In contrast to the original linked-list-based
!> implementation, the current version stores the matrix row by row in
!> dynamically resizable arrays:
!> - one sparse row object per first-index value,
!> - sorted column indices within each sparse row,
!> - numerical values associated with the stored column indices.
!>
!> The public contract intentionally remains compatible with the rest of
!> the project:
!> - initialization through \ref initialize_sparse,
!> - additive assembly through \ref add,
!> - cleanup through \ref clear_matrix,
!> - conversion to dense storage through \ref to_dense_matrix,
!> - conversion to MUMPS triplet storage through
!>   \ref to_mumps_format and \ref to_mumps_format_transposed.
!>
!> The helper routines \ref find, \ref find_line, \ref clear_entry, and
!> \ref clear_line are kept private for compatibility with the historical
!> module structure, but they are no longer part of the supported public
!> API.
!
!------------------------------------------------------------------------------
module sparse

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Historical sparse-entry node type.
!>
!> @details
!> This type is retained for compatibility with the historical internal
!> API. The optimized sparse backend no longer uses linked lists as the
!> primary storage representation.
!
!---------------------------------------------------------------------------
   type sparse_matrix_entry
!> @brief First matrix index of the stored entry.
      integer(kind=4) :: x = 0
!> @brief Second matrix index of the stored entry.
      integer(kind=4) :: y = 0
!> @brief Numerical value of the stored entry.
      real(kind=8) :: val = 0.0d0
!> @brief Pointer to the next entry in the same sparse line.
      type(sparse_matrix_entry), pointer :: next => NULL()
   end type sparse_matrix_entry

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Historical sparse-line node type.
!>
!> @details
!> This type is retained for compatibility with the historical internal
!> API. The optimized sparse backend no longer uses linked lists of lines
!> as the primary storage representation.
!
!---------------------------------------------------------------------------
   type sparse_matrix_line
!> @brief First matrix index associated with the sparse line.
      integer(kind=4) :: x = 0
!> @brief Pointer to the next sparse line.
      type(sparse_matrix_line), pointer :: next => NULL()
!> @brief Pointer to the first sparse entry stored in the sparse line.
      type(sparse_matrix_entry), pointer :: first => NULL()
   end type sparse_matrix_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Internal row-wise sparse storage.
!>
!> @details
!> Each row stores:
!> - the number of explicitly stored entries,
!> - the currently allocated capacity,
!> - the sorted list of column indices,
!> - the values associated with the stored column indices.
!>
!> Entries inside a row are always kept sorted by increasing column
!> index.
!
!---------------------------------------------------------------------------
   type sparse_row_storage
!> @brief Number of explicitly stored entries in the row.
      integer(kind=4) :: nnz = 0
!> @brief Allocated row capacity.
      integer(kind=4) :: cap = 0
!> @brief Stored column indices.
      integer(kind=4), allocatable :: col(:)
!> @brief Stored numerical values.
      real(kind=8), allocatable :: val(:)
   end type sparse_row_storage

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Top-level sparse matrix container.
!>
!> @details
!> This derived type stores:
!> - the nominal matrix dimensions \p x and \p y,
!> - the total number of explicitly stored unique entries,
!> - row-wise sparse storage.
!>
!> The matrix is assembled additively. Repeated insertion into the same
!> matrix position updates the existing entry value and does not increase
!> \p total_entries.
!
!---------------------------------------------------------------------------
   type sparse_matrix
!> @brief Size of the matrix in the first index direction.
      integer(kind=4) :: x = 0
!> @brief Size of the matrix in the second index direction.
      integer(kind=4) :: y = 0
!> @brief Total number of explicitly stored unique sparse entries.
      integer(kind=8) :: total_entries = 0_8
!> @brief Row-wise sparse storage.
      type(sparse_row_storage), allocatable :: rows(:)
!> @brief Cache of the most recently used row index.
      integer(kind=4) :: last_row = -1
   end type sparse_matrix

private :: clear_entry, find, find_line, clear_line
private :: sparse_matrix_entry, sparse_matrix_line
public  :: initialize_sparse, clear_matrix, add, to_mumps_format
public  :: to_dense_matrix, to_mumps_format_transposed
public  :: sparse_matrix

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates and initializes an empty sparse matrix.
!>
!> @details
!> This routine allocates the sparse matrix object, stores its nominal
!> dimensions, initializes the entry counter to zero, resets the row
!> cache, and allocates row storage.
!
! Input:
! ------
!> @param[in] x
!> Size of the matrix in the first index direction.
!>
!> @param[in] y
!> Size of the matrix in the second index direction.
!
! Output:
! -------
!> @param[out] matrix
!> Newly allocated empty sparse matrix.
!
!---------------------------------------------------------------------------
subroutine initialize_sparse(x, y, matrix)
   implicit none
!> @brief Size of the matrix in the first direction.
   integer(kind=4), intent(in) :: x
!> @brief Size of the matrix in the second direction.
   integer(kind=4), intent(in) :: y
!> @brief Newly allocated sparse matrix.
   type(sparse_matrix), pointer, intent(out) :: matrix

   allocate(matrix)
   matrix%x = x
   matrix%y = y
   matrix%total_entries = 0_8
   matrix%last_row = -1

   if (x > 0) then
      allocate(matrix%rows(0:x-1))
   else
      allocate(matrix%rows(1:0))
   end if
end subroutine initialize_sparse

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Ensures that a sparse row has at least the requested capacity.
!>
!> @details
!> The row capacity is grown geometrically to limit the total number of
!> reallocations during sparse assembly.
!
! Input:
! ------
!> @param[in] needed
!> Required minimum row capacity.
!
! Input/Output:
! -------------
!> @param[inout] row
!> Sparse row whose storage may be reallocated.
!
!---------------------------------------------------------------------------
subroutine ensure_row_capacity(row, needed)
   implicit none
!> @brief Sparse row whose storage may be reallocated.
   type(sparse_row_storage), intent(inout) :: row
!> @brief Required minimum row capacity.
   integer(kind=4), intent(in) :: needed
!> @brief New capacity after growth.
   integer(kind=4) :: newcap
!> @brief Reallocated column-index buffer.
   integer(kind=4), allocatable :: newcol(:)
!> @brief Reallocated value buffer.
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
!> @brief Locates a column index inside one sparse row.
!>
!> @details
!> This routine performs a binary search in the sorted list of column
!> indices stored in one row. It returns:
!> - whether the requested column index already exists,
!> - the existing position if found,
!> - the insertion position if not found.
!
! Input:
! ------
!> @param[in] y
!> Column index to be located.
!>
!> @param[in] row
!> Sparse row searched by the routine.
!
! Output:
! -------
!> @param[out] found
!> Logical flag indicating whether the requested column index was found.
!>
!> @param[out] pos
!> Existing position if found, otherwise insertion position preserving
!> sorted order.
!
!---------------------------------------------------------------------------
subroutine locate_in_row(row, y, found, pos)
   implicit none
!> @brief Sparse row searched by the routine.
   type(sparse_row_storage), intent(in) :: row
!> @brief Column index to be located.
   integer(kind=4), intent(in) :: y
!> @brief Logical flag indicating whether the column index was found.
   logical, intent(out) :: found
!> @brief Existing or insertion position.
   integer(kind=4), intent(out) :: pos
!> @brief Lower bound of the binary-search window.
   integer(kind=4) :: lo
!> @brief Upper bound of the binary-search window.
   integer(kind=4) :: hi
!> @brief Midpoint of the binary-search window.
   integer(kind=4) :: mid

   found = .false.

   if (row%nnz == 0) then
      pos = 1
      return
   end if

   if (y < row%col(1)) then
      pos = 1
      return
   end if

   if (y > row%col(row%nnz)) then
      pos = row%nnz + 1
      return
   end if

   if (y == row%col(row%nnz)) then
      found = .true.
      pos = row%nnz
      return
   end if

   lo = 1
   hi = row%nnz

   do while (lo <= hi)
      mid = (lo + hi) / 2

      if (row%col(mid) == y) then
         found = .true.
         pos = mid
         return
      else if (row%col(mid) < y) then
         lo = mid + 1
      else
         hi = mid - 1
      end if
   end do

   pos = lo
end subroutine locate_in_row

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Inserts a new entry into one sparse row.
!>
!> @details
!> This routine inserts a new matrix entry into the specified row at the
!> requested position. The sorted column order is preserved.
!
! Input:
! ------
!> @param[in] pos
!> Insertion position preserving sorted order.
!>
!> @param[in] y
!> Column index of the inserted entry.
!>
!> @param[in] value
!> Numerical value of the inserted entry.
!
! Input/Output:
! -------------
!> @param[inout] row
!> Sparse row modified by the insertion.
!
!---------------------------------------------------------------------------
subroutine insert_into_row(row, pos, y, value)
   implicit none
!> @brief Sparse row modified by the insertion.
   type(sparse_row_storage), intent(inout) :: row
!> @brief Insertion position preserving sorted order.
   integer(kind=4), intent(in) :: pos
!> @brief Column index of the inserted entry.
   integer(kind=4), intent(in) :: y
!> @brief Numerical value of the inserted entry.
   real(kind=8), intent(in) :: value

   call ensure_row_capacity(row, row%nnz + 1)

   if (pos <= row%nnz) then
      row%col(pos+1:row%nnz+1) = row%col(pos:row%nnz)
      row%val(pos+1:row%nnz+1) = row%val(pos:row%nnz)
   end if

   row%col(pos) = y
   row%val(pos) = value
   row%nnz = row%nnz + 1
end subroutine insert_into_row

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Adds a numerical contribution to a sparse matrix entry.
!>
!> @details
!> If the entry \f$(x,y)\f$ is already present, the value is accumulated
!> into the existing storage location. Otherwise, a new sparse entry is
!> inserted and the global unique-entry counter is incremented.
!>
!> Entries within one row are kept sorted by increasing second index.
!>
!> The historical contract of the module is preserved: this routine does
!> not perform defensive range validation of \p x and \p y.
!
! Input:
! ------
!> @param[in] x
!> First matrix index of the modified entry.
!>
!> @param[in] y
!> Second matrix index of the modified entry.
!>
!> @param[in] val
!> Numerical contribution to be added.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix modified by the additive update.
!
!---------------------------------------------------------------------------
subroutine add(matrix, x, y, val)
   implicit none
!> @brief Sparse matrix modified by the additive update.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief First matrix index of the modified entry.
   integer(kind=4), intent(in) :: x
!> @brief Second matrix index of the modified entry.
   integer(kind=4), intent(in) :: y
!> @brief Numerical contribution to be added.
   real(kind=8), intent(in) :: val
!> @brief Indicates whether the requested column index already exists.
   logical :: found
!> @brief Existing or insertion position inside the sparse row.
   integer(kind=4) :: pos

   call locate_in_row(matrix%rows(x), y, found, pos)

   if (found) then
      matrix%rows(x)%val(pos) = matrix%rows(x)%val(pos) + val
   else
      call insert_into_row(matrix%rows(x), pos, y, val)
      matrix%total_entries = matrix%total_entries + 1_8
   end if

   matrix%last_row = x
end subroutine add

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Finds a sparse entry and returns its current value.
!>
!> @details
!> This private helper routine is retained for historical compatibility.
!> It returns a newly allocated temporary entry object containing the
!> current value stored at \f$(x,y)\f$. If the entry does not exist, the
!> returned value is zero.
!
! Input:
! ------
!> @param[in] x
!> First matrix index of the requested entry.
!>
!> @param[in] y
!> Second matrix index of the requested entry.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix searched by the routine.
!
! Output:
! -------
!> @param[out] entr
!> Newly allocated temporary entry containing the requested value.
!
!---------------------------------------------------------------------------
subroutine find(matrix, x, y, entr)
   implicit none
!> @brief Sparse matrix searched by the routine.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief First matrix index of the requested entry.
   integer(kind=4), intent(in) :: x
!> @brief Second matrix index of the requested entry.
   integer(kind=4), intent(in) :: y
!> @brief Newly allocated temporary entry containing the requested value.
   type(sparse_matrix_entry), pointer, intent(out) :: entr
!> @brief Indicates whether the requested column index exists.
   logical :: found
!> @brief Existing or insertion position inside the sparse row.
   integer(kind=4) :: pos

   allocate(entr)
   entr%x = x
   entr%y = y
   entr%next => NULL()

   call locate_in_row(matrix%rows(x), y, found, pos)

   if (found) then
      entr%val = matrix%rows(x)%val(pos)
   else
      entr%val = 0.0d0
   end if
end subroutine find

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Historical helper returning a temporary sparse-line object.
!>
!> @details
!> This private helper is retained only for compatibility with the former
!> module structure. The optimized backend does not use linked lists of
!> sparse lines internally.
!
! Input:
! ------
!> @param[in] x
!> First matrix index associated with the returned temporary line.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix referenced only for compatibility.
!
! Output:
! -------
!> @param[out] line
!> Newly allocated temporary sparse-line object.
!
!---------------------------------------------------------------------------
subroutine find_line(matrix, x, line)
   implicit none
!> @brief Sparse matrix referenced only for compatibility.
   type(sparse_matrix), intent(inout) :: matrix
!> @brief First matrix index associated with the returned temporary line.
   integer(kind=4), intent(in) :: x
!> @brief Newly allocated temporary sparse-line object.
   type(sparse_matrix_line), pointer, intent(out) :: line

   allocate(line)
   line%x = x
   line%next => NULL()
   line%first => NULL()

   if (matrix%last_row == -huge(1)) then
      line%x = x
   end if
end subroutine find_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Clears a temporary sparse-entry object.
!>
!> @details
!> This private helper deallocates the temporary object returned by
!> \ref find. It is retained for compatibility with the former module
!> structure.
!
! Input/Output:
! -------------
!> @param[inout] entr
!> Temporary sparse-entry object to be deallocated.
!
!---------------------------------------------------------------------------
subroutine clear_entry(entr)
   implicit none
!> @brief Temporary sparse-entry object to be deallocated.
   type(sparse_matrix_entry), pointer, intent(inout) :: entr

   if (.not. associated(entr)) return
   deallocate(entr)
   nullify(entr)
end subroutine clear_entry

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Clears a temporary sparse-line object.
!>
!> @details
!> This private helper deallocates the temporary object returned by
!> \ref find_line. It is retained for compatibility with the former
!> module structure.
!
! Input/Output:
! -------------
!> @param[inout] line
!> Temporary sparse-line object to be deallocated.
!
!---------------------------------------------------------------------------
subroutine clear_line(line)
   implicit none
!> @brief Temporary sparse-line object to be deallocated.
   type(sparse_matrix_line), pointer, intent(inout) :: line

   if (.not. associated(line)) return
   deallocate(line)
   nullify(line)
end subroutine clear_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Releases all storage owned by a sparse matrix.
!>
!> @details
!> This routine deallocates all row buffers, deallocates the row-storage
!> array, deallocates the sparse matrix object itself, and nullifies the
!> caller's pointer.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix to be deallocated.
!
!---------------------------------------------------------------------------
subroutine clear_matrix(matrix)
   implicit none
!> @brief Sparse matrix to be deallocated.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief Row-loop index.
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
!>
!> @details
!> This routine fills the dense matrix with zeros and then copies all
!> explicitly stored sparse entries into their corresponding dense
!> positions. Because the sparse backend stores unique entries only,
!> direct assignment is sufficient.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted.
!
! Output:
! -------
!> @param[out] dmatrix
!> Dense matrix representation of the sparse matrix.
!
!---------------------------------------------------------------------------
subroutine to_dense_matrix(matrix, dmatrix)
   implicit none
!> @brief Sparse matrix to be converted.
   type(sparse_matrix), intent(in) :: matrix
!> @brief Dense matrix representation of the sparse matrix.
   real(kind=8), dimension(0:matrix%x - 1, 0:matrix%y - 1), intent(out) :: dmatrix
!> @brief Row-loop index.
   integer(kind=4) :: i
!> @brief Entry-loop index inside one sparse row.
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
!> @brief Converts a sparse matrix to the MUMPS triplet format.
!>
!> @details
!> The routine fills the MUMPS arrays:
!> - \p irn with one-based row indices,
!> - \p jcn with one-based column indices,
!> - \p a   with numerical values.
!>
!> The entries are exported row by row and, within each row, in
!> increasing column order.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted.
!
! Input/Output:
! -------------
!> @param[inout] mumps_par
!> MUMPS structure receiving the triplet representation.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
!> @brief Sparse matrix to be converted.
   type(sparse_matrix), pointer, intent(in) :: matrix
!> @brief MUMPS structure receiving the triplet representation.
   type(dmumps_struc), intent(inout) :: mumps_par
!> @brief Row-loop index.
   integer(kind=4) :: i
!> @brief Entry-loop index inside one sparse row.
   integer(kind=4) :: k
!> @brief Running counter of exported nonzero entries.
   integer(kind=8) :: nz_counter

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if

   mumps_par%NZ = matrix%total_entries

   if (associated(mumps_par%irn)) deallocate(mumps_par%irn)
   if (associated(mumps_par%jcn)) deallocate(mumps_par%jcn)
   if (associated(mumps_par%a))   deallocate(mumps_par%a)

   allocate(mumps_par%irn(mumps_par%NZ))
   allocate(mumps_par%jcn(mumps_par%NZ))
   allocate(mumps_par%a(mumps_par%NZ))

   nz_counter = 0_8
   do i = 0, matrix%x - 1
      do k = 1, matrix%rows(i)%nnz
         nz_counter = nz_counter + 1_8
         mumps_par%irn(nz_counter) = i + 1
         mumps_par%jcn(nz_counter) = matrix%rows(i)%col(k) + 1
         mumps_par%a(nz_counter)   = matrix%rows(i)%val(k)
      end do
   end do
end subroutine to_mumps_format

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to the transposed MUMPS triplet
!> format.
!>
!> @details
!> This routine preserves the historical project contract of
!> \ref to_mumps_format_transposed. In particular, it exports:
!> - \p irn as the stored column index,
!> - \p jcn as the stored row index,
!> - \p a   as the stored numerical value,
!>
!> without adding the one-based shift used by \ref to_mumps_format.
!> Although this behavior differs from the standard one-based MUMPS
!> convention, it matches the legacy implementation and the current test
!> contract of the code base.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted.
!
! Input/Output:
! -------------
!> @param[inout] mumps_par
!> MUMPS structure receiving the transposed triplet representation.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format_transposed(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
!> @brief Sparse matrix to be converted.
   type(sparse_matrix), pointer, intent(in) :: matrix
!> @brief MUMPS structure receiving the transposed triplet representation.
   type(dmumps_struc), intent(inout) :: mumps_par
!> @brief Row-loop index.
   integer(kind=4) :: i
!> @brief Entry-loop index inside one sparse row.
   integer(kind=4) :: k
!> @brief Running counter of exported nonzero entries.
   integer(kind=8) :: nz_counter

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if

   mumps_par%NZ = matrix%total_entries

   if (associated(mumps_par%irn)) deallocate(mumps_par%irn)
   if (associated(mumps_par%jcn)) deallocate(mumps_par%jcn)
   if (associated(mumps_par%a))   deallocate(mumps_par%a)

   allocate(mumps_par%irn(mumps_par%NZ))
   allocate(mumps_par%jcn(mumps_par%NZ))
   allocate(mumps_par%a(mumps_par%NZ))

   nz_counter = 0_8
   do i = 0, matrix%x - 1
      do k = 1, matrix%rows(i)%nnz
         nz_counter = nz_counter + 1_8
         mumps_par%irn(nz_counter) = matrix%rows(i)%col(k)
         mumps_par%jcn(nz_counter) = i
         mumps_par%a(nz_counter)   = matrix%rows(i)%val(k)
      end do
   end do
end subroutine to_mumps_format_transposed

end module sparse