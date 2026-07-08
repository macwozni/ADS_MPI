!------------------------------------------------------------------------------
!
! MODULE: sparse
!
! DESCRIPTION:
!> @file sparse.F90
!> @brief Module providing linked-list sparse-matrix storage and helper
!> conversion routines.
!>
!> @details
!> This module defines a lightweight sparse-matrix infrastructure based
!> on linked lists of:
!> - sparse entries stored within one row-like structure,
!> - sparse lines grouped by the first matrix index,
!> - a top-level sparse matrix object containing dimensions and the total
!>   number of stored nonzero entries.
!>
!> The provided functionality includes:
!> - sparse-matrix initialization through \ref initialize_sparse,
!> - lookup or insertion of row structures and entries through
!>   \ref find_line and \ref find,
!> - additive assembly through \ref add,
!> - recursive cleanup through \ref clear_matrix,
!>   \ref clear_line, and \ref clear_entry,
!> - conversion to dense storage through \ref to_dense_matrix,
!> - conversion to MUMPS triplet storage through
!>   \ref to_mumps_format and \ref to_mumps_format_transposed.
!>
!> Within the project, this module is used primarily by the projection
!> and solver layers to assemble one-dimensional operators and pass them
!> to direct linear solvers.
!
!------------------------------------------------------------------------------
module sparse

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Linked-list node representing one sparse matrix entry.
!>
!> @details
!> This derived type stores:
!> - the first matrix index \p x,
!> - the second matrix index \p y,
!> - the numerical value \p val,
!> - a pointer to the next entry in the same sparse line.
!>
!> Entries inside one line are maintained in ascending order with
!> respect to the second index.
!
!---------------------------------------------------------------------------
   type sparse_matrix_entry
!> @brief First matrix index of the stored entry.
      integer(kind=4) :: x
!> @brief Second matrix index of the stored entry.
      integer(kind=4) :: y
!> @brief Numerical value of the stored entry.
      real(kind=8) :: val
!> @brief Pointer to the next entry in the same sparse line.
      type(sparse_matrix_entry), pointer :: next
   end type sparse_matrix_entry

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Linked-list node representing one sparse matrix line.
!>
!> @details
!> This derived type groups all sparse entries associated with one value
!> of the first matrix index. It stores:
!> - the line index \p x,
!> - a pointer to the next sparse line,
!> - a pointer to the first entry belonging to the current line.
!>
!> Sparse lines are maintained in ascending order with respect to the
!> line index.
!
!---------------------------------------------------------------------------
   type sparse_matrix_line
!> @brief First matrix index associated with the current sparse line.
      integer(kind=4) :: x
!> @brief Pointer to the next sparse line.
      type(sparse_matrix_line), pointer :: next
!> @brief Pointer to the first sparse entry stored in the current line.
      type(sparse_matrix_entry), pointer :: first
   end type sparse_matrix_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Top-level sparse matrix container.
!>
!> @details
!> This derived type stores:
!> - the nominal matrix dimensions \p x and \p y,
!> - the number of explicitly stored entries \p total_entries,
!> - a pointer to the first sparse line.
!>
!> The actual storage is organized as an ordered linked list of sparse
!> lines, each of which contains an ordered linked list of sparse
!> entries.
!
!---------------------------------------------------------------------------
   type sparse_matrix
!> @brief Size of the matrix in the first index direction.
      integer(kind=4) :: x
!> @brief Size of the matrix in the second index direction.
      integer(kind=4) :: y
!> @brief Total number of explicitly stored sparse entries.
      integer(kind=8) :: total_entries
!> @brief Pointer to the first sparse line.
      type(sparse_matrix_line), pointer :: first
   end type sparse_matrix

private :: clear_entry,   find, find_line,  clear_line, to_dense_matrix
private :: to_mumps_format_transposed, sparse_matrix_entry, sparse_matrix_line
public  :: initialize_sparse, clear_matrix, add, to_mumps_format
public  :: sparse_matrix

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Allocates and initializes an empty sparse matrix.
!>
!> @details
!> This routine allocates the sparse matrix object, stores its nominal
!> dimensions, initializes the entry counter to zero, and sets the
!> pointer to the first sparse line to `NULL()`.
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

   allocate (matrix)
   matrix%x = x
   matrix%y = y
   matrix%total_entries = 0
   matrix%first => NULL()
end subroutine initialize_sparse

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Finds or creates a sparse entry at the specified matrix
!> coordinates.
!>
!> @details
!> This routine first locates the sparse line associated with the first
!> index \p x by calling \ref find_line. It then searches the linked list
!> of entries in that line for the second index \p y.
!>
!> If the entry already exists, the corresponding pointer is returned.
!> Otherwise, a new entry is inserted in sorted order with initial value
!> zero, and the global entry counter \p matrix%total_entries is
!> incremented.
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
!> Sparse matrix searched and possibly extended.
!
! Output:
! -------
!> @param[out] entr
!> Pointer to the located or newly inserted sparse entry.
!
! Notes:
! ------
!> @note
!> Entries inside one sparse line are kept ordered by increasing
!> second-index value.
!
!---------------------------------------------------------------------------
subroutine find(matrix, x, y, entr)
   implicit none
!> @brief Sparse matrix searched and possibly extended.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief First matrix index of the requested entry.
   integer(kind=4), intent(in) :: x
!> @brief Second matrix index of the requested entry.
   integer(kind=4), intent(in) :: y
!> @brief Pointer to the located or inserted entry.
   type(sparse_matrix_entry), pointer, intent(out) :: entr
!> @brief Sparse line corresponding to the first index.
   type(sparse_matrix_line), pointer :: line
!> @brief Current entry pointer during traversal.
   type(sparse_matrix_entry), pointer :: tmp
!> @brief Auxiliary pointer used during insertion.
   type(sparse_matrix_entry), pointer :: tmp2

   call find_line(matrix, x, line)

   if (.NOT. associated(line%first)) then
      allocate (line%first)
      line%first%x = x
      line%first%y = y
      line%first%next => NULL()
      line%first%val = 0.d0
      entr => line%first
      matrix%total_entries = matrix%total_entries + 1
      return
   end if

   if (line%first%y .GT. y) then
      tmp2 => line%first
      allocate (line%first)
      line%first%x = x
      line%first%y = y
      line%first%next => tmp2
      line%first%val = 0.d0
      entr => line%first
      matrix%total_entries = matrix%total_entries + 1
      return
   end if

   tmp => line%first

   do while (associated(tmp%next))
      if (tmp%y .EQ. y) then
         entr => tmp
         return
      end if
      if (tmp%next%y .GT. y) then
         tmp2 => tmp%next
         allocate (tmp%next)
         tmp%next%x = x
         tmp%next%y = y
         tmp%next%next => tmp2
         tmp%next%val = 0.d0
         entr => tmp%next
         matrix%total_entries = matrix%total_entries + 1
         return
      end if
      tmp => tmp%next
   end do

   if (tmp%y .EQ. y) then
      entr => tmp
      return
   end if

   allocate (tmp%next)
   tmp%next%x = x
   tmp%next%y = y
   tmp%next%next => NULL()
   tmp%next%val = 0.d0
   entr => tmp%next
   matrix%total_entries = matrix%total_entries + 1
end subroutine find

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Finds or creates the sparse line associated with a specified
!> first matrix index.
!>
!> @details
!> This routine searches the linked list of sparse lines for the line
!> whose first-index value equals \p x.
!>
!> If such a line already exists, its pointer is returned. Otherwise, a
!> new line is inserted in sorted order and initialized with an empty
!> entry list.
!
! Input:
! ------
!> @param[in] x
!> First matrix index of the requested sparse line.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix searched and possibly extended.
!
! Output:
! -------
!> @param[out] line
!> Pointer to the located or newly inserted sparse line.
!
! Notes:
! ------
!> @note
!> Sparse lines are kept ordered by increasing first-index value.
!
!---------------------------------------------------------------------------
subroutine find_line(matrix, x, line)
   implicit none
!> @brief Sparse matrix searched and possibly extended.
   type(sparse_matrix), intent(inout) :: matrix
!> @brief First matrix index of the requested sparse line.
   integer(kind=4), intent(in) :: x
!> @brief Pointer to the located or inserted sparse line.
   type(sparse_matrix_line), pointer, intent(out) :: line
!> @brief Current sparse line pointer during traversal.
   type(sparse_matrix_line), pointer :: tmp
!> @brief Auxiliary pointer used during insertion.
   type(sparse_matrix_line), pointer :: tmp2

   if (.NOT. associated(matrix%first)) then
      allocate (matrix%first)
      matrix%first%x = x
      matrix%first%next => NULL()
      matrix%first%first => NULL()
      line => matrix%first
      return
   end if

   if (matrix%first%x .GT. x) then
      tmp2 => matrix%first
      allocate (matrix%first)
      matrix%first%x = x
      matrix%first%next => tmp2
      matrix%first%first => NULL()
      line => matrix%first
      return
   end if

   tmp => matrix%first

   do while (associated(tmp%next))
      if (tmp%x .EQ. x) then
         line => tmp
         return
      end if
      if (tmp%next%x .GT. x) then
         tmp2 => tmp%next
         allocate (tmp%next)
         tmp%next%x = x
         tmp%next%next => tmp2
         tmp%next%first => NULL()
         line => tmp%next
         return
      end if
      tmp => tmp%next
   end do

   if (tmp%x .EQ. x) then
      line => tmp
      return
   end if

   allocate (tmp%next)
   tmp%next%x = x
   tmp%next%next => NULL()
   tmp%next%first => NULL()
   line => tmp%next
end subroutine find_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Adds a value to a sparse matrix entry.
!>
!> @details
!> This routine locates the sparse entry at coordinates \p (x,y) by
!> calling \ref find and then increments its stored value by \p val.
!> If the entry does not exist yet, it is created with initial value zero
!> before the addition is applied.
!
! Input:
! ------
!> @param[in] x
!> First matrix index of the updated entry.
!>
!> @param[in] y
!> Second matrix index of the updated entry.
!>
!> @param[in] val
!> Value added to the sparse entry.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix updated in place.
!
!---------------------------------------------------------------------------
subroutine add(matrix, x, y, val)
   implicit none
!> @brief Sparse matrix updated in place.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief First matrix index of the updated entry.
   integer(kind=4), intent(in) :: x
!> @brief Second matrix index of the updated entry.
   integer(kind=4), intent(in) :: y
!> @brief Value added to the sparse entry.
   real(kind=8), intent(in) :: val
!> @brief Pointer to the located sparse entry.
   type(sparse_matrix_entry), pointer :: entr

   call find(matrix, x, y, entr)
   entr%val = entr%val + val
end subroutine add

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Releases the memory associated with a sparse matrix.
!>
!> @details
!> This routine deallocates the top-level sparse matrix object and then
!> recursively releases all sparse lines and entries previously linked to
!> it by calling \ref clear_line.
!
! Input/Output:
! -------------
!> @param[inout] matrix
!> Sparse matrix to be destroyed.
!
! Notes:
! ------
!> @note
!> The routine returns immediately if the input pointer is not
!> associated.
!
!---------------------------------------------------------------------------
subroutine clear_matrix(matrix)
   implicit none
!> @brief Sparse matrix to be destroyed.
   type(sparse_matrix), pointer, intent(inout) :: matrix
!> @brief Temporary pointer to the first sparse line.
   type(sparse_matrix_line), pointer :: line

   if (.NOT. associated(matrix)) return
   line => matrix%first
   deallocate (matrix)
   call clear_line(line)
end subroutine clear_matrix

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Recursively releases a linked list of sparse lines.
!>
!> @details
!> This routine deallocates the current sparse line, recursively clears
!> all sparse entries stored in that line through \ref clear_entry, and
!> then continues with the next sparse line.
!
! Input/Output:
! -------------
!> @param[inout] line
!> Pointer to the first sparse line of the list to be destroyed.
!
!---------------------------------------------------------------------------
recursive subroutine clear_line(line)
   implicit none
!> @brief Pointer to the sparse line to be destroyed.
   type(sparse_matrix_line), pointer, intent(inout) :: line
!> @brief Pointer to the next sparse line.
   type(sparse_matrix_line), pointer :: tmp

   if (.NOT. associated(line)) return
   tmp => line%next
   call clear_entry(line%first)
   deallocate (line)
   call clear_line(tmp)
end subroutine clear_line

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Recursively releases a linked list of sparse entries.
!>
!> @details
!> This routine deallocates the current sparse entry and then continues
!> recursively with the next entry in the same sparse line.
!
! Input/Output:
! -------------
!> @param[inout] entr
!> Pointer to the first sparse entry of the list to be destroyed.
!
!---------------------------------------------------------------------------
recursive subroutine clear_entry(entr)
   implicit none
!> @brief Pointer to the sparse entry to be destroyed.
   type(sparse_matrix_entry), pointer, intent(inout) :: entr
!> @brief Pointer to the next sparse entry.
   type(sparse_matrix_entry), pointer :: tmp

   if (.NOT. associated(entr)) return
   tmp => entr%next
   deallocate (entr)
   call clear_entry(tmp)
end subroutine clear_entry

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to dense two-dimensional storage.
!>
!> @details
!> This routine initializes the dense matrix \p dmatrix with zeros and
!> then traverses all sparse lines and entries stored in \p matrix. Each
!> stored sparse value is written to the corresponding position in the
!> dense output array.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted.
!
! Output:
! -------
!> @param[out] dmatrix
!> Dense matrix representation of the sparse input matrix.
!
! Notes:
! ------
!> @note
!> The dense matrix uses zero-based indexing in both dimensions in order
!> to match the indexing convention of the sparse storage.
!
!---------------------------------------------------------------------------
subroutine to_dense_matrix(matrix, dmatrix)
   implicit none
!> @brief Sparse matrix to be converted.
   type(sparse_matrix), intent(in) :: matrix
!> @brief Unused local integer variables retained in the original implementation.
   integer :: x, y
!> @brief Dense output matrix.
   real(kind=8), dimension(0:matrix%x - 1, 0:matrix%y - 1), intent(out) :: dmatrix
!> @brief Current sparse line during traversal.
   type(sparse_matrix_line), pointer :: line
!> @brief Current sparse entry during traversal.
   type(sparse_matrix_entry), pointer :: entr

   dmatrix = 0.d0
   line => matrix%first

   do while (associated(line))
      entr => line%first
      do while (associated(entr))
         dmatrix(entr%x, entr%y) = entr%val
         entr => entr%next
      end do
      line => line%next
   end do
end subroutine to_dense_matrix

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to MUMPS coordinate-triplet storage.
!>
!> @details
!> This routine traverses all sparse entries and writes them into the
!> arrays expected by the MUMPS direct solver:
!> - `irn` for row indices,
!> - `jcn` for column indices,
!> - `a`   for numerical values.
!>
!> The matrix order `N` stored in \p mumps_par is chosen as the larger of
!> the two nominal matrix dimensions. The number of explicitly stored
!> entries is copied to `NZ`.
!>
!> The row and column indices are shifted by one in order to comply with
!> the one-based indexing convention expected by MUMPS.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted.
!
! Input/Output:
! -------------
!> @param[inout] mumps_par
!> MUMPS structure filled with triplet storage arrays.
!
! Notes:
! ------
!> @note
!> The procedure allocates the arrays `irn`, `jcn`, and `a` inside
!> \p mumps_par.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
!> @brief Sparse matrix to be converted.
   type(sparse_matrix), pointer, intent(in) :: matrix
!> @brief MUMPS structure filled with triplet storage arrays.
   type(dmumps_struc), intent(inout) :: mumps_par
!> @brief Current sparse entry during traversal.
   type(sparse_matrix_entry), pointer :: entr
!> @brief Current sparse line during traversal.
   type(sparse_matrix_line), pointer :: line
!> @brief Running triplet-entry index.
   integer(kind=4) :: i

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if
   mumps_par%NZ = matrix%total_entries
   allocate (mumps_par%irn(matrix%total_entries))
   allocate (mumps_par%jcn(matrix%total_entries))
   allocate (mumps_par%a(matrix%total_entries))

   i = 1
   line => matrix%first
   do while (associated(line))
      entr => line%first
      do while (associated(entr))
         mumps_par%irn(i) = entr%x + 1
         mumps_par%jcn(i) = entr%y + 1
         mumps_par%a(i) = entr%val
         i = i + 1
         entr => entr%next
      end do
      line => line%next
   end do
end subroutine to_mumps_format

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts a sparse matrix to transposed MUMPS coordinate-triplet
!> storage.
!>
!> @details
!> This routine is analogous to \ref to_mumps_format, but it writes the
!> sparse entries with swapped matrix indices so that the transposed
!> operator is exported to MUMPS triplet storage.
!>
!> The matrix order `N` stored in \p mumps_par is chosen as the larger of
!> the two nominal matrix dimensions. The number of explicitly stored
!> entries is copied to `NZ`.
!
! Input:
! ------
!> @param[in] matrix
!> Sparse matrix to be converted in transposed form.
!
! Input/Output:
! -------------
!> @param[inout] mumps_par
!> MUMPS structure filled with transposed triplet storage arrays.
!
! Notes:
! ------
!> @note
!> The procedure allocates the arrays `irn`, `jcn`, and `a` inside
!> \p mumps_par.
!
!> @warning
!> The transposed index assignment follows the original implementation
!> verbatim.
!
!---------------------------------------------------------------------------
subroutine to_mumps_format_transposed(matrix, mumps_par)
   implicit none
   include 'dmumps_struc.h'
!> @brief Sparse matrix to be converted in transposed form.
   type(sparse_matrix), pointer, intent(in) :: matrix
!> @brief MUMPS structure filled with transposed triplet storage arrays.
   type(dmumps_struc), intent(inout) :: mumps_par
!> @brief Current sparse entry during traversal.
   type(sparse_matrix_entry), pointer :: entr
!> @brief Current sparse line during traversal.
   type(sparse_matrix_line), pointer :: line
!> @brief Running triplet-entry index.
   integer(kind=4) :: i

   if (matrix%x > matrix%y) then
      mumps_par%N = matrix%x
   else
      mumps_par%N = matrix%y
   end if
   mumps_par%NZ = matrix%total_entries
   allocate (mumps_par%irn(matrix%total_entries))
   allocate (mumps_par%jcn(matrix%total_entries))
   allocate (mumps_par%a(matrix%total_entries))

   i = 1
   line => matrix%first
   do while (associated(line))
      entr => line%first
      do while (associated(entr))
         mumps_par%irn(i) = entr%y
         mumps_par%jcn(i) = entr%x
         mumps_par%a(i) = entr%val
         i = i + 1
         entr => entr%next
      end do
      line => line%next
   end do
end subroutine to_mumps_format_transposed

end module sparse