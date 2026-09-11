module sparse
   implicit none

   type sparse_matrix
      integer(kind=4) :: x = 0
      integer(kind=4) :: y = 0
   end type sparse_matrix

contains

   subroutine initialize_sparse(x, y, matrix)
      integer(kind=4), intent(in) :: x, y
      type(sparse_matrix), pointer, intent(out) :: matrix

      allocate(matrix)
      matrix%x = x
      matrix%y = y
   end subroutine initialize_sparse


   subroutine clear_matrix(matrix)
      type(sparse_matrix), pointer, intent(inout) :: matrix

      if (associated(matrix)) deallocate(matrix)
      nullify(matrix)
   end subroutine clear_matrix

end module sparse
