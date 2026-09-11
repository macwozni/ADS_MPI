module projection_engine
   use sparse, only: sparse_matrix, initialize_sparse
   implicit none

   integer(kind=4) :: compute_calls = 0

contains

   subroutine reset_projection_spy()
      compute_calls = 0
   end subroutine reset_projection_spy


   subroutine ComputeMatrix(U1, p1, n1, nelem1, U2, p2, n2, nelem2, &
                            mixA, mixB, mixBT, equ, sprsmtrx)
      integer(kind=4), intent(in) :: p1, n1, nelem1, p2, n2, nelem2
      real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
      real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
      real(kind=8), intent(in) :: mixA(4), mixB(4), mixBT(4)
      logical, intent(in) :: equ
      type(sparse_matrix), pointer, intent(out) :: sprsmtrx

      compute_calls = compute_calls + 1
      call initialize_sparse(n1 + 1, n1 + 1, sprsmtrx)
   end subroutine ComputeMatrix

end module projection_engine
