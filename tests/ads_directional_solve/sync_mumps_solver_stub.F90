module mumps_solver
   use sparse, only: sparse_matrix
   implicit none

   integer(kind=4) :: solve_calls = 0
   integer(kind=4) :: injected_solver_status = 0

contains

   subroutine configure_solver_failure(status)
      integer(kind=4), intent(in) :: status

      solve_calls = 0
      injected_solver_status = status
   end subroutine configure_solver_failure


   subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx, ierr)
      real(kind=8), intent(inout) :: RHS(:, :)
      integer(kind=4), intent(in) :: eqnum, n, p
      type(sparse_matrix), pointer, intent(in) :: sprsmtrx
      integer(kind=4), intent(out) :: ierr

      solve_calls = solve_calls + 1
      ierr = injected_solver_status
      if (ierr == 0) RHS = 2.d0*RHS
   end subroutine SolveOneDirection

end module mumps_solver
