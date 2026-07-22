module mumps_solver
   use sparse, only: sparse_matrix
   use directional_test_support, only: solve_calls, solve_contract_ok, &
      expected_matrix_size, expected_packed_rhs, exact_matrix
   implicit none

contains

   subroutine SolveOneDirection(RHS, eqnum, n, p, sprsmtrx)
      real(kind=8), intent(inout) :: RHS(:, :)
      integer(kind=4), intent(in) :: eqnum, n, p
      type(sparse_matrix), pointer, intent(in) :: sprsmtrx

      solve_calls = solve_calls + 1
      solve_contract_ok = solve_contract_ok .and. &
         eqnum == size(expected_packed_rhs, 2)
      solve_contract_ok = solve_contract_ok .and. &
         n == expected_matrix_size - 1 .and. p == expected_matrix_size - 1
      solve_contract_ok = solve_contract_ok .and. &
         exact_matrix(RHS, expected_packed_rhs)

      if (associated(sprsmtrx)) then
         solve_contract_ok = solve_contract_ok .and. &
            sprsmtrx%x == expected_matrix_size .and. &
            sprsmtrx%y == expected_matrix_size
      else
         solve_contract_ok = .false.
      end if

      RHS = 2.d0*RHS
   end subroutine SolveOneDirection

end module mumps_solver
