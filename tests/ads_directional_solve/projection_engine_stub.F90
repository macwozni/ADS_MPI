module projection_engine
   use sparse, only: sparse_matrix, initialize_sparse
   use directional_test_support, only: compute_calls, compute_contract_ok, &
      expected_p_test, expected_n_test, expected_nelem_test, &
      expected_p_trial, expected_n_trial, expected_nelem_trial, &
      expected_matrix_size, expected_equ, expected_mixA, expected_mixB, &
      expected_mixBT, expected_u_test, expected_u_trial, &
      injected_matrix_size_delta
   implicit none

contains

   subroutine ComputeMatrix(U1, p1, n1, nelem1, U2, p2, n2, nelem2, &
                            mixA, mixB, mixBT, equ, sprsmtrx)
      integer(kind=4), intent(in) :: p1, n1, nelem1, p2, n2, nelem2
      real(kind=8), intent(in) :: U1(0:n1 + p1 + 1)
      real(kind=8), intent(in) :: U2(0:n2 + p2 + 1)
      real(kind=8), intent(in) :: mixA(4), mixB(4), mixBT(4)
      logical, intent(in) :: equ
      type(sparse_matrix), pointer, intent(out) :: sprsmtrx

      compute_calls = compute_calls + 1
      compute_contract_ok = compute_contract_ok .and. &
         p1 == expected_p_test .and. n1 == expected_n_test .and. &
         nelem1 == expected_nelem_test
      compute_contract_ok = compute_contract_ok .and. &
         p2 == expected_p_trial .and. n2 == expected_n_trial .and. &
         nelem2 == expected_nelem_trial
      compute_contract_ok = compute_contract_ok .and. &
         (equ .eqv. expected_equ)
      compute_contract_ok = compute_contract_ok .and. &
         all(mixA == expected_mixA) .and. all(mixB == expected_mixB) .and. &
         all(mixBT == expected_mixBT)

      if (size(U1) == size(expected_u_test)) then
         compute_contract_ok = compute_contract_ok .and. &
            all(U1 == expected_u_test)
      else
         compute_contract_ok = .false.
      end if
      if (size(U2) == size(expected_u_trial)) then
         compute_contract_ok = compute_contract_ok .and. &
            all(U2 == expected_u_trial)
      else
         compute_contract_ok = .false.
      end if

      call initialize_sparse(expected_matrix_size + injected_matrix_size_delta, &
                             expected_matrix_size + injected_matrix_size_delta, &
                             sprsmtrx)
   end subroutine ComputeMatrix

end module projection_engine
