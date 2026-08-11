program test_mumps_solver
   use sparse, only: sparse_matrix, initialize_sparse, clear_matrix, add
   use mumps_solver, only: SolveOneDirection
   use mumps_test_support, only: configure_mumps_stub, expected_jobs, contract_recorded, &
      recorded_comm, recorded_sym, recorded_par, recorded_n, recorded_nz, &
      recorded_icntl, recorded_irn, recorded_jcn, recorded_a
   use mpi, only: MPI_COMM_SELF
   implicit none

   type(sparse_matrix), pointer :: matrix
   integer(kind=4) :: checks, failures

   checks = 0
   failures = 0
   call initialize_sparse(2, 2, matrix)
   call add(matrix, 0, 0, 2.d0)
   call add(matrix, 1, 1, 3.d0)

   call test_success(matrix, failures)
   call test_legacy_success_interface(matrix, failures)
   call test_initialization_error(matrix, failures)
   call test_initialization_error_without_instance(matrix, failures)
   call test_analysis_error(matrix, failures)
   call test_factorization_error(matrix, failures)
   call test_solve_error(matrix, failures)
   call test_later_rhs_error(matrix, failures)
   call test_finalization_error(matrix, failures)
   call test_first_error_survives_finalization(matrix, failures)
   call test_zero_right_hand_sides(matrix, failures)
   call test_success(matrix, failures)

   call clear_matrix(matrix)

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' MUMPS solver error-path checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' MUMPS solver error-path checks)'
      stop 1
   end if

contains

   subroutine test_success(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2), expected(2, 2)
      integer(kind=4) :: status

      rhs = reshape((/1.d0, 2.d0, 3.d0, 4.d0/), shape(rhs))
      expected = 2.d0*rhs
      call configure_mumps_stub(999, 0)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('success returns zero', status == 0, failure_count)
      call assert_true('success copies every solved RHS', all(rhs == expected), &
                       failure_count)
      call assert_true('success executes every MUMPS phase', &
                       expected_jobs((/-1, 1, 2, 3, 3, -2/)), failure_count)
      call assert_true('analysis receives a recorded MUMPS contract', &
                       contract_recorded, failure_count)
      call assert_true('each rank solves through MPI_COMM_SELF', &
                       recorded_comm == MPI_COMM_SELF, failure_count)
      call assert_true('MUMPS receives unsymmetric host participation', &
                       recorded_sym == 0 .and. recorded_par == 1, failure_count)
      call assert_true('MUMPS receives exact matrix dimensions', &
                       recorded_n == 2 .and. recorded_nz == 2, failure_count)
      call assert_true('MUMPS receives exact one-based row triplets', &
                       all(recorded_irn(1:2) == (/1, 2/)), failure_count)
      call assert_true('MUMPS receives exact one-based column triplets', &
                       all(recorded_jcn(1:2) == (/1, 2/)), failure_count)
      call assert_true('MUMPS receives exact matrix values', &
                       all(recorded_a(1:2) == (/2.d0, 3.d0/)), failure_count)
      call assert_true('MUMPS receives required analysis controls', &
                       recorded_icntl(1) == 1 .and. recorded_icntl(2) == 0 .and. &
                       recorded_icntl(3) == 0 .and. recorded_icntl(4) == 0 .and. &
                       recorded_icntl(5) == 0 .and. recorded_icntl(7) == 7 .and. &
                       recorded_icntl(11) == 2 .and. recorded_icntl(12) == 0 .and. &
                       recorded_icntl(13) == 1 .and. recorded_icntl(14) == 50 .and. &
                       recorded_icntl(18) == 0 .and. recorded_icntl(19) == 0 .and. &
                       recorded_icntl(20) == 0, failure_count)
   end subroutine test_success


   subroutine test_legacy_success_interface(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 1), expected(2, 1)

      rhs(:, 1) = (/5.d0, 7.d0/)
      expected = 2.d0*rhs
      call configure_mumps_stub(999, 0)
      call SolveOneDirection(rhs, 1, 1, 1, sparse_matrix_ptr)

      call assert_true('legacy five-argument interface still solves', &
                       all(rhs == expected), failure_count)
      call assert_true('legacy success still finalizes MUMPS', &
                       expected_jobs((/-1, 1, 2, 3, -2/)), failure_count)
   end subroutine test_legacy_success_interface


   subroutine test_initialization_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2), before(2, 2)
      integer(kind=4) :: status

      rhs = 1.d0
      before = rhs
      call configure_mumps_stub(-1, -16)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('initialization error is returned', status == -16, &
                       failure_count)
      call assert_true('initialization error leaves RHS untouched', &
                       all(rhs == before), failure_count)
      call assert_true('failed initialization is finalized without later phases', &
                       expected_jobs((/-1, -2/)), failure_count)
   end subroutine test_initialization_error


   subroutine test_initialization_error_without_instance(sparse_matrix_ptr, &
                                                         failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2), before(2, 2)
      integer(kind=4) :: status

      rhs = reshape((/1.d0, 2.d0, 3.d0, 4.d0/), shape(rhs))
      before = rhs
      call configure_mumps_stub(-1, -23, start_instance=.false.)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('pre-instance initialization error is returned', &
                       status == -23, failure_count)
      call assert_true('pre-instance initialization error leaves RHS untouched', &
                       all(rhs == before), failure_count)
      call assert_true('an instance which never started is not finalized', &
                       expected_jobs((/-1/)), failure_count)
   end subroutine test_initialization_error_without_instance


   subroutine test_analysis_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2)
      integer(kind=4) :: status

      rhs = 1.d0
      call configure_mumps_stub(1, -17)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('analysis error is returned', status == -17, failure_count)
      call assert_true('analysis error still finalizes MUMPS', &
                       expected_jobs((/-1, 1, -2/)), failure_count)
   end subroutine test_analysis_error


   subroutine test_factorization_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2)
      integer(kind=4) :: status

      rhs = 1.d0
      call configure_mumps_stub(2, -18, use_infog=.true.)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('global factorization error is returned', status == -18, &
                       failure_count)
      call assert_true('factorization error skips solves and finalizes', &
                       expected_jobs((/-1, 1, 2, -2/)), failure_count)
   end subroutine test_factorization_error


   subroutine test_solve_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2), before(2, 2)
      integer(kind=4) :: status

      rhs = reshape((/1.d0, 2.d0, 3.d0, 4.d0/), shape(rhs))
      before = rhs
      call configure_mumps_stub(3, -19)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('solve error is returned', status == -19, failure_count)
      call assert_true('failed solve is not copied into RHS', all(rhs == before), &
                       failure_count)
      call assert_true('solve error skips remaining RHS and finalizes', &
                       expected_jobs((/-1, 1, 2, 3, -2/)), failure_count)
   end subroutine test_solve_error


   subroutine test_later_rhs_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2), before(2, 2)
      integer(kind=4) :: status

      rhs = reshape((/1.d0, 2.d0, 3.d0, 4.d0/), shape(rhs))
      before = rhs
      call configure_mumps_stub(3, -22, occurrence=2)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('later RHS solve error is returned', status == -22, &
                       failure_count)
      call assert_true('completed RHS remains solved after a later error', &
                       all(rhs(:, 1) == 2.d0*before(:, 1)), failure_count)
      call assert_true('failed later RHS is not copied back', &
                       all(rhs(:, 2) == before(:, 2)), failure_count)
      call assert_true('later RHS error still finalizes once', &
                       expected_jobs((/-1, 1, 2, 3, 3, -2/)), failure_count)
   end subroutine test_later_rhs_error


   subroutine test_finalization_error(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2)
      integer(kind=4) :: status

      rhs = 1.d0
      call configure_mumps_stub(999, 0, finalize_code=-20)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('finalization error is returned', status == -20, &
                       failure_count)
      call assert_true('finalization error occurs after all solves', &
                       expected_jobs((/-1, 1, 2, 3, 3, -2/)), failure_count)
   end subroutine test_finalization_error


   subroutine test_first_error_survives_finalization(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 2)
      integer(kind=4) :: status

      rhs = 1.d0
      call configure_mumps_stub(2, -21, finalize_code=-99)
      call SolveOneDirection(rhs, 2, 1, 1, sparse_matrix_ptr, status)

      call assert_true('finalization cannot overwrite the first error', &
                       status == -21, failure_count)
      call assert_true('first-error path finalizes exactly once', &
                       expected_jobs((/-1, 1, 2, -2/)), failure_count)
   end subroutine test_first_error_survives_finalization


   subroutine test_zero_right_hand_sides(sparse_matrix_ptr, failure_count)
      type(sparse_matrix), pointer, intent(in) :: sparse_matrix_ptr
      integer(kind=4), intent(inout) :: failure_count
      real(kind=8) :: rhs(2, 0)
      integer(kind=4) :: status

      call configure_mumps_stub(999, 0)
      call SolveOneDirection(rhs, 0, 1, 1, sparse_matrix_ptr, status)

      call assert_true('zero RHS columns return success', status == 0, &
                       failure_count)
      call assert_true('zero RHS columns skip solve jobs but finalize', &
                       expected_jobs((/-1, 1, 2, -2/)), failure_count)
   end subroutine test_zero_right_hand_sides


   subroutine assert_true(label, condition, failure_count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      integer(kind=4), intent(inout) :: failure_count

      checks = checks + 1
      if (.not. condition) then
         failure_count = failure_count + 1
         write (*, '(A,A)') 'FAIL: ', trim(label)
      end if
   end subroutine assert_true

end program test_mumps_solver
