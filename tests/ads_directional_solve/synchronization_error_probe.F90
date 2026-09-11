program synchronization_error_probe
   use Setup, only: ADS_Setup
   use mpi, only: MPI_COMM_WORLD, FAIL_ON_LOCAL_ERROR, allreduce_calls, &
                  bcast_calls, abort_calls, &
                  abort_code, abort_comm, bcast_root, collective_contract_ok, &
                  configure_collective_failures
   use my_mpi, only: configure_transport_failure, transport_operation_count, &
                     gather_call_count, scatter_call_count
   use projection_engine, only: compute_calls, reset_projection_spy
   use mumps_solver, only: solve_calls, configure_solver_failure
   use reorderRHS, only: reorder_calls, reset_reorder_spy
   use ads_directional_solve, only: solve_problem
   implicit none

   integer(kind=4), parameter :: ALLREDUCE_INITIAL_ERROR = 7201
   integer(kind=4), parameter :: ALLREDUCE_AFTER_SOLVER_ERROR = 7202
   integer(kind=4), parameter :: BCAST_ERROR = 7203
   integer(kind=4), parameter :: SOLVER_ERROR_FOR_REDUCTION = 7301
   integer(kind=4), parameter :: SOLVER_ERROR_FOR_BCAST = 7302
   real(kind=8), parameter :: OUTPUT_SENTINEL = -987654321.25d0
   type(ADS_Setup) :: ads_test, ads_trial
   integer(kind=4) :: checks, failures
   character(len=32) :: test_case

   checks = 0
   failures = 0
   call prepare_space(ads_trial)
   call prepare_space(ads_test)
   call get_command_argument(1, test_case)

   select case (trim(test_case))
   case ('initial-allreduce')
      call test_initial_allreduce_failure()
   case ('solver-allreduce')
      call test_allreduce_failure_preserves_solver_error()
   case ('solver-bcast')
      call test_bcast_failure_preserves_solver_error()
   case ('')
      call test_initial_allreduce_failure()
      call test_allreduce_failure_preserves_solver_error()
      call test_bcast_failure_preserves_solver_error()
   case default
      write (*, '(A)') 'usage: synchronization_error_probe '// &
                       '{initial-allreduce|solver-allreduce|solver-bcast}'
      stop 2
   end select

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' synchronization-error checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' synchronization-error checks)'
      stop 1
   end if

contains

   subroutine test_initial_allreduce_failure()
      real(kind=8), allocatable :: F(:, :), F2(:, :), Ft(:, :), Ft2(:, :)
      integer(kind=4) :: status

      call prepare_case(F, F2)
      call configure_collective_failures(1, ALLREDUCE_INITIAL_ERROR, 0)
      call configure_transport_failure(0, 0)
      call configure_solver_failure(0)
      call reset_projection_spy()
      call reset_reorder_spy()

      call run_solve(F, F2, Ft, Ft2, status)

      call assert_true('Allreduce failure returns its exact MPI status', &
                       status == ALLREDUCE_INITIAL_ERROR)
      call assert_true('Allreduce failure aborts once with its MPI status', &
                       abort_calls == 1 .and. &
                       abort_code == ALLREDUCE_INITIAL_ERROR .and. &
                       abort_comm == MPI_COMM_WORLD)
      call assert_true('Allreduce failure never enters Bcast or later workflow', &
                       allreduce_calls == 1 .and. bcast_calls == 0 .and. &
                       transport_operation_count == 0 .and. &
                       compute_calls == 0 .and. solve_calls == 0 .and. &
                       scatter_call_count == 0 .and. reorder_calls == 0)
      call assert_true('Allreduce failure leaves caller buffers unchanged', &
                       F(1, 1) == 41.d0 .and. F2(1, 1) == OUTPUT_SENTINEL)
      call assert_true('Allreduce failure uses the expected collective contract', &
                       collective_contract_ok)
   end subroutine test_initial_allreduce_failure


   subroutine test_allreduce_failure_preserves_solver_error()
      real(kind=8), allocatable :: F(:, :), F2(:, :), Ft(:, :), Ft2(:, :)
      integer(kind=4) :: status

      call prepare_case(F, F2)
      call configure_collective_failures(FAIL_ON_LOCAL_ERROR, &
                                         ALLREDUCE_AFTER_SOLVER_ERROR, 0)
      call configure_transport_failure(0, 0)
      call configure_solver_failure(SOLVER_ERROR_FOR_REDUCTION)
      call reset_projection_spy()
      call reset_reorder_spy()

      call run_solve(F, F2, Ft, Ft2, status)

      call assert_true('failed Allreduce preserves an earlier solver status', &
                       status == SOLVER_ERROR_FOR_REDUCTION)
      call assert_true('failed Allreduce aborts once with the MPI status', &
                       abort_calls == 1 .and. &
                       abort_code == ALLREDUCE_AFTER_SOLVER_ERROR .and. &
                       abort_comm == MPI_COMM_WORLD)
      call assert_true('failed Allreduce skips Bcast, scatter, and reorder', &
                       allreduce_calls == 3 .and. bcast_calls == 0 .and. &
                       gather_call_count == 1 .and. compute_calls == 1 .and. &
                       solve_calls == 1 .and. scatter_call_count == 0 .and. &
                       reorder_calls == 0)
      call assert_true('failed Allreduce does not publish a partial solve', &
                       F(1, 1) == 41.d0 .and. F2(1, 1) == OUTPUT_SENTINEL)
      call assert_true('post-solver Allreduce uses the expected contract', &
                       collective_contract_ok)
   end subroutine test_allreduce_failure_preserves_solver_error


   subroutine test_bcast_failure_preserves_solver_error()
      real(kind=8), allocatable :: F(:, :), F2(:, :), Ft(:, :), Ft2(:, :)
      integer(kind=4) :: status

      call prepare_case(F, F2)
      call configure_collective_failures(0, 0, BCAST_ERROR)
      call configure_transport_failure(0, 0)
      call configure_solver_failure(SOLVER_ERROR_FOR_BCAST)
      call reset_projection_spy()
      call reset_reorder_spy()

      call run_solve(F, F2, Ft, Ft2, status)

      call assert_true('Bcast failure preserves the earlier solver status', &
                       status == SOLVER_ERROR_FOR_BCAST)
      call assert_true('Bcast failure aborts once with the MPI status', &
                       abort_calls == 1 .and. abort_code == BCAST_ERROR .and. &
                       abort_comm == MPI_COMM_WORLD)
      call assert_true('Bcast failure stops before scatter and reorder', &
                       allreduce_calls == 3 .and. bcast_calls == 1 .and. &
                       bcast_root == 0 .and. gather_call_count == 1 .and. &
                       compute_calls == 1 .and. solve_calls == 1 .and. &
                       scatter_call_count == 0 .and. reorder_calls == 0)
      call assert_true('Bcast failure does not publish a partial solve', &
                       F(1, 1) == 41.d0 .and. F2(1, 1) == OUTPUT_SENTINEL)
      call assert_true('Bcast failure uses the expected collective contract', &
                       collective_contract_ok)
   end subroutine test_bcast_failure_preserves_solver_error


   subroutine run_solve(F, F2, Ft, Ft2, status)
      real(kind=8), allocatable, intent(inout) :: F(:, :), F2(:, :)
      real(kind=8), allocatable, intent(inout) :: Ft(:, :), Ft2(:, :)
      integer(kind=4), intent(out) :: status
      integer(kind=4), parameter :: direction(3) = (/0, 0, 0/)
      real(kind=8), parameter :: mixA(4) = (/1.d0, 0.d0, 0.d0, 0.d0/)
      real(kind=8), parameter :: mixB(4) = 0.d0
      real(kind=8), parameter :: mixBT(4) = 0.d0

      call solve_problem(ads_test, ads_trial, 1, 2, 3, mixA, mixB, mixBT, &
                         direction, .false., F, F2, Ft, Ft2, status)
   end subroutine run_solve


   subroutine prepare_case(F, F2)
      real(kind=8), allocatable, intent(out) :: F(:, :), F2(:, :)

      allocate(F(1, 1), F2(1, 1))
      F = 41.d0
      F2 = OUTPUT_SENTINEL
   end subroutine prepare_case


   subroutine prepare_space(ads)
      type(ADS_Setup), intent(out) :: ads

      ads%n = 0
      ads%p = 0
      ads%s = 1
      ads%nelem = 1
      ads%nrcpp = 1
      ads%ibeg = 1
      ads%iend = 1
      allocate(ads%Ux(0:1), ads%Uy(0:1), ads%Uz(0:1))
      ads%Ux = (/0.d0, 1.d0/)
      ads%Uy = ads%Ux
      ads%Uz = ads%Ux
      allocate(ads%dimensionsX(1), ads%dimensionsY(1), ads%dimensionsZ(1))
      allocate(ads%shiftsX(1), ads%shiftsY(1), ads%shiftsZ(1))
      ads%dimensionsX = 1
      ads%dimensionsY = 1
      ads%dimensionsZ = 1
      ads%shiftsX = 0
      ads%shiftsY = 0
      ads%shiftsZ = 0
   end subroutine prepare_space


   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (condition) then
         write (*, '(A)') 'PASS '//trim(label)
      else
         failures = failures + 1
         write (*, '(A)') 'FAIL '//trim(label)
      end if
   end subroutine assert_true

end program synchronization_error_probe
