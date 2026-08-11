program transport_error_probe
   use Setup, only: ADS_Setup
   use parallelism, only: InitializeParallelism, Cleanup_Parallelism
   use communicators, only: CreateCommunicators, Cleanup_Communicators
   use ads_directional_solve, only: solve_problem
   use directional_test_support, only: configure_directional_spies, &
      compute_calls, solve_calls, compute_contract_ok, solve_contract_ok
   use my_mpi, only: configure_transport_failure, transport_operation_count, &
      gather_call_count, scatter_call_count
   implicit none

   real(kind=8), parameter :: SENTINEL = -987654321.25d0
   type(ADS_Setup) :: ads_test, ads_trial
   integer(kind=4) :: checks, failures, ierr

   checks = 0
   failures = 0
   call InitializeParallelism(1, 1, 1, ierr)
   call assert_true('single-rank parallelism initializes', ierr == 0)
   call CreateCommunicators(ierr)
   call assert_true('single-rank communicators initialize', ierr == 0)

   call prepare_space(ads_trial, (/1, 1, 1/), 1000)
   call prepare_space(ads_test, (/2, 1, 1/), 2000)

   call exercise_transport_failure('trial Gather', 1, -611, 1, 0, 0, 0)
   call exercise_transport_failure('iGRM Gather', 2, -622, 2, 0, 0, 0)
   call exercise_transport_failure('iGRM Scatter', 3, -633, 2, 1, 1, 1)
   call exercise_transport_failure('trial Scatter', 4, -644, 2, 2, 1, 1)

   call release_space(ads_test)
   call release_space(ads_trial)
   call Cleanup_Communicators(ierr)
   call assert_true('single-rank communicators clean up', ierr == 0)
   call Cleanup_Parallelism(ierr)
   call assert_true('single-rank parallelism cleans up', ierr == 0)

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' directional transport-error checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' directional transport-error checks)'
      stop 1
   end if

contains

   subroutine exercise_transport_failure(label, operation, injected_status, &
                                         expected_gathers, expected_scatters, &
                                         expected_computes, expected_solves)
      character(len=*), intent(in) :: label
      integer(kind=4), intent(in) :: operation, injected_status
      integer(kind=4), intent(in) :: expected_gathers, expected_scatters
      integer(kind=4), intent(in) :: expected_computes, expected_solves
      integer(kind=4), parameter :: direction(3) = (/1, 0, 0/)
      real(kind=8), parameter :: mixA(4) = (/1.d0, 2.d0, 3.d0, 4.d0/)
      real(kind=8), parameter :: mixB(4) = (/5.d0, 6.d0, 7.d0, 8.d0/)
      real(kind=8), parameter :: mixBT(4) = (/9.d0, 10.d0, 11.d0, 12.d0/)
      real(kind=8), allocatable :: F(:, :), F2(:, :), Ft(:, :), Ft2(:, :)
      real(kind=8), allocatable :: F_before(:, :), Ft_before(:, :)
      real(kind=8), allocatable :: expected_packed(:, :)
      integer(kind=4) :: status
      logical :: contract_ok

      call allocate_buffers(F, F2, Ft, Ft2)
      allocate(F_before(size(F, 1), size(F, 2)))
      allocate(Ft_before(size(Ft, 1), size(Ft, 2)))
      F_before = F
      Ft_before = Ft
      allocate(expected_packed(size(F, 1) + size(Ft, 1), size(F, 2)))
      expected_packed(1:size(F, 1), :) = F
      expected_packed(size(F, 1) + 1:size(expected_packed, 1), :) = Ft

      call configure_directional_spies(ads_test, ads_trial, 1, mixA, mixB, &
                                       mixBT, .false., expected_packed)
      call configure_transport_failure(operation, injected_status)

      call solve_problem(ads_test, ads_trial, 1, 2, 3, mixA, mixB, mixBT, &
                         direction, .true., F, F2, Ft, Ft2, status)

      call assert_true(trim(label)//' preserves the first transport error', &
                       status == injected_status)
      call assert_true(trim(label)//' stops at the exact transport stage', &
                       transport_operation_count == operation .and. &
                       gather_call_count == expected_gathers .and. &
                       scatter_call_count == expected_scatters)
      call assert_true(trim(label)//' skips all later compute/solve stages', &
                       compute_calls == expected_computes .and. &
                       solve_calls == expected_solves)
      contract_ok = expected_computes == 0 .or. compute_contract_ok
      contract_ok = contract_ok .and. (expected_solves == 0 .or. solve_contract_ok)
      call assert_true(trim(label)//' keeps successful earlier-stage contracts exact', &
                       contract_ok)
      call assert_true(trim(label)//' does not publish or reorder a partial result', &
                       all(F == F_before) .and. all(Ft == Ft_before) .and. &
                       all(F2 == SENTINEL) .and. all(Ft2 == SENTINEL))
   end subroutine exercise_transport_failure


   subroutine allocate_buffers(F, F2, Ft, Ft2)
      real(kind=8), allocatable, intent(out) :: F(:, :), F2(:, :)
      real(kind=8), allocatable, intent(out) :: Ft(:, :), Ft2(:, :)
      integer(kind=4) :: column, row

      allocate(F(ads_trial%s(1), ads_trial%s(2)*ads_trial%s(3)))
      allocate(F2(ads_trial%s(2), ads_trial%s(1)*ads_trial%s(3)))
      allocate(Ft(ads_test%s(1), ads_trial%s(2)*ads_trial%s(3)))
      allocate(Ft2(ads_trial%s(2), ads_test%s(1)*ads_trial%s(3)))

      do column = 1, size(F, 2)
         do row = 1, size(F, 1)
            F(row, column) = real(100*column + row, kind=8)
         end do
      end do
      do column = 1, size(Ft, 2)
         do row = 1, size(Ft, 1)
            Ft(row, column) = real(1000 + 100*column + row, kind=8)
         end do
      end do
      F2 = SENTINEL
      Ft2 = SENTINEL
   end subroutine allocate_buffers


   subroutine prepare_space(ads, n, marker)
      type(ADS_Setup), intent(out) :: ads
      integer(kind=4), intent(in) :: n(3), marker
      integer(kind=4) :: i

      ads%n = n
      ads%p = 0
      ads%nelem = n + 1
      ads%m = ads%n + ads%p + 1
      ads%nrcpp = n + 1
      ads%ibeg = 1
      ads%iend = n + 1
      ads%s = n + 1

      allocate(ads%dimensionsX(1), ads%dimensionsY(1), ads%dimensionsZ(1))
      allocate(ads%shiftsX(1), ads%shiftsY(1), ads%shiftsZ(1))
      ads%dimensionsX = ads%s(1)*ads%s(2)*ads%s(3)
      ads%dimensionsY = ads%dimensionsX
      ads%dimensionsZ = ads%dimensionsX
      ads%shiftsX = 0
      ads%shiftsY = 0
      ads%shiftsZ = 0

      allocate(ads%Ux(0:ads%n(1) + ads%p(1) + 1))
      allocate(ads%Uy(0:ads%n(2) + ads%p(2) + 1))
      allocate(ads%Uz(0:ads%n(3) + ads%p(3) + 1))
      do i = 0, ubound(ads%Ux, 1)
         ads%Ux(i) = real(marker + 100 + i, kind=8)
      end do
      do i = 0, ubound(ads%Uy, 1)
         ads%Uy(i) = real(marker + 200 + i, kind=8)
      end do
      do i = 0, ubound(ads%Uz, 1)
         ads%Uz(i) = real(marker + 300 + i, kind=8)
      end do
   end subroutine prepare_space


   subroutine release_space(ads)
      type(ADS_Setup), intent(inout) :: ads

      if (allocated(ads%Ux)) deallocate(ads%Ux)
      if (allocated(ads%Uy)) deallocate(ads%Uy)
      if (allocated(ads%Uz)) deallocate(ads%Uz)
      if (allocated(ads%dimensionsX)) deallocate(ads%dimensionsX)
      if (allocated(ads%dimensionsY)) deallocate(ads%dimensionsY)
      if (allocated(ads%dimensionsZ)) deallocate(ads%dimensionsZ)
      if (allocated(ads%shiftsX)) deallocate(ads%shiftsX)
      if (allocated(ads%shiftsY)) deallocate(ads%shiftsY)
      if (allocated(ads%shiftsZ)) deallocate(ads%shiftsZ)
   end subroutine release_space


   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         write (*, '(A,A)') 'FAIL: ', trim(label)
      end if
   end subroutine assert_true

end program transport_error_probe
