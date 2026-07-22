program test_ads_directional_solve
   use mpi
   use Setup, only: ADS_Setup
   use parallelism, only: MYRANK, MYRANKX, MYRANKY, MYRANKZ, &
      NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints, FillDimVector, &
      InitializeParallelism, Cleanup_Parallelism
   use communicators, only: CreateCommunicators, Cleanup_Communicators
   use ads_directional_solve, only: solve_problem
   use directional_test_support, only: configure_directional_spies, &
      compute_calls, solve_calls, compute_contract_ok, solve_contract_ok
   implicit none

   integer(kind=4), parameter :: TRIAL_TAG = 1
   integer(kind=4), parameter :: TEST_TAG = 2
   real(kind=8), parameter :: SENTINEL = -987654321.25d0

   type(ADS_Setup) :: ads_test, ads_trial
   integer(kind=4) :: process_grid(3), direction(3)
   integer(kind=4) :: checks, failures, ierr, axis, enrichment_axis
   character(len=96) :: label

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), &
                              process_grid(3), ierr)
   call CreateCommunicators(ierr)

   checks = 0
   failures = 0
   call prepare_space(ads_trial, (/4, 5, 6/), (/1, 2, 3/), TRIAL_TAG)
   call prepare_space(ads_test, (/4, 5, 6/), (/2, 4, 5/), TEST_TAG)

   direction = 0
   do axis = 1, 3
      write (label, '(A,I0)') 'ADS axis ', axis
      call exercise_direction(trim(label), axis, direction, .false., &
                              failures)
   end do

   do enrichment_axis = 1, 3
      direction = 0
      direction(enrichment_axis) = 1
      do axis = 1, 3
         if (axis == enrichment_axis) then
            write (label, '(A,I0)') 'iGRM active axis ', axis
         else
            write (label, '(A,I0,A,I0)') 'iGRM enrichment ', &
               enrichment_axis, ' transverse to solve axis ', axis
         end if
         call exercise_direction(trim(label), axis, direction, .true., &
                                 failures)
      end do
   end do

   if (MYRANK == 0) then
      if (failures == 0) then
         write (*, '(A,I0,A)') 'OK (', checks, &
                               ' ADS directional solve checks)'
      else
         write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                    ' ADS directional solve checks)'
      end if
   end if

   call release_space(ads_test)
   call release_space(ads_trial)
   call Cleanup_Communicators(ierr)
   call Cleanup_Parallelism(ierr)

   if (failures /= 0) stop 1

contains

   subroutine exercise_direction(case_label, a, selected_direction, &
                                 igrm, failure_count)
      character(len=*), intent(in) :: case_label
      integer(kind=4), intent(in) :: a, selected_direction(3)
      logical, intent(in) :: igrm
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: b, c, current_order(3), next_order(3)
      integer(kind=4) :: mixed_beg(3), mixed_extents(3)
      integer(kind=4) :: expected_local_calls, solve_ierr
      real(kind=8) :: mixA(4), mixB(4), mixBT(4)
      real(kind=8), allocatable :: F(:, :), F2(:, :), Ft(:, :), Ft2(:, :)
      real(kind=8), allocatable :: F_before(:, :), Ft_before(:, :)
      real(kind=8), allocatable :: expected_F2(:, :), expected_Ft2(:, :)
      real(kind=8), allocatable :: expected_packed(:, :)
      logical :: equ, local_ok

      call axis_orders(a, b, c, current_order, next_order)
      call build_encoded_buffer(ads_trial%ibeg, ads_trial%s, current_order, &
                                TRIAL_TAG, F)
      allocate(F_before(size(F, 1), size(F, 2)))
      F_before = F
      call build_encoded_buffer(ads_trial%ibeg, ads_trial%s, next_order, &
                                TRIAL_TAG, expected_F2)
      expected_F2 = 2.d0*expected_F2
      allocate(F2(size(expected_F2, 1), size(expected_F2, 2)))
      F2 = SENTINEL

      mixed_beg = merge(ads_test%ibeg, ads_trial%ibeg, &
                        selected_direction == 1)
      mixed_extents = merge(ads_test%s, ads_trial%s, &
                            selected_direction == 1)
      if (igrm) then
         call build_encoded_buffer(mixed_beg, mixed_extents, current_order, &
                                   TEST_TAG, Ft)
         allocate(Ft_before(size(Ft, 1), size(Ft, 2)))
         Ft_before = Ft
         call build_encoded_buffer(mixed_beg, mixed_extents, next_order, &
                                   TEST_TAG, expected_Ft2)
         expected_Ft2 = 2.d0*expected_Ft2
         allocate(Ft2(size(expected_Ft2, 1), size(expected_Ft2, 2)))
         Ft2 = SENTINEL
      end if

      call build_expected_packed(a, b, c, selected_direction, igrm, &
                                 expected_packed)
      equ = selected_direction(a) /= 1
      call case_mixes(a, selected_direction, mixA, mixB, mixBT)
      call configure_directional_spies(ads_test, ads_trial, a, mixA, mixB, &
                                       mixBT, equ, expected_packed)

      solve_ierr = -1
      call solve_problem(ads_test, ads_trial, a, b, c, mixA, mixB, mixBT, &
                         selected_direction, igrm, F, F2, Ft, Ft2, &
                         solve_ierr)

      local_ok = solve_ierr == 0
      local_ok = local_ok .and. exact_array(F, F_before)
      local_ok = local_ok .and. exact_array(F2, expected_F2)
      if (igrm) then
         local_ok = local_ok .and. exact_array(Ft, Ft_before)
         local_ok = local_ok .and. exact_array(Ft2, expected_Ft2)
      else
         local_ok = local_ok .and. .not. allocated(Ft)
         local_ok = local_ok .and. .not. allocated(Ft2)
      end if
      call assert_global(case_label//' returns exact reordered buffers', &
                         local_ok, failure_count)

      local_ok = compute_contract_ok .and. solve_contract_ok
      call assert_global(case_label//' packs and forwards the exact contract', &
                         local_ok, failure_count)

      expected_local_calls = 0
      if (axis_rank(a) == 0) expected_local_calls = 1
      local_ok = compute_calls == expected_local_calls .and. &
                 solve_calls == expected_local_calls
      call assert_global(case_label//' solves only on the processor face', &
                         local_ok, failure_count)
   end subroutine exercise_direction


   subroutine build_expected_packed(a, b, c, selected_direction, igrm, &
                                    packed)
      integer(kind=4), intent(in) :: a, b, c, selected_direction(3)
      logical, intent(in) :: igrm
      real(kind=8), allocatable, intent(out) :: packed(:, :)
      integer(kind=4) :: mixed_beg(3), mixed_extents(3), mixed_rows
      integer(kind=4) :: order(3)
      real(kind=8), allocatable :: trial_block(:, :), test_block(:, :)

      order = (/a, b, c/)
      call build_gathered_block(ads_trial%n(a) + 1, ads_trial%ibeg, &
                                ads_trial%s, order, TRIAL_TAG, trial_block)

      if (.not. igrm) then
         allocate(packed(size(trial_block, 1), size(trial_block, 2)))
         packed = trial_block
         return
      end if

      mixed_beg = merge(ads_test%ibeg, ads_trial%ibeg, &
                        selected_direction == 1)
      mixed_extents = merge(ads_test%s, ads_trial%s, &
                            selected_direction == 1)
      mixed_rows = (1 - selected_direction(a))*(ads_trial%n(a) + 1) + &
                   selected_direction(a)*(ads_test%n(a) + 1)
      call build_gathered_block(mixed_rows, mixed_beg, mixed_extents, &
                                order, TEST_TAG, test_block)

      if (selected_direction(a) == 1) then
         allocate(packed(size(trial_block, 1) + size(test_block, 1), &
                         size(trial_block, 2)))
         packed(1:size(trial_block, 1), :) = trial_block
         packed(size(trial_block, 1) + 1:size(packed, 1), :) = test_block
      else
         allocate(packed(size(trial_block, 1), &
                         size(trial_block, 2) + size(test_block, 2)))
         packed(:, 1:size(trial_block, 2)) = trial_block
         packed(:, size(trial_block, 2) + 1:size(packed, 2)) = test_block
      end if
   end subroutine build_expected_packed


   subroutine build_gathered_block(global_rows, begins, extents, order, &
                                   tag, values)
      integer(kind=4), intent(in) :: global_rows, begins(3), extents(3)
      integer(kind=4), intent(in) :: order(3), tag
      real(kind=8), allocatable, intent(out) :: values(:, :)
      integer(kind=4) :: coords(3), row, j, k, column

      allocate(values(global_rows, extents(order(2))*extents(order(3))))
      do k = 1, extents(order(3))
         coords(order(3)) = begins(order(3)) + k - 1
         do j = 1, extents(order(2))
            coords(order(2)) = begins(order(2)) + j - 1
            column = j + (k - 1)*extents(order(2))
            do row = 1, global_rows
               coords(order(1)) = row
               values(row, column) = encoded_value(tag, coords)
            end do
         end do
      end do
   end subroutine build_gathered_block


   subroutine build_encoded_buffer(begins, extents, order, tag, values)
      integer(kind=4), intent(in) :: begins(3), extents(3), order(3), tag
      real(kind=8), allocatable, intent(out) :: values(:, :)
      integer(kind=4) :: coords(3), i, j, k, column

      allocate(values(extents(order(1)), &
                      extents(order(2))*extents(order(3))))
      do k = 1, extents(order(3))
         coords(order(3)) = begins(order(3)) + k - 1
         do j = 1, extents(order(2))
            coords(order(2)) = begins(order(2)) + j - 1
            column = j + (k - 1)*extents(order(2))
            do i = 1, extents(order(1))
               coords(order(1)) = begins(order(1)) + i - 1
               values(i, column) = encoded_value(tag, coords)
            end do
         end do
      end do
   end subroutine build_encoded_buffer


   real(kind=8) function encoded_value(tag, coords) result(value)
      integer(kind=4), intent(in) :: tag, coords(3)

      value = real(tag*1000000 + coords(1)*10000 + coords(2)*100 + &
                   coords(3), kind=8)
   end function encoded_value


   subroutine prepare_space(ads, element_count, degree, tag)
      type(ADS_Setup), intent(out) :: ads
      integer(kind=4), intent(in) :: element_count(3), degree(3), tag
      integer(kind=4) :: ranks(3), process_counts(3), axis, i

      ads%nelem = element_count
      ads%p = degree
      ads%n = element_count + degree - 1
      ads%m = ads%n + ads%p + 1
      ranks = (/MYRANKX, MYRANKY, MYRANKZ/)
      process_counts = (/NRPROCX, NRPROCY, NRPROCZ/)

      do axis = 1, 3
         call ComputeEndpoints(ranks(axis), process_counts(axis), &
            ads%n(axis), ads%p(axis), ads%nrcpp(axis), ads%ibeg(axis), &
            ads%iend(axis), ads%mine(axis), ads%maxe(axis))
      end do
      ads%s = ads%iend - ads%ibeg + 1
      ads%lnelem = ads%maxe - ads%mine + 1

      call FillDimVector(ads%dimensionsX, ads%shiftsX, ads%nrcpp(1), &
                         ads%s(2)*ads%s(3), ads%n(1), NRPROCX)
      call FillDimVector(ads%dimensionsY, ads%shiftsY, ads%nrcpp(2), &
                         ads%s(1)*ads%s(3), ads%n(2), NRPROCY)
      call FillDimVector(ads%dimensionsZ, ads%shiftsZ, ads%nrcpp(3), &
                         ads%s(1)*ads%s(2), ads%n(3), NRPROCZ)

      allocate(ads%Ux(0:ads%n(1) + ads%p(1) + 1))
      allocate(ads%Uy(0:ads%n(2) + ads%p(2) + 1))
      allocate(ads%Uz(0:ads%n(3) + ads%p(3) + 1))
      do i = 0, ubound(ads%Ux, 1)
         ads%Ux(i) = real(tag*1000 + 100 + i, kind=8)
      end do
      do i = 0, ubound(ads%Uy, 1)
         ads%Uy(i) = real(tag*1000 + 200 + i, kind=8)
      end do
      do i = 0, ubound(ads%Uz, 1)
         ads%Uz(i) = real(tag*1000 + 300 + i, kind=8)
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


   subroutine axis_orders(axis, b, c, current_order, next_order)
      integer(kind=4), intent(in) :: axis
      integer(kind=4), intent(out) :: b, c, current_order(3), next_order(3)

      select case (axis)
      case (1)
         b = 2
         c = 3
         current_order = (/1, 2, 3/)
         next_order = (/2, 1, 3/)
      case (2)
         b = 1
         c = 3
         current_order = (/2, 1, 3/)
         next_order = (/3, 1, 2/)
      case (3)
         b = 1
         c = 2
         current_order = (/3, 1, 2/)
         next_order = (/1, 2, 3/)
      case default
         stop 72
      end select
   end subroutine axis_orders


   subroutine case_mixes(axis, selected_direction, mixA, mixB, mixBT)
      integer(kind=4), intent(in) :: axis, selected_direction(3)
      real(kind=8), intent(out) :: mixA(4), mixB(4), mixBT(4)
      real(kind=8) :: marker

      marker = real(axis + 10*dot_product(selected_direction, (/1, 2, 3/)), &
                    kind=8)
      mixA = marker + (/0.125d0, 0.25d0, 0.5d0, 0.75d0/)
      mixB = marker + (/10.125d0, 10.25d0, 10.5d0, 10.75d0/)
      mixBT = marker + (/20.125d0, 20.25d0, 20.5d0, 20.75d0/)
   end subroutine case_mixes


   integer(kind=4) function axis_rank(axis) result(rank)
      integer(kind=4), intent(in) :: axis

      select case (axis)
      case (1)
         rank = MYRANKX
      case (2)
         rank = MYRANKY
      case (3)
         rank = MYRANKZ
      case default
         stop 73
      end select
   end function axis_rank


   logical function exact_array(actual, expected) result(matches)
      real(kind=8), intent(in) :: actual(:, :), expected(:, :)

      matches = all(shape(actual) == shape(expected))
      if (.not. matches) return
      matches = all(actual == expected)
   end function exact_array


   subroutine assert_global(assertion_label, local_ok, failure_count)
      character(len=*), intent(in) :: assertion_label
      logical, intent(in) :: local_ok
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: local_flag, global_flag, reduce_ierr

      checks = checks + 1
      local_flag = 0
      if (local_ok) local_flag = 1
      call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, &
                         MPI_COMM_WORLD, reduce_ierr)
      if (global_flag == 0) then
         failure_count = failure_count + 1
         if (MYRANK == 0) write (*, '(A,A)') 'FAIL: ', trim(assertion_label)
      end if
   end subroutine assert_global


   subroutine read_process_grid(grid)
      integer(kind=4), intent(out) :: grid(3)
      character(len=64) :: argument
      integer(kind=4) :: i, parse_status

      if (command_argument_count() /= 3) then
         write (*, '(A)') 'usage: test_ads_directional_solve procx procy procz'
         stop 2
      end if
      do i = 1, 3
         grid(i) = 0
         call get_command_argument(i, argument)
         read (argument, *, iostat=parse_status) grid(i)
         if (parse_status /= 0) then
            write (*, '(A,I0)') 'invalid process-grid argument ', i
            stop 2
         end if
         if (grid(i) <= 0) then
            write (*, '(A,I0)') 'invalid process-grid argument ', i
            stop 2
         end if
      end do
   end subroutine read_process_grid

end program test_ads_directional_solve
