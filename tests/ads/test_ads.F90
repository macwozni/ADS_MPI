program test_ads
   use Setup, only: ADS_Setup, ADS_compute_data
   use ADSS, only: Sub_Step, Step, MultiStep
   use workflow_test_support, only: reset_workflow_spy, fail_solve, &
      solve_call_count, form_un_call_count, form_un_substeps, &
      form_rhs_call_count, normalize_call_count, distribute_call_count, forcing
   implicit none

   integer(kind=4) :: checks, failures

   checks = 0
   failures = 0

   call test_substep_stops_on_first_solve_error()
   call test_substep_stops_on_second_solve_error()
   call test_substep_stops_on_third_solve_error()
   call test_substep_success_path()
   call test_step_success_publishes_numeric_result()
   call test_step_preserves_error_and_cleans_transients()
   call test_multistep_preserves_error_and_cleans_transients()
   call test_multistep_stops_in_third_substep()
   call test_multistep_recovers_after_error()

   if (failures == 0) then
      write (*, '(A,I0,A)') 'OK (', checks, ' ADS error-propagation checks)'
   else
      write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                 ' ADS error-propagation checks)'
      stop 1
   end if

contains

   subroutine test_substep_stops_on_first_solve_error()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_directional_buffers(data, ads%s)
      call reset_workflow_spy()
      call fail_solve(1, -101)
      call fail_solve(2, -999)

      call invoke_substep(ads, data, status)

      call assert_true('Sub_Step returns the first directional-solve error', &
                       status == -101)
      call assert_true('Sub_Step stops after a failing first solve', &
                       solve_call_count == 1)
      call assert_true('Sub_Step skips normalization and distribution after first error', &
                       normalize_call_count == 0 .and. distribute_call_count == 0)
      call assert_true('Sub_Step forms its RHS exactly once before the error', &
                       form_rhs_call_count == 1)
   end subroutine test_substep_stops_on_first_solve_error


   subroutine test_substep_stops_on_second_solve_error()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_directional_buffers(data, ads%s)
      call reset_workflow_spy()
      call fail_solve(2, -202)
      call fail_solve(3, -999)

      call invoke_substep(ads, data, status)

      call assert_true('Sub_Step preserves an error from the second solve', &
                       status == -202)
      call assert_true('Sub_Step does not start the third solve after an error', &
                       solve_call_count == 2)
      call assert_true('Sub_Step does not publish a partially solved RHS', &
                       normalize_call_count == 0 .and. distribute_call_count == 0)
   end subroutine test_substep_stops_on_second_solve_error


   subroutine test_substep_stops_on_third_solve_error()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_directional_buffers(data, ads%s)
      call reset_workflow_spy()
      call fail_solve(3, -203)

      call invoke_substep(ads, data, status)

      call assert_true('Sub_Step preserves an error from the third solve', &
                       status == -203 .and. solve_call_count == 3)
      call assert_true('Sub_Step does not publish a failed third solve', &
                       normalize_call_count == 0 .and. distribute_call_count == 0)
   end subroutine test_substep_stops_on_third_solve_error


   subroutine test_substep_success_path()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_directional_buffers(data, ads%s)
      call reset_workflow_spy()

      call invoke_substep(ads, data, status)

      call assert_true('Sub_Step still completes all three solves on success', &
                       status == 0 .and. solve_call_count == 3)
      call assert_true('Sub_Step publishes only a fully solved result', &
                       normalize_call_count == 1 .and. distribute_call_count == 1)
   end subroutine test_substep_success_path


   subroutine test_step_success_publishes_numeric_result()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_step_state(data)
      call reset_workflow_spy()

      call Step(6, forcing, ads, data, status)

      call assert_true('Step success runs all directional solves', &
                       status == 0 .and. solve_call_count == 3)
      call assert_true('Step publishes the completely solved numeric buffer', &
                       allocated(data%FF) .and. all(data%FF == 3.0d0))
      call assert_true('Step releases all temporary solve buffers on success', &
                       .not. allocated(data%F) .and. .not. allocated(data%F2) .and. &
                       .not. allocated(data%F3))
      call assert_true('Step normalizes and distributes exactly once', &
                       normalize_call_count == 1 .and. distribute_call_count == 1)
   end subroutine test_step_success_publishes_numeric_result


   subroutine test_step_preserves_error_and_cleans_transients()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_step_state(data)
      allocate(data%Ft(1, 1), data%Ft2(1, 1), data%Ft3(1, 1))
      allocate(data%FF(1, 1), data%FFt(1, 1))
      data%Ft = 7.0d0
      data%Ft2 = 8.0d0
      data%Ft3 = 9.0d0
      data%FF = 111.0d0
      data%FFt = 222.0d0

      call reset_workflow_spy()
      call fail_solve(1, -301)
      call fail_solve(2, -999)
      call Step(7, forcing, ads, data, status)

      call assert_true('Step returns the first solve error unchanged', &
                       status == -301 .and. solve_call_count == 1)
      call assert_true('Step cleans every transient directional buffer on error', &
                       transient_buffers_are_clean(data))
      call assert_true('Step keeps the last complete solution on error', &
                       allocated(data%FF) .and. allocated(data%FFt) .and. &
                       all(data%FF == 111.0d0) .and. all(data%FFt == 222.0d0))
      call assert_true('Step does not normalize or distribute a partial result', &
                       normalize_call_count == 0 .and. distribute_call_count == 0)
      call assert_true('Step reconstructs once but starts no later solve', &
                       form_un_call_count == 1 .and. form_un_substeps(1) == 1)
   end subroutine test_step_preserves_error_and_cleans_transients


   subroutine test_multistep_preserves_error_and_cleans_transients()
      type(ADS_Setup) :: ads_test, ads_trial
      type(ADS_compute_data) :: data
      real(kind=8) :: mix(4, 3), alpha_step(7, 3)
      integer(kind=4) :: rhs_state(6, 3), status

      call initialize_setup_stub(ads_test)
      call initialize_setup_stub(ads_trial)
      mix = 1.0d0
      alpha_step = 2.0d0
      rhs_state = 3

      call reset_workflow_spy()
      call fail_solve(5, -405)
      call fail_solve(6, -999)
      call MultiStep(11, mix, forcing, ads_test, ads_trial, data, 2, &
                     alpha_step, status, rhs_du_state=rhs_state)

      call assert_true('MultiStep returns an error from a later substep unchanged', &
                       status == -405)
      call assert_true('MultiStep stops before a solve following the first error', &
                       solve_call_count == 5)
      call assert_true('MultiStep cleans all current-substep transient buffers', &
                       transient_buffers_are_clean(data))
      call assert_true('MultiStep retains the last fully completed intermediate state', &
                       allocated(data%FF) .and. allocated(data%FFt) .and. &
                       all(data%FF == 3.0d0) .and. all(data%FFt == 3.0d0))
      call assert_true('MultiStep publishes only its fully completed first substep', &
                       normalize_call_count == 1 .and. distribute_call_count == 1)
      call assert_true('MultiStep stops reconstruction at the failed substep', &
                       form_un_call_count == 2 .and. &
                       all(form_un_substeps(1:2) == (/1, 2/)))
      call assert_true('MultiStep clears the transient RHS derivative selector', &
                       all(data%rhs_du_state == 0))
   end subroutine test_multistep_preserves_error_and_cleans_transients


   subroutine test_multistep_stops_in_third_substep()
      type(ADS_Setup) :: ads_test, ads_trial
      type(ADS_compute_data) :: data
      real(kind=8) :: mix(4, 3), alpha_step(7, 3)
      integer(kind=4) :: status

      call initialize_setup_stub(ads_test)
      call initialize_setup_stub(ads_trial)
      mix = 1.0d0
      alpha_step = 1.0d0
      call reset_workflow_spy()
      call fail_solve(9, -409)

      call MultiStep(13, mix, forcing, ads_test, ads_trial, data, 3, &
                     alpha_step, status)

      call assert_true('MultiStep preserves a third-solve error in its third substep', &
                       status == -409 .and. solve_call_count == 9)
      call assert_true('MultiStep publishes only two complete substeps before that error', &
                       normalize_call_count == 2 .and. distribute_call_count == 2)
      call assert_true('MultiStep retains the second complete intermediate state', &
                       allocated(data%FF) .and. allocated(data%FFt) .and. &
                       all(data%FF == 6.0d0) .and. all(data%FFt == 6.0d0))
      call assert_true('MultiStep cleans failed third-substep temporaries', &
                       transient_buffers_are_clean(data))
   end subroutine test_multistep_stops_in_third_substep


   subroutine test_multistep_recovers_after_error()
      type(ADS_Setup) :: ads_test, ads_trial
      type(ADS_compute_data) :: data
      real(kind=8) :: mix(4, 3), alpha_step(7, 3)
      integer(kind=4) :: status

      call initialize_setup_stub(ads_test)
      call initialize_setup_stub(ads_trial)
      mix = 1.0d0
      alpha_step = 1.0d0
      call reset_workflow_spy()

      call MultiStep(12, mix, forcing, ads_test, ads_trial, data, 3, &
                     alpha_step, status)

      call assert_true('MultiStep remains reusable after an earlier injected error', &
                       status == 0 .and. solve_call_count == 9)
      call assert_true('Successful MultiStep completes and publishes all substeps', &
                       normalize_call_count == 3 .and. distribute_call_count == 3)
      call assert_true('Successful MultiStep leaves only persistent result buffers', &
                       transient_buffers_are_clean(data) .and. &
                       allocated(data%FF) .and. allocated(data%FFt))
      call assert_true('Successful MultiStep publishes the ninth solve result', &
                       all(data%FF == 9.0d0) .and. all(data%FFt == 9.0d0))
   end subroutine test_multistep_recovers_after_error


   subroutine invoke_substep(ads, data, status)
      type(ADS_Setup), intent(in) :: ads
      type(ADS_compute_data), intent(inout) :: data
      integer(kind=4), intent(out) :: status
      real(kind=8) :: mix(4), alpha_step(7, 3)
      integer(kind=4) :: direction(3), abc(3, 3)

      mix = (/1.0d0, 2.0d0, 3.0d0, 4.0d0/)
      alpha_step = 1.0d0
      direction = 0
      abc(:, 1) = (/1, 2, 3/)
      abc(:, 2) = (/2, 1, 3/)
      abc(:, 3) = (/3, 1, 2/)
      call Sub_Step(ads, ads, 1, mix, direction, 1, abc, 1, &
                    alpha_step, forcing, data, status)
   end subroutine invoke_substep


   subroutine initialize_setup_stub(ads)
      type(ADS_Setup), intent(out) :: ads

      ads%s = (/2, 2, 2/)
   end subroutine initialize_setup_stub


   subroutine allocate_directional_buffers(data, extents)
      type(ADS_compute_data), intent(inout) :: data
      integer(kind=4), intent(in) :: extents(3)

      allocate(data%F(extents(1), extents(2)*extents(3)))
      allocate(data%F2(extents(2), extents(3)*extents(1)))
      allocate(data%F3(extents(3), extents(1)*extents(2)))
      allocate(data%Ft(extents(1), extents(2)*extents(3)))
      allocate(data%Ft2(extents(2), extents(3)*extents(1)))
      allocate(data%Ft3(extents(3), extents(1)*extents(2)))
   end subroutine allocate_directional_buffers


   subroutine allocate_step_state(data)
      type(ADS_compute_data), intent(inout) :: data

      allocate(data%Un13(1, 1, 1, 1, 1, 1))
      allocate(data%Un23(1, 1, 1, 1, 1, 1))
      data%Un13 = 13.0d0
      data%Un23 = 23.0d0
   end subroutine allocate_step_state


   logical function transient_buffers_are_clean(data) result(clean)
      type(ADS_compute_data), intent(in) :: data

      clean = .not. allocated(data%F) .and. .not. allocated(data%F2) .and. &
              .not. allocated(data%F3) .and. .not. allocated(data%Ft) .and. &
              .not. allocated(data%Ft2) .and. .not. allocated(data%Ft3)
   end function transient_buffers_are_clean


   subroutine assert_true(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition

      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         write (*, '(A)') 'FAIL: '//label
      end if
   end subroutine assert_true

end program test_ads
