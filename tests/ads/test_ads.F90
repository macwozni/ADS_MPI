program test_ads
   use Setup, only: ADS_Setup, ADS_compute_data
   use ADSS, only: Sub_Step, Step, MultiStep
   use workflow_test_support, only: reset_workflow_spy, fail_solve, &
      solve_call_count, form_un_call_count, form_un_substeps, &
      form_rhs_call_count, normalize_call_count, distribute_call_count, forcing, &
      custom_rhs_point, recorded_solve_axes, recorded_solve_directions, &
      recorded_solve_igrm, recorded_solve_mixA, recorded_solve_mixB, &
      recorded_solve_mixBT, recorded_rhs_point_present, recorded_rhs_point_value
   implicit none

   integer(kind=4) :: checks, failures

   checks = 0
   failures = 0

   call test_substep_stops_on_first_solve_error()
   call test_substep_stops_on_second_solve_error()
   call test_substep_stops_on_third_solve_error()
   call test_substep_success_path()
   call test_substep_uses_mass_mix_in_second_solve()
   call test_substep_uses_mass_mix_in_third_solve()
   call test_step_success_publishes_numeric_result()
   call test_step_forwards_rhs_point()
   call test_step_preserves_error_and_cleans_transients()
   call test_multistep_stops_in_first_substep()
   call test_multistep_preserves_error_and_cleans_transients()
   call test_multistep_stops_in_third_substep()
   call test_multistep_recovers_after_error()
   call test_multistep_forwards_rhs_point_and_lhs_mix()

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
      call assert_true('Sub_Step forwards the exact legacy mix to every solve', &
                       three_solve_contract_matches((/1.0d0, 2.0d0, 3.0d0, 4.0d0/)))
   end subroutine test_substep_success_path


   subroutine test_substep_uses_mass_mix_in_second_solve()
      integer(kind=4), parameter :: abc(3, 3) = reshape((/ &
         2, 1, 3,  1, 2, 3,  3, 1, 2 /), (/3, 3/))

      call verify_substep_mass_mix_position(abc, 2, &
         'Sub_Step selects the mass operator in its second solve')
   end subroutine test_substep_uses_mass_mix_in_second_solve


   subroutine test_substep_uses_mass_mix_in_third_solve()
      integer(kind=4), parameter :: abc(3, 3) = reshape((/ &
         2, 1, 3,  3, 1, 2,  1, 2, 3 /), (/3, 3/))

      call verify_substep_mass_mix_position(abc, 3, &
         'Sub_Step selects the mass operator in its third solve')
   end subroutine test_substep_uses_mass_mix_in_third_solve


   subroutine verify_substep_mass_mix_position(abc, mass_position, label)
      integer(kind=4), intent(in) :: abc(3, 3), mass_position
      character(len=*), intent(in) :: label
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      real(kind=8), parameter :: mass_mix(4) = (/1.0d0, 0.0d0, 0.0d0, 0.0d0/)
      real(kind=8) :: mix(4), alpha_step(7, 3), lhs_mix(4, 3)
      integer(kind=4) :: call_index, direction(3), status
      logical :: contract_matches

      call initialize_setup_stub(ads)
      call allocate_directional_buffers(data, ads%s)
      mix = -999.0d0
      alpha_step = 1.0d0
      direction = (/1, 0, 0/)
      lhs_mix(:, 1) = (/11.0d0, 12.0d0, 13.0d0, 14.0d0/)
      lhs_mix(:, 2) = (/21.0d0, 22.0d0, 23.0d0, 24.0d0/)
      lhs_mix(:, 3) = (/31.0d0, 32.0d0, 33.0d0, 34.0d0/)
      call reset_workflow_spy()

      call Sub_Step(ads, ads, 1, mix, direction, 1, abc, 1, &
                    alpha_step, forcing, data, status, lhs_mix=lhs_mix)

      contract_matches = status == 0 .and. solve_call_count == 3
      do call_index = 1, 3
         contract_matches = contract_matches .and. &
            all(recorded_solve_axes(:, call_index) == abc(:, call_index))
         contract_matches = contract_matches .and. &
            all(recorded_solve_mixB(:, call_index) == &
                lhs_mix(:, abc(1, call_index)))
         contract_matches = contract_matches .and. &
            all(recorded_solve_mixBT(:, call_index) == &
                lhs_mix(:, abc(1, call_index)))
         if (call_index == mass_position) then
            contract_matches = contract_matches .and. &
               all(recorded_solve_mixA(:, call_index) == mass_mix)
         else
            contract_matches = contract_matches .and. &
               all(recorded_solve_mixA(:, call_index) == &
                   lhs_mix(:, abc(1, call_index)))
         end if
      end do

      call assert_true(label, contract_matches)
   end subroutine verify_substep_mass_mix_position


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
      call assert_true('Step leaves the optional point callback absent', &
                       .not. any(recorded_rhs_point_present))
   end subroutine test_step_success_publishes_numeric_result


   subroutine test_step_forwards_rhs_point()
      type(ADS_Setup) :: ads
      type(ADS_compute_data) :: data
      integer(kind=4) :: status

      call initialize_setup_stub(ads)
      call allocate_step_state(data)
      call reset_workflow_spy()

      call Step(8, forcing, ads, data, status, custom_rhs_point)

      call assert_true('Step forwards an optional point callback to RHS assembly', &
                       status == 0 .and. form_rhs_call_count == 1 .and. &
                       recorded_rhs_point_present(1))
      call assert_true('Step forwards the exact point callback implementation', &
                       recorded_rhs_point_value(1) == 1001.0d0)
   end subroutine test_step_forwards_rhs_point


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


   subroutine test_multistep_stops_in_first_substep()
      type(ADS_Setup) :: ads_test, ads_trial
      type(ADS_compute_data) :: data
      real(kind=8) :: mix(4, 3), alpha_step(7, 3)
      integer(kind=4) :: status

      call initialize_setup_stub(ads_test)
      call initialize_setup_stub(ads_trial)
      mix = 1.0d0
      alpha_step = 1.0d0
      call reset_workflow_spy()
      call fail_solve(1, -401)

      call MultiStep(9, mix, forcing, ads_test, ads_trial, data, 1, &
                     alpha_step, status)

      call assert_true('MultiStep preserves an error from its first substep', &
                       status == -401 .and. solve_call_count == 1)
      call assert_true('MultiStep publishes nothing after a first-substep error', &
                       normalize_call_count == 0 .and. distribute_call_count == 0 .and. &
                       .not. allocated(data%FF) .and. .not. allocated(data%FFt))
      call assert_true('MultiStep cleans first-substep transient buffers', &
                       transient_buffers_are_clean(data))
      call assert_true('MultiStep stops reconstruction after its first substep', &
                       form_un_call_count == 1 .and. form_un_substeps(1) == 1)
      call assert_true('MultiStep resets RHS derivative state after an early error', &
                       all(data%rhs_du_state == 0))
   end subroutine test_multistep_stops_in_first_substep


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


   subroutine test_multistep_forwards_rhs_point_and_lhs_mix()
      type(ADS_Setup) :: ads_test, ads_trial
      type(ADS_compute_data) :: data
      real(kind=8) :: mix(4, 3), alpha_step(7, 3), lhs_mix(4, 3, 3)
      integer(kind=4) :: axis, coefficient, status, substep

      call initialize_setup_stub(ads_test)
      call initialize_setup_stub(ads_trial)
      mix = -999.0d0
      alpha_step = 0.0d0
      do substep = 1, 3
         do axis = 1, 3
            do coefficient = 1, 4
               lhs_mix(coefficient, axis, substep) = &
                  real(100*substep + 10*axis + coefficient, kind=8)
            end do
         end do
      end do
      call reset_workflow_spy()

      call MultiStep(15, mix, forcing, ads_test, ads_trial, data, 4, &
                     alpha_step, status, RHS_point=custom_rhs_point, &
                     lhs_mix=lhs_mix)

      call assert_true('MultiStep forwards RHS_point through all three assemblies', &
                       status == 0 .and. form_rhs_call_count == 3 .and. &
                       all(recorded_rhs_point_present))
      call assert_true('MultiStep forwards the exact RHS_point callback', &
                       all(recorded_rhs_point_value == (/1001.0d0, 1002.0d0, 1003.0d0/)))
      call assert_true('MultiStep forwards exact directional LHS coefficients', &
                       multistep_lhs_contract_matches(lhs_mix))
   end subroutine test_multistep_forwards_rhs_point_and_lhs_mix


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


   logical function three_solve_contract_matches(expected_mix) result(matches)
      real(kind=8), intent(in) :: expected_mix(4)
      integer(kind=4), parameter :: expected_axes(3, 3) = reshape((/ &
         1, 2, 3, &
         2, 1, 3, &
         3, 1, 2 /), (/3, 3/))
      integer(kind=4) :: call_index

      matches = solve_call_count == 3
      do call_index = 1, 3
         matches = matches .and. all(recorded_solve_axes(:, call_index) == &
                                     expected_axes(:, call_index))
         matches = matches .and. all(recorded_solve_directions(:, call_index) == 0)
         matches = matches .and. .not. recorded_solve_igrm(call_index)
         matches = matches .and. all(recorded_solve_mixA(:, call_index) == expected_mix)
         matches = matches .and. all(recorded_solve_mixB(:, call_index) == expected_mix)
         matches = matches .and. all(recorded_solve_mixBT(:, call_index) == expected_mix)
      end do
   end function three_solve_contract_matches


   logical function multistep_lhs_contract_matches(lhs_mix) result(matches)
      real(kind=8), intent(in) :: lhs_mix(4, 3, 3)
      integer(kind=4), parameter :: expected_axes(3, 9) = reshape((/ &
         1, 2, 3,  2, 1, 3,  3, 1, 2, &
         2, 1, 3,  3, 1, 2,  1, 2, 3, &
         3, 1, 2,  1, 2, 3,  2, 1, 3 /), (/3, 9/))
      real(kind=8), parameter :: mass_mix(4) = (/1.0d0, 0.0d0, 0.0d0, 0.0d0/)
      integer(kind=4) :: axis, call_index, substep

      matches = solve_call_count == 9
      call_index = 0
      do substep = 1, 3
         do axis = 1, 3
            call_index = call_index + 1
            matches = matches .and. all(recorded_solve_axes(:, call_index) == &
                                        expected_axes(:, call_index))
            matches = matches .and. all(recorded_solve_directions(:, call_index) == &
                                        merge(1, 0, (/1, 2, 3/) == substep))
            matches = matches .and. recorded_solve_igrm(call_index)
            matches = matches .and. all(recorded_solve_mixB(:, call_index) == &
                                        lhs_mix(:, expected_axes(1, call_index), substep))
            matches = matches .and. all(recorded_solve_mixBT(:, call_index) == &
                                        lhs_mix(:, expected_axes(1, call_index), substep))
            if (expected_axes(1, call_index) == substep) then
               matches = matches .and. all(recorded_solve_mixA(:, call_index) == mass_mix)
            else
               matches = matches .and. all(recorded_solve_mixA(:, call_index) == &
                                           lhs_mix(:, expected_axes(1, call_index), substep))
            end if
         end do
      end do
   end function multistep_lhs_contract_matches


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
