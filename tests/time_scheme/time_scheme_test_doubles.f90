module projection_engine
   implicit none

   integer(kind=4) :: validation_call_count = 0
   integer(kind=4) :: recorded_test_p(3) = 0
   integer(kind=4) :: recorded_test_n(3) = 0
   integer(kind=4) :: recorded_test_nelem(3) = 0
   integer(kind=4) :: recorded_test_size(3) = 0
   integer(kind=4) :: recorded_trial_p(3) = 0
   integer(kind=4) :: recorded_trial_n(3) = 0
   integer(kind=4) :: recorded_trial_nelem(3) = 0
   integer(kind=4) :: recorded_trial_size(3) = 0
   real(kind=8) :: recorded_test_sum(3) = 0.0d0
   real(kind=8) :: recorded_trial_sum(3) = 0.0d0

contains

   subroutine reset_validation_spy()
      validation_call_count = 0
      recorded_test_p = 0
      recorded_test_n = 0
      recorded_test_nelem = 0
      recorded_test_size = 0
      recorded_trial_p = 0
      recorded_trial_n = 0
      recorded_trial_nelem = 0
      recorded_trial_size = 0
      recorded_test_sum = 0.0d0
      recorded_trial_sum = 0.0d0
   end subroutine reset_validation_spy


   subroutine ValidateIGRMMesh(U_test, p_test, n_test, nelem_test, &
                               U_trial, p_trial, n_trial, nelem_trial)
      integer(kind=4), intent(in) :: p_test, n_test, nelem_test
      integer(kind=4), intent(in) :: p_trial, n_trial, nelem_trial
      real(kind=8), intent(in) :: U_test(0:n_test + p_test + 1)
      real(kind=8), intent(in) :: U_trial(0:n_trial + p_trial + 1)
      integer(kind=4) :: call_index

      validation_call_count = validation_call_count + 1
      call_index = validation_call_count
      if (call_index > size(recorded_test_p)) return

      recorded_test_p(call_index) = p_test
      recorded_test_n(call_index) = n_test
      recorded_test_nelem(call_index) = nelem_test
      recorded_test_size(call_index) = size(U_test)
      recorded_trial_p(call_index) = p_trial
      recorded_trial_n(call_index) = n_trial
      recorded_trial_nelem(call_index) = nelem_trial
      recorded_trial_size(call_index) = size(U_trial)
      recorded_test_sum(call_index) = sum(U_test)
      recorded_trial_sum(call_index) = sum(U_trial)
   end subroutine ValidateIGRMMesh

end module projection_engine


module ADSS
   use Setup, only: ADS_Setup, ADS_compute_data
   use Interfaces, only: forcing_fun, rhs_point_fun
   implicit none

   integer(kind=4), parameter :: STEP_STATUS = 731
   integer(kind=4), parameter :: MULTISTEP_STATUS = 947

   integer(kind=4) :: step_call_count = 0
   integer(kind=4) :: multistep_call_count = 0
   integer(kind=4) :: recorded_iter = 0
   integer(kind=4) :: recorded_history_index = 0
   integer(kind=4) :: recorded_test_n(3) = 0
   integer(kind=4) :: recorded_trial_n(3) = 0
   logical :: recorded_rhs_point_present = .false.
   logical :: recorded_lhs_mix_present = .false.
   logical :: recorded_rhs_du_state_present = .false.
   real(kind=8) :: recorded_forcing_value = 0.0d0
   real(kind=8) :: recorded_rhs_point_value = 0.0d0
   real(kind=8) :: recorded_ads_time = 0.0d0
   real(kind=8) :: recorded_mix(4, 3) = 0.0d0
   real(kind=8) :: recorded_alpha_step(7, 3) = 0.0d0
   real(kind=8) :: recorded_lhs_mix(4, 3, 3) = 0.0d0
   integer(kind=4) :: recorded_rhs_du_state(6, 3) = 0

contains

   subroutine reset_adss_spy()
      step_call_count = 0
      multistep_call_count = 0
      recorded_iter = 0
      recorded_history_index = 0
      recorded_test_n = 0
      recorded_trial_n = 0
      recorded_rhs_point_present = .false.
      recorded_lhs_mix_present = .false.
      recorded_rhs_du_state_present = .false.
      recorded_forcing_value = 0.0d0
      recorded_rhs_point_value = 0.0d0
      recorded_ads_time = 0.0d0
      recorded_mix = 0.0d0
      recorded_alpha_step = 0.0d0
      recorded_lhs_mix = 0.0d0
      recorded_rhs_du_state = 0
   end subroutine reset_adss_spy


   subroutine Step(iter, RHS_fun, ads, ads_data, mierr, RHS_point)
      integer(kind=4), intent(in) :: iter
      procedure(forcing_fun) :: RHS_fun
      type(ADS_Setup), intent(in) :: ads
      type(ADS_compute_data), intent(inout) :: ads_data
      integer(kind=4), intent(out) :: mierr
      procedure(rhs_point_fun), optional :: RHS_point
      real(kind=8) :: du(3), X(3)
      real(kind=8) :: alpha_step(7, 3)
      integer(kind=4) :: indices(3)

      step_call_count = step_call_count + 1
      recorded_iter = iter
      recorded_test_n = ads%n
      recorded_trial_n = ads%n
      recorded_rhs_point_present = present(RHS_point)
      recorded_ads_time = ads_data%t

      du = (/ 1.0d0, 2.0d0, 3.0d0 /)
      X = (/ 4.0d0, 5.0d0, 6.0d0 /)
      recorded_forcing_value = RHS_fun(7.0d0, du, X)
      if (present(RHS_point)) then
         indices = 0
         alpha_step = 0.0d0
         call RHS_point(ads, X, indices, indices, indices, du, 0, 8.0d0, 0.0d0, 0.0d0, &
                        ads_data, 1.0d0, 1.0d0, indices, 1, alpha_step, RHS_fun, &
                        recorded_rhs_point_value)
      end if
      mierr = STEP_STATUS
   end subroutine Step


   subroutine MultiStep(iter, mix, RHS_fun, ads_test, ads_trial, ads_data, n, alpha_step, mierr, &
                        RHS_point, lhs_mix, rhs_du_state)
      integer(kind=4), intent(in) :: iter, n
      real(kind=8), intent(in) :: mix(4, 3)
      procedure(forcing_fun) :: RHS_fun
      type(ADS_Setup), intent(in) :: ads_test, ads_trial
      type(ADS_compute_data), intent(inout) :: ads_data
      real(kind=8), intent(in) :: alpha_step(7, 3)
      integer(kind=4), intent(out) :: mierr
      procedure(rhs_point_fun), optional :: RHS_point
      real(kind=8), intent(in), optional :: lhs_mix(4, 3, 3)
      integer(kind=4), intent(in), optional :: rhs_du_state(6, 3)
      real(kind=8) :: du(3), X(3)
      integer(kind=4) :: indices(3)

      multistep_call_count = multistep_call_count + 1
      recorded_iter = iter
      recorded_history_index = n
      recorded_test_n = ads_test%n
      recorded_trial_n = ads_trial%n
      recorded_rhs_point_present = present(RHS_point)
      recorded_lhs_mix_present = present(lhs_mix)
      recorded_rhs_du_state_present = present(rhs_du_state)
      recorded_ads_time = ads_data%t
      recorded_mix = mix
      recorded_alpha_step = alpha_step
      if (present(lhs_mix)) recorded_lhs_mix = lhs_mix
      if (present(rhs_du_state)) recorded_rhs_du_state = rhs_du_state

      du = (/ 1.0d0, 2.0d0, 3.0d0 /)
      X = (/ 4.0d0, 5.0d0, 6.0d0 /)
      recorded_forcing_value = RHS_fun(7.0d0, du, X)
      if (present(RHS_point)) then
         indices = 0
         call RHS_point(ads_test, X, indices, indices, indices, du, n, 8.0d0, 0.0d0, 0.0d0, &
                        ads_data, 1.0d0, 1.0d0, indices, 1, alpha_step, RHS_fun, &
                        recorded_rhs_point_value)
      end if
      mierr = MULTISTEP_STATUS
   end subroutine MultiStep

end module ADSS
