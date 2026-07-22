module rhs_eq_probe_callbacks
   implicit none

contains

   function zero_forcing(un, derivative, point) result(value)
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: derivative(3), point(3)
      real(kind=8) :: value

      value = 0.d0*un + 0.d0*sum(derivative) + 0.d0*sum(point)
   end function zero_forcing

end module rhs_eq_probe_callbacks


program rhs_eq_probe
   use ISO_FORTRAN_ENV, only: ERROR_UNIT
   use Setup, only: ADS_Setup, ADS_compute_data
   use RHS_eq, only: ComputePointForRHS
   use rhs_eq_probe_callbacks, only: zero_forcing
   implicit none

   character(len=32) :: mode
   type(ADS_Setup) :: ads
   type(ADS_compute_data) :: ads_data
   integer(kind=4) :: a(3), direction(3), e(3), k(3), substep
   real(kind=8) :: alpha_step(7, 3), du(3), ret, X(3)

   if (command_argument_count() /= 1) stop 90
   call get_command_argument(1, mode)

   a = 0
   e = 1
   k = 1
   direction = 0
   du = (/ 1.d0, 2.d0, 3.d0 /)
   X = (/ 4.d0, 5.d0, 6.d0 /)
   alpha_step = 0.d0
   ads%tau = 1.d0
   allocate(ads%NNx(0:1, 0:0, 1:1, 1:1))
   allocate(ads%NNy(0:1, 0:0, 1:1, 1:1))
   allocate(ads%NNz(0:1, 0:0, 1:1, 1:1))
   ads%NNx = 1.d0
   ads%NNy = 1.d0
   ads%NNz = 1.d0
   ads_data%state_mine = 1
   ads_data%rhs_du_state = 0

   select case (trim(mode))
   case ("substep-zero")
      substep = 0
   case ("substep-four")
      substep = 4
   case ("state-negative")
      substep = 1
      ads_data%rhs_du_state(1, 1) = -1
   case ("state-four")
      substep = 1
      ads_data%rhs_du_state(6, 1) = 4
   case default
      write(ERROR_UNIT, '(A)') "unknown probe mode"
      stop 90
   end select

   call ComputePointForRHS(ads, X, k, e, a, du, 0, 7.d0, 11.d0, 13.d0, &
                           ads_data, 1.d0, 1.d0, direction, substep, alpha_step, &
                           zero_forcing, ret)

   write(ERROR_UNIT, '(A)') "validation unexpectedly succeeded"
   stop 98

end program rhs_eq_probe
