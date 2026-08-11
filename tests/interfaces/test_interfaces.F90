program test_interfaces
   use Setup, only: ADS_Setup, ADS_compute_data
   use Interfaces, only: forcing_fun, rhs_point_fun
   implicit none

   procedure(forcing_fun), pointer :: forcing_callback
   procedure(rhs_point_fun), pointer :: rhs_callback
   type(ADS_Setup) :: ads
   type(ADS_compute_data) :: data
   real(kind=8) :: X(3), du(3), alpha_step(7, 3), result, expected
   integer(kind=4) :: k(3), e(3), a(3), direction(3)

   forcing_callback => sample_forcing
   rhs_callback => sample_rhs_point

   ads%tau = 0.5d0
   data%t = 0.25d0
   X = (/1.d0, 2.d0, 3.d0/)
   du = (/4.d0, 5.d0, 6.d0/)
   k = (/1, 2, 3/)
   e = (/4, 5, 6/)
   a = (/0, 1, 2/)
   direction = (/1, 0, 0/)
   alpha_step = 0.d0
   alpha_step(1, 1) = 8.d0

   call rhs_callback(ads, X, k, e, a, du, 7, 1.5d0, 2.5d0, 3.5d0, &
                     data, 2.d0, 0.25d0, direction, 2, alpha_step, &
                     forcing_callback, result)

   expected = 2.d0*0.25d0*sample_forcing(1.5d0, du, X) + &
              real(sum(k + e + a + direction) + 7 + 2, kind=8) + &
              2.5d0 + 3.5d0 + 8.d0 + ads%tau + data%t

   if (abs(result - expected) > 1.d-12) then
      write (*, '(A,2ES24.15)') 'FAIL: callback contract mismatch: ', result, expected
      stop 1
   end if

   write (*, '(A)') 'OK (forcing_fun and rhs_point_fun callback contracts)'

contains

   function sample_forcing(un, gradient, point) result(value)
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: gradient(3), point(3)
      real(kind=8) :: value

      value = un + sum(gradient) + sum(point)
   end function sample_forcing

   subroutine sample_rhs_point(ads_arg, point, quadrature, element, basis_index, gradient, &
                               history_index, un11, un13, un23, data_arg, J, W, &
                               direction_arg, substep, alpha, forcing, value)
      use Setup, only: ADS_Setup, ADS_compute_data
      use Interfaces, only: forcing_fun
      type(ADS_Setup), intent(in) :: ads_arg
      real(kind=8), intent(in) :: point(3), gradient(3)
      integer(kind=4), intent(in) :: quadrature(3), element(3), basis_index(3)
      integer(kind=4), intent(in) :: history_index, direction_arg(3), substep
      real(kind=8), intent(in) :: un11, un13, un23, J, W, alpha(7, 3)
      type(ADS_compute_data), intent(in) :: data_arg
      procedure(forcing_fun) :: forcing
      real(kind=8), intent(out) :: value

      value = J*W*forcing(un11, gradient, point) + &
              real(sum(quadrature + element + basis_index + direction_arg) + &
                   history_index + substep, kind=8) + &
              un13 + un23 + alpha(1, 1) + ads_arg%tau + data_arg%t
   end subroutine sample_rhs_point

end program test_interfaces
