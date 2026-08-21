module input_data
   implicit none

   integer(kind=4) :: body_force_calls = 0
   real(kind=8) :: recorded_X(3) = 0.d0

contains

   function body_force(X) result(value)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value(3)

      body_force_calls = body_force_calls + 1
      recorded_X = X
      value(1) = 17.d0 + X(1) + 2.d0*X(2) + 3.d0*X(3)
      value(2) = 31.d0 - 2.d0*X(1) + X(2) - X(3)
      value(3) = -5.d0 + 4.d0*X(1) - 3.d0*X(2) + 2.d0*X(3)
   end function body_force

end module input_data
