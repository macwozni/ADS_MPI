module input_data
   implicit none

   integer(kind=4) :: exact_forcing_calls = 0
   real(kind=8) :: recorded_X(3) = 0.d0

contains

   function exact_forcing(X) result(value)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      exact_forcing_calls = exact_forcing_calls + 1
      recorded_X = X
      value = 17.d0 + X(1) + 2.d0*X(2) + 3.d0*X(3)
   end function exact_forcing

end module input_data
