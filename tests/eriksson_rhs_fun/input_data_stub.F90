module input_data
   implicit none

   real(kind=8) :: t = 0.d0
   integer(kind=4) :: initial_state_calls = 0
   real(kind=8) :: recorded_coordinates(3) = 0.d0

contains

   function initial_state(x, y, z) result(value)
      real(kind=8), intent(in) :: x, y, z
      real(kind=8) :: value

      initial_state_calls = initial_state_calls + 1
      recorded_coordinates = (/ x, y, z /)
      value = 1.d0 + 2.d0*x - 3.d0*y + 4.d0*z
   end function initial_state

end module input_data
