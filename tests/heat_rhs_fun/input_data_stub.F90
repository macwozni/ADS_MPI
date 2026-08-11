module input_data
   implicit none

   real(kind=8) :: t = 0.0d0

contains

   function initial_state(x, y, z) result(value)
      real(kind=8), intent(in) :: x, y, z
      real(kind=8) :: value

      value = 1.0d0 + 2.0d0*x - 3.0d0*y + 4.0d0*z
   end function initial_state

end module input_data
