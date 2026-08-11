module math
   implicit none

contains

   function lerp(t, x, y) result(value)
      real(kind=8), intent(in) :: t, x, y
      real(kind=8) :: value

      value = (1.d0 - t)*x + t*y
   end function lerp

end module math
