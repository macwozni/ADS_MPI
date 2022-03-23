module Interfaces

   interface
      function forcing_fun(un, du, X) result(ret)
         implicit none
         real(kind=8), intent(in) :: un
         real(kind=8), intent(in), dimension(3) :: du
         real(kind=8), intent(in), dimension(3) :: X
         real(kind=8) :: ret
      end function

   end interface

end module Interfaces
