module Interfaces

   interface
      function forcing_fun (un, X) result (ret)
         implicit none
         real (kind = 8), intent(in), dimension(:) :: un
         real (kind = 8), intent(in), dimension(3) :: X
         real (kind = 8) :: ret
      end function

   end interface

end module Interfaces