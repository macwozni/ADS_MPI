module math
   implicit none

   integer(kind=4) :: bump3d_calls = 0
   real(kind=8) :: recorded_r = 0.d0
   real(kind=8) :: recorded_Rr = 0.d0
   real(kind=8) :: recorded_x = 0.d0
   real(kind=8) :: recorded_y = 0.d0
   real(kind=8) :: recorded_z = 0.d0

contains

   subroutine reset_bump3d_stub
      bump3d_calls = 0
      recorded_r = 0.d0
      recorded_Rr = 0.d0
      recorded_x = 0.d0
      recorded_y = 0.d0
      recorded_z = 0.d0
   end subroutine reset_bump3d_stub


   function bump3d(r, Rr, x, y, z) result(val)
      real(kind=8), intent(in) :: r, Rr, x, y, z
      real(kind=8) :: val

      bump3d_calls = bump3d_calls + 1
      recorded_r = r
      recorded_Rr = Rr
      recorded_x = x
      recorded_y = y
      recorded_z = z
      val = 7.25d0
   end function bump3d

end module math
