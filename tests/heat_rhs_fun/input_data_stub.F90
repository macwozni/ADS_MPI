module input_data
   implicit none

   real(kind=8) :: t = 0.0d0

contains

   function initial_state(x, y, z) result(value)
      real(kind=8), intent(in) :: x, y, z
      real(kind=8) :: dx, dy, dz, r2, value

      ! Profile used by iga-ads examples/heat/heat_3d.hpp.  Keeping the
      ! compatibility datum in this stub gives the callback test an exact,
      ! self-contained numerical oracle.
      dx = x - 0.5d0
      dy = y - 0.5d0
      dz = z - 0.5d0
      r2 = min(8.0d0*(dx*dx + dy*dy + dz*dz), 1.0d0)
      value = (r2 - 1.0d0)**2*(r2 + 1.0d0)**2
   end function initial_state

end module input_data
