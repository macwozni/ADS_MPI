! Minimal implementation of the index conversion used by legacy Norm_L2.
module projection_engine
   implicit none

   integer(kind=4) :: global2local_call_count = 0

contains

   subroutine reset_global2local_spy()
      global2local_call_count = 0
   end subroutine reset_global2local_spy


   subroutine global2local(ind, n, x, y, z)
      integer(kind=4), intent(in) :: ind
      integer(kind=4), intent(in), dimension(3) :: n
      integer(kind=4), intent(out) :: x, y, z

      integer(kind=4) :: remainder

      global2local_call_count = global2local_call_count + 1
      z = ind/((n(1) + 1)*(n(2) + 1))
      remainder = ind - z*(n(1) + 1)*(n(2) + 1)
      y = remainder/(n(1) + 1)
      x = remainder - y*(n(1) + 1)
   end subroutine global2local

end module projection_engine
