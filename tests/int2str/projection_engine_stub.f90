! Fail fast if the unrelated legacy Norm_L2 implementation is ever entered.
module projection_engine
   implicit none

contains

   subroutine global2local(ind, n, x, y, z)
      use ISO_FORTRAN_ENV, only: ERROR_UNIT
      integer(kind=4), intent(in) :: ind
      integer(kind=4), intent(in), dimension(3) :: n
      integer(kind=4), intent(out) :: x, y, z

      x = 0
      y = 0
      z = 0
      write (ERROR_UNIT, *) 'Unexpected call to legacy global2local:', ind, n
      stop 99
   end subroutine global2local

end module projection_engine
