module basis
   implicit none

contains

   function EvalSpline(d, &
                       Ux, px, nx, nelemx, &
                       Uy, py, ny, nelemy, &
                       Uz, pz, nz, nelemz, &
                       coeffs, x, y, z) result(value)
      integer(kind=4), intent(in) :: d
      integer(kind=4), intent(in) :: px, nx, nelemx
      integer(kind=4), intent(in) :: py, ny, nelemy
      integer(kind=4), intent(in) :: pz, nz, nelemz
      real(kind=8), intent(in) :: Ux(0:nx + px + 1)
      real(kind=8), intent(in) :: Uy(0:ny + py + 1)
      real(kind=8), intent(in) :: Uz(0:nz + pz + 1)
      real(kind=8), intent(in) :: coeffs(0:nx, 0:ny, 0:nz)
      real(kind=8), intent(in) :: x, y, z
      real(kind=8) :: value

      value = coeffs(0, 0, 0) + 100.d0*x + 10.d0*y + z
      value = value + 0.d0*real(d + px + nx + nelemx + py + ny + nelemy + &
                                pz + nz + nelemz, kind=8)
      value = value + 0.d0*(Ux(0) + Uy(0) + Uz(0))
   end function EvalSpline

end module basis
