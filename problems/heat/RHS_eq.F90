module RHS_eq

implicit none

contains

! -------------------------------------------------------------------
! Right-hand side of the equation.
!
! Input:
! ------
! ads             - ADS setup structure
! X_              - quadrature points
! k_              - indexes for quadrature points
! e_              - indexes for elements
! a_              - indexes of basis functions
! ads_data        - data structures for ADS
! J               - jacobian
! W               - weight for quadratures
!
! Output:
! -------
! ret             - value of RHS function at given point
!
! -------------------------------------------------------------------

subroutine ComputePointForRHS( &
ads, &
X, &
k, &
e, &
a, &
ads_data, J, W, ret)
use Setup, ONLY: ADS_Setup,ADS_compute_data
use input_data
use projection_engine, ONLY: global2local
implicit none
type (ADS_setup), intent(in) :: ads
real   (kind=8), intent(in), dimension(3)  :: X
integer(kind=4), intent(in), dimension(3)  :: k
integer(kind=4), intent(in), dimension(3)  :: e
integer(kind=4), intent(in), dimension(3)  :: a
type (ADS_compute_data), intent(in) :: ads_data
real   (kind=8), intent(in)  :: J,W
real (kind = 8), intent(out) :: ret
real   (kind=8) :: fval,kqval
real   (kind=8) :: dvx,dvy,dvz,rhs,v
integer(kind = 4) :: rx, ry, rz, ix, iy, iz, sx, sy, sz
integer(kind = 4) :: bx, by, bz
real   (kind=8), dimension(3)  :: du
real (kind = 8) :: dux, duy, duz
integer(kind = 4) :: mx, my, mz, ngx, ngy, ngz
integer(kind = 4) :: ind, ind1, ind23, indx, indy, indz
integer(kind = 4) :: indbx, indby, indbz
real (kind = 8) :: Uval, ucoeff

Uval = 0
dux = 0
duy = 0
duz = 0
do bx = 0, ads % p(1)
   do by = 0, ads % p(2)
      do bz = 0, ads % p(3)
         ind = (ads % Ox(e(1)) + bx) + (ads % Oy(e(2)) + by)*(ads % n(1) + 1) + (ads % Oz(e(3)) + bz)* &
         (ads % n(2) + 1)*(ads % n(1) + 1)
         call global2local(ind, ads % n, indbx, indby, indbz)

         rx = 2
         ry = 2
         rz = 2
         if (indbx < ads % ibeg(1) - 1) rx = 1
         if (indbx > ads % iend(1) - 1) rx = 3
         if (indby < ads % ibeg(2) - 1) ry = 1
         if (indby > ads % iend(2) - 1) ry = 3
         if (indbz < ads % ibeg(3) - 1) rz = 1
         if (indbz > ads % iend(3) - 1) rz = 3

         ix = indbx - ads % ibegsx(rx) + 1
         iy = indby - ads % ibegsy(ry) + 1
         iz = indbz - ads % ibegsz(rz) + 1
         sx = ads % iendsx(rx) - ads % ibegsx(rx) + 1
         sy = ads % iendsy(ry) - ads % ibegsy(ry) + 1
         sz = ads % iendsz(rz) - ads % ibegsz(rz) + 1
         ind = ix + sx * (iy + sy * iz)

#ifdef IDEBUG
         if (ind < 0 .or. ind > ads % nrcpp(3) * ads % nrcpp(1) * ads % nrcpp(2) - 1) then
            write(ERROR_UNIT, *) PRINTRANK, 'Oh crap', ix, iy, iz
            write(ERROR_UNIT, *) PRINTRANK, 'r', rx, ry, rz
            write(ERROR_UNIT, *) PRINTRANK, 'x', ads % ibeg(1), ads % iend(1)
            write(ERROR_UNIT, *) PRINTRANK, 'y', ads % ibeg(2), ads % iend(2)
            write(ERROR_UNIT, *) PRINTRANK, 'z', ads % ibeg(3), ads % iend(3)
            write(ERROR_UNIT, *) PRINTRANK, 'idxs', indx, indy, indz
            write(ERROR_UNIT, *) PRINTRANK, 'sizes=', sx, sy, sz
            write(ERROR_UNIT, *) PRINTRANK, 'begsx=', ads % ibegsx
            write(ERROR_UNIT, *) PRINTRANK, 'endsx=', ads % iendsx
            write(ERROR_UNIT, *) PRINTRANK, 'begsy=', ads % ibegsy
            write(ERROR_UNIT, *) PRINTRANK, 'endsy=', ads % iendsy
            write(ERROR_UNIT, *) PRINTRANK, 'begsz=', ads % ibegsz
            write(ERROR_UNIT, *) PRINTRANK, 'endsz=', ads % iendsz
         endif
#endif
                                 
         Ucoeff = ads_data % R(ind + 1, rx, ry, rz)
         v = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
         dvx = ads % NNx(1, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
         dvy = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(1, by, k(2), e(2)) * ads % NNz(0, bz, k(3), e(3))
         dvz = ads % NNx(0, bx, k(1), e(1)) * ads % NNy(0, by, k(2), e(2)) * ads % NNz(1, bz, k(3), e(3))

         Uval = Uval + Ucoeff * v
         dux = dux + Ucoeff * dvx
         duy = duy + Ucoeff * dvy
         duz = duz + Ucoeff * dvz
      enddo
   enddo
enddo

du = (/ dux, duy, duz /)

v   = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3))
dvx = ads % NNx(1,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvy = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(1,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvz = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(1,a(3),k(3),e(3)) 

kqval = 1.d0
fval = 0.d0
!--- Real
if (t > 0.0) then
  rhs = Dt * ( - kqval  * (du(1)*dvx + du(2)*dvy + du(3)*dvz) + v * fval)
  ret = J*W*(v * Uval + rhs)

  l2norm = l2norm + J*W*v*Uval*Uval
else
  fval = initial_state(X(1),X(2),X(3))
  ret= J*W*v*fval
  l2norm = l2norm + J*W*v*fval*fval
endif

end subroutine


end module