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
! du_             - value of derivative from previous time step
! Uval            - previous solution coefficient at given point
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
du, &
Uval, &
ads_data, J, W, l2norm, ret)
use Setup, ONLY: ADS_Setup,ADS_compute_data
use input_data
use projection_engine, ONLY: global2local
implicit none
type (ADS_setup), intent(in) :: ads
real   (kind=8), intent(in), dimension(3)  :: X
integer(kind=4), intent(in), dimension(3)  :: k
integer(kind=4), intent(in), dimension(3)  :: e
integer(kind=4), intent(in), dimension(3)  :: a
real   (kind=8), intent(in), dimension(3)  :: du
real (kind = 8), intent(in) :: Uval
type (ADS_compute_data), intent(in) :: ads_data
real   (kind=8), intent(in)  :: J,W
real (kind = 8), intent(out) :: l2norm
real (kind = 8), intent(out) :: ret
real   (kind=8) :: fval,kqval
real   (kind=8) :: dvx,dvy,dvz,rhs,v


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

  l2norm = J*W*v*Uval*Uval
else
  fval = initial_state(X(1),X(2),X(3))
  ret= J*W*v*fval
  l2norm = J*W*v*fval*fval
endif

end subroutine ComputePointForRHS


end module RHS_eq