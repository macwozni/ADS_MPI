module RHS_eq

implicit none

contains

! -------------------------------------------------------------------
! Right-hand side of the equation.
!
! Input:
! ------
! X_              - quadrature points
! k_              - inde(1)es for quadrature points
! e_              - inde(1)es for elements
! p_              - degrees of approximation
! nelem_          - number of subintervals
! a_, b_          - inde(1)es of basis functions
! NN_             - values of basis functions in quadrature points
! du_             - value of derivative from previous time step
! ibeg_, iend_    - piece of domain associated with this process
! ibegs_, iends_  - pieces of domain surrounding this process' piece
! mine_, ma(1)e_    - indices of first and last elements in each direction
! Uval            - previous solution coefficient at given point
! J               - jacobian
! W               - weight for quadratures
!
! Output:
! -------
! F               - value of RHS function at given point
!
! -------------------------------------------------------------------

subroutine ComputePointForRHS( &
ads, &
X, &
k, &
e, &
a, &
b, &
du, &
NNx, NNy, NNz, &
Uval, J, W, ret)
use Setup, ONLY: ADS_Setup
use input_data
implicit none
type (ADS_setup), intent(in) :: ads
real   (kind=8), intent(in), dimension(3)  :: X
integer(kind=4), intent(in), dimension(3)  :: k
integer(kind=4), intent(in), dimension(3)  :: e
integer(kind=4), intent(in), dimension(3)  :: a
integer(kind=4), intent(in), dimension(3)  :: b
real   (kind=8), intent(in), dimension(3)  :: du
real   (kind=8), intent(in)  :: Uval
real   (kind=8), intent(in)  :: J,W
real (kind = 8), intent(in) :: &
NNx(0:1, 0:ads%p(1), ads%p(1) + 1, ads%nelem(1)), &
NNy(0:1, 0:ads%p(2), ads%p(2) + 1, ads%nelem(2)), &
NNz(0:1, 0:ads%p(3), ads%p(3) + 1, ads%nelem(3))
real (kind = 8), intent(out) :: ret
real   (kind=8) :: fval,vpump,vdrain,kqval
real   (kind=8) :: Umax = -1d10, Umin = 1d10
real   (kind=8) :: dvx,dvy,dvz,rhs,v


vpump = pumping(X(1),X(2),X(3))    
Umax = max(Umax, Uval)
Umin = min(Umin, Uval)

v   = NNx(0,a(1),k(1),e(1)) * NNy(0,a(2),k(2),e(2)) * NNz(0,a(3),k(3),e(3))
dvx = NNx(1,a(1),k(1),e(1)) * NNy(0,a(2),k(2),e(2)) * NNz(0,a(3),k(3),e(3)) 
dvy = NNx(0,a(1),k(1),e(1)) * NNy(1,a(2),k(2),e(2)) * NNz(0,a(3),k(3),e(3)) 
dvz = NNx(0,a(1),k(1),e(1)) * NNy(0,a(2),k(2),e(2)) * NNz(1,a(3),k(3),e(3)) 

kqval = Kqvals(k(1),k(2),k(3),e(1)-ads%mine(1)+1,e(2)-ads%mine(2)+1,e(3)-ads%mine(3)+1)
vdrain = max(0.d0, draining(Uval, X(1),X(2),X(3)))
fval = vpump - vdrain
!--- Real
if (t > 0.0) then
  rhs = Dt * ( - kqval * exp(mi * Uval) * (du(1)*dvx + du(2)*dvy + du(3)*dvz) + v * fval)
  ret = J*W*(v * Uval + rhs)

  drained = drained + J*W*v*Dt*vdrain
  l2norm = l2norm + J*W*v*Uval*Uval
else
  fval = initial_state(X(1),X(2),X(3))
  ret= J*W*v*fval
  l2norm = l2norm + J*W*v*fval*fval
endif

end subroutine


end module