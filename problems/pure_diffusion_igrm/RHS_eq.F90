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
! n               - nuber of previous time steps
! Un              - U_n, previous solution coefficient at given point
! Un13            - U_n+1/3
! Un23            - U_n+2/3
! ads_data        - data structures for ADS
! J               - jacobian
! W               - weight for quadratures
! directon        -
! substep         - 
!
! Output:
! -------
! ret             - value of RHS function at given point
! l2norm          -
!
! -------------------------------------------------------------------

subroutine ComputePointForRHS( &
ads, &
X, &
k, &
e, &
a, &
du, &
n, &
un, &
un13, &
un23, &
ads_data, J, W, direction, substep, l2norm, ret)
use Setup, ONLY: ADS_Setup,ADS_compute_data
use projection_engine, ONLY: global2local
use input_data
implicit none
type (ADS_setup), intent(in) :: ads
real   (kind=8), intent(in), dimension(3)  :: X
integer(kind=4), intent(in), dimension(3)  :: k
integer(kind=4), intent(in), dimension(3)  :: e
integer(kind=4), intent(in), dimension(3)  :: a
real   (kind=8), intent(in), dimension(3)  :: du
integer (kind = 4), intent(in) :: n
real (kind = 8), intent(in), dimension(n)  :: un
real (kind = 8), intent(in) :: un13,un23
type (ADS_compute_data), intent(in) :: ads_data
real   (kind=8), intent(in)  :: J,W
integer (kind=4), intent(in) :: direction,substep
real (kind = 8), intent(out) :: l2norm
real (kind = 8), intent(out) :: ret
real   (kind=8) :: vforce
real   (kind=8) :: Umax = -1d10, Umin = 1d10
real   (kind=8) :: dvx,dvy,dvz,rhs,v
logical :: boundary



boundary = .FALSE.

vforce = forcing(X(1),X(2),X(3))    

v   = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3))
dvx = ads % NNx(1,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvy = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(1,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvz = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(1,a(3),k(3),e(3)) 

if (X(1) .EQ. ads%Ux(0)) boundary = .TRUE
if (X(1) .EQ. ads%Ux(size(ads%Ux))) boundary = .TRUE
if (X(2) .EQ. ads%Uy(0)) boundary = .TRUE
if (X(2) .EQ. ads%Uy(size(ads%Uy))) boundary = .TRUE
if (X(3) .EQ. ads%Uz(0)) boundary = .TRUE
if (X(3) .EQ. ads%Uz(size(ads%Uz))) boundary = .TRUE

ret = vforce * v

if (boundary) then
    bound_norm = ! TODO
    bound = bound_norm * g(X(1),X(2),X(3)) * v
    ret = ret + bound
endif


end subroutine ComputePointForRHS


end module RHS_eq