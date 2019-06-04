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
un, &
un13, &
un23, &
ads_data, J, W, direction, substep, l2norm, ret)
use Setup, ONLY: ADS_Setup,ADS_compute_data
use input_data
implicit none
type (ADS_setup), intent(in) :: ads
real   (kind=8), intent(in), dimension(3)  :: X
integer(kind=4), intent(in), dimension(3)  :: k
integer(kind=4), intent(in), dimension(3)  :: e
integer(kind=4), intent(in), dimension(3)  :: a
real   (kind=8), intent(in), dimension(3)  :: du
real (kind = 8), intent(in) :: un,un13,un23
type (ADS_compute_data), intent(in) :: ads_data
real   (kind=8), intent(in)  :: J,W
integer (kind=4), intent(in) :: direction,substep
real (kind = 8), intent(out) :: l2norm
real (kind = 8), intent(out) :: ret
real   (kind=8) :: fval
real   (kind=8) :: v

v   = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3))

fval = 1.d0 !initial_state(X(1),X(2),X(3))
ret= J*W*v*fval
l2norm = J*W*v*fval*fval
end subroutine ComputePointForRHS


end module RHS_eq