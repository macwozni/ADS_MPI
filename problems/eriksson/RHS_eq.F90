!------------------------------------------------------------------------------
!
! MODULE: RHS_eq
!
! DESCRIPTION:
!> @file RHS_eq.F90
!> @brief Legacy quadrature-level RHS callback for the Eriksson example.
!>
!> @details
!> This module keeps the historical ADS RHS interface used by older
!> drivers. The current migrated driver uses the pointwise forcing callback
!> from \ref input_data.
!
!------------------------------------------------------------------------------
module RHS_eq

implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Forms the legacy local RHS contribution for one quadrature point.
!>
!> @details
!> On the initial step the routine contributes the initial state. On later
!> steps it combines the previous solution value with a simple diffusion
!> term evaluated against derivatives of the local test function.
!
! Input:
! ------
!> @param[in] ads
!> ADS setup structure containing basis values.
!>
!> @param[in] X
!> Physical quadrature-point coordinates.
!>
!> @param[in] k
!> Quadrature-point indices in the three directions.
!>
!> @param[in] e
!> Element indices in the three directions.
!>
!> @param[in] a
!> Local basis-function indices in the three directions.
!>
!> @param[in] du
!> Previous-solution gradient at the quadrature point.
!>
!> @param[in] n
!> Number of stored previous states.
!>
!> @param[in] un
!> Previous solution values at the quadrature point.
!>
!> @param[in] un13
!> Intermediate one-third state value.
!>
!> @param[in] un23
!> Intermediate two-third state value.
!>
!> @param[in] ads_data
!> ADS runtime-data structure.
!>
!> @param[in] J
!> Element Jacobian factor.
!>
!> @param[in] W
!> Quadrature weight.
!>
!> @param[in] direction
!> Directional substep selector.
!>
!> @param[in] substep
!> Substep number in the legacy ADS workflow.
!
! Output:
! -------
!> @param[out] l2norm
!> Local contribution to the squared L2 norm.
!>
!> @param[out] ret
!> Local right-hand-side contribution.
!
!---------------------------------------------------------------------------

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
use input_data
use projection_engine, ONLY: global2local
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
  ret = J*W*(v * un(1) + rhs)

  l2norm = J*W*v*un(1)*un(1)
else
  fval = initial_state(X(1),X(2),X(3))
  ret= J*W*v*fval
  l2norm = J*W*v*fval*fval
endif

end subroutine ComputePointForRHS


end module RHS_eq
