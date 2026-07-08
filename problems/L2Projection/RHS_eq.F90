!------------------------------------------------------------------------------
!
! MODULE: RHS_eq
!
! DESCRIPTION:
!> @file RHS_eq.F90
!> @brief Legacy quadrature-level RHS callback for the L2 projection
!> problem.
!>
!> @details
!> This module contains the historical \c ComputePointForRHS interface used
!> by older ADS drivers. The current problem driver uses the newer pointwise
!> forcing callback from \ref input_data, but this file is kept documented as
!> part of the legacy problem sources.
!
!------------------------------------------------------------------------------
module RHS_eq

implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the legacy L2 projection right-hand side at one
!> quadrature point.
!>
!> @details
!> The routine multiplies the constant source value by the local test
!> function, quadrature weight, and element Jacobian, and also forms the
!> associated local contribution to the squared L2 norm.
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
real   (kind=8) :: fval
real   (kind=8) :: v

v   = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3))

fval = 1.d0 !initial_state(X(1),X(2),X(3))
ret= J*W*v*fval
l2norm = J*W*v*fval*fval
end subroutine ComputePointForRHS


end module RHS_eq
