!------------------------------------------------------------------------------
!
! MODULE: RHS_eq
!
! DESCRIPTION:
!> @file RHS_eq.F90
!> @brief Legacy quadrature-level RHS callback for the pure-diffusion iGRM
!> example.
!>
!> @details
!> This module preserves the older boundary-aware RHS interface. The
!> current migrated driver uses the pointwise forcing callback from
!> \ref input_data.
!
!------------------------------------------------------------------------------
module RHS_eq

implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Forms the legacy local RHS contribution with boundary terms.
!>
!> @details
!> The routine was intended to add volume forcing and optional boundary
!> influx contributions at quadrature points located on the domain
!> boundary. It is retained as legacy source alongside the migrated driver.
!
! Input:
! ------
!> @param[in] ads
!> ADS setup structure containing basis values and boundary knot data.
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
real (kind = 9) :: epsilon_zero


boundary = .FALSE.
epsilon_zero = 10E-14

vforce = forcing(X(1),X(2),X(3))    

v   = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3))
dvx = ads % NNx(1,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvy = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(1,a(2),k(2),e(2)) * ads % NNz(0,a(3),k(3),e(3)) 
dvz = ads % NNx(0,a(1),k(1),e(1)) * ads % NNy(0,a(2),k(2),e(2)) * ads % NNz(1,a(3),k(3),e(3)) 

if (abs((X(1) - ads%Ux(0)) boundary) .LE. epsilon_zero) = .TRUE
if (abs((X(1) - ads%Ux(size(ads%Ux)))) .LE. epsilon_zero) boundary = .TRUE
if (abs((X(2) - ads%Uy(0))) .LE. epsilon_zero) boundary = .TRUE
if (abs((X(2) - ads%Uy(size(ads%Uy)))) .LE. epsilon_zero) boundary = .TRUE
if (abs((X(3) - ads%Uz(0))) .LE. epsilon_zero) boundary = .TRUE
if (abs((X(3) - ads%Uz(size(ads%Uz)))) .LE. epsilon_zero) boundary = .TRUE

ret = vforce * v

if (boundary) then
! influx
    bound_norm = abs(b(X(1),X(2),X(3)) * n(X(1),X(2),X(3)))
    bound_norm = bound_norm - b(X(1),X(2),X(3) * n(X(1),X(2),X(3))
    bound_norm = bound_norm / 2.d0
! g - source on boundary
    bound = bound_norm * g(X(1),X(2),X(3)) * v
    ret = ret + bound
endif


end subroutine ComputePointForRHS


end module RHS_eq
