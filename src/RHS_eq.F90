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
      un11, &
      un13, &
      un23, &
      ads_data, J, W, direction, substep, &
      alpha_step, &
      forcing, &
      l2norm, ret)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use Interfaces, ONLY: forcing_fun
      implicit none
      type(ADS_setup), intent(in) :: ads
      real(kind=8), intent(in), dimension(3)  :: X
      integer(kind=4), intent(in), dimension(3)  :: k
      integer(kind=4), intent(in), dimension(3)  :: e
      integer(kind=4), intent(in), dimension(3)  :: a
      real(kind=8), intent(in), dimension(3)  :: du
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: un11, un13, un23
      type(ADS_compute_data), intent(in) :: ads_data
      real(kind=8), intent(in)  :: J, W
      integer(kind=4), dimension(3), intent(in) :: direction
      integer(kind=4), intent(in) :: substep
      procedure(forcing_fun) :: forcing
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
      real(kind=8), intent(out) :: l2norm
      real(kind=8), intent(out) :: ret
      real(kind=8) :: fval
      real(kind=8), dimension(3) :: dv
      real(kind=8) :: rhs, v, u
      real(kind=8), dimension(7) :: alpha

      alpha = alpha_step(:, substep)
      if (substep .EQ. 1) u = un11
      if (substep .EQ. 2) u = un13
      if (substep .EQ. 3) u = un23

      v = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
      dv(1) = ads%NNx(1, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
      dv(2) = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(1, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
      dv(3) = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(1, a(3), k(3), e(3))

      rhs = forcing(un11, du, X)

      fval = un11*v
      fval = fval + alpha(1)*du(1)*dv(1)
      fval = fval + alpha(2)*du(1)*v
      fval = fval + alpha(3)*du(2)*dv(2)
      fval = fval + alpha(4)*du(2)*v
      fval = fval + alpha(5)*du(3)*dv(3)
      fval = fval + alpha(6)*du(3)*v
      fval = fval*ads%tau
      fval = fval + alpha(3)*rhs*v
      fval = fval + u*v
   end subroutine ComputePointForRHS

end module RHS_eq
