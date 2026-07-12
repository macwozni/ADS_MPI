!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise callbacks for the iGRM heat example.
!>
!> @details
!> This module provides the scalar forcing callback and a full
!> quadrature-point weak-form callback used by the iGRM multi-step path.
!> The first pseudo-step projects the heat initial condition, while later
!> calls assemble the diffusion-only heat equation terms selected by the
!> configured iGRM time scheme.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the scalar forcing for the iGRM heat example.
!>
!> @details
!> At the initial pseudo-step the callback returns the initial state. For
!> later physical time steps it returns zero, so the evolution is driven
!> only by the heat diffusion operator.
!
! Input:
! ------
!> @param[in] un
!> Previous solution value at the quadrature point.
!>
!> @param[in] du
!> Previous solution gradient at the quadrature point.
!>
!> @param[in] X
!> Physical coordinates of the quadrature point.
!
! Output:
! -------
!> @return ret
!> Pointwise forcing value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(ret)
      use input_data, ONLY: t, initial_state
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: ret

      if (t > 0.d0) then
         ret = 0.d0
      else
         ret = initial_state(X(1), X(2), X(3))
      end if

   end function forcing

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the iGRM heat weak-form contribution at one point.
!>
!> @details
!> For \f$t=0\f$ the callback returns only the mass-projection
!> contribution \f$(u_0,v)\f$. The driver pairs this with a mass-only iGRM
!> setup to initialize the trial-space coefficients.
!>
!> For \f$t>0\f$ the callback follows the configured iGRM time-scheme
!> tables. It keeps the mass term, the directional diffusion rows, and
!> the scalar forcing row:
!> \f[
!>   (u,v) + \tau\sum_i \alpha_i(\partial_i u, \partial_i v)
!>          + \tau\alpha_f(f,v).
!> \f]
!> Transport rows are intentionally omitted, matching the heat equation.
!
! Input:
! ------
!> @param[in] ads
!> Mixed iGRM setup structure with basis and quadrature data.
!>
!> @param[in] X
!> Physical coordinates of the quadrature point.
!>
!> @param[in] k
!> Quadrature indices in the three directions.
!>
!> @param[in] e
!> Element indices in the three directions.
!>
!> @param[in] a
!> Local basis-function indices in the three directions.
!>
!> @param[in] du
!> Reconstructed gradient of the currently selected state.
!>
!> @param[in] un11
!> Previous-step solution value at the quadrature point.
!>
!> @param[in] un13
!> First intermediate solution value.
!>
!> @param[in] un23
!> Second intermediate solution value.
!>
!> @param[in] ads_data
!> Runtime ADS data containing derivative-state buffers and selectors.
!>
!> @param[in] J
!> Element Jacobian factor.
!>
!> @param[in] W
!> Product of quadrature weights.
!>
!> @param[in] substep
!> Current iGRM substep number.
!>
!> @param[in] alpha_step
!> Coefficient table of the selected iGRM time scheme.
!>
!> @param[in] forcing_cb
!> Scalar forcing callback.
!
! Output:
! -------
!> @param[out] ret
!> Weighted pointwise contribution accumulated into the RHS.
!
!---------------------------------------------------------------------------
   subroutine heat_igrm_rhs_point( &
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
      forcing_cb, &
      ret)
      use Setup, ONLY: ADS_Setup, ADS_compute_data
      use Interfaces, ONLY: forcing_fun
      use input_data, ONLY: t
      implicit none
      type(ADS_setup), intent(in) :: ads
      real(kind=8), intent(in), dimension(3) :: X
      integer(kind=4), intent(in), dimension(3) :: k
      integer(kind=4), intent(in), dimension(3) :: e
      integer(kind=4), intent(in), dimension(3) :: a
      real(kind=8), intent(in), dimension(3) :: du
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: un11, un13, un23
      type(ADS_compute_data), intent(in) :: ads_data
      real(kind=8), intent(in) :: J, W
      integer(kind=4), intent(in), dimension(3) :: direction
      integer(kind=4), intent(in) :: substep
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
      procedure(forcing_fun) :: forcing_cb
      real(kind=8), intent(out) :: ret
      real(kind=8), dimension(3) :: dv
      real(kind=8), dimension(6) :: rhs_du
      real(kind=8), dimension(7) :: alpha
      real(kind=8) :: v, u, source
      integer(kind=4) :: idir, statex, statey, statez, term

      v = ads%NNx(0, a(1), k(1), e(1))* &
          ads%NNy(0, a(2), k(2), e(2))* &
          ads%NNz(0, a(3), k(3), e(3))
      dv(1) = ads%NNx(1, a(1), k(1), e(1))* &
              ads%NNy(0, a(2), k(2), e(2))* &
              ads%NNz(0, a(3), k(3), e(3))
      dv(2) = ads%NNx(0, a(1), k(1), e(1))* &
              ads%NNy(1, a(2), k(2), e(2))* &
              ads%NNz(0, a(3), k(3), e(3))
      dv(3) = ads%NNx(0, a(1), k(1), e(1))* &
              ads%NNy(0, a(2), k(2), e(2))* &
              ads%NNz(1, a(3), k(3), e(3))

      source = forcing_cb(un11, du, X)
      if (t <= 0.d0) then
         ret = J*W*v*source
         return
      end if

      alpha = alpha_step(:, substep)
      select case (substep)
      case (1)
         u = un11
      case (2)
         u = un13
      case (3)
         u = un23
      case default
         stop "wrong substep"
      end select

      statex = e(1) - ads_data%state_mine(1) + 1
      statey = e(2) - ads_data%state_mine(2) + 1
      statez = e(3) - ads_data%state_mine(3) + 1
      do term = 1, 6
         idir = (term + 1)/2
         select case (ads_data%rhs_du_state(term, substep))
         case (0)
            rhs_du(term) = du(idir)
         case (1)
            rhs_du(term) = ads_data%dUn0(statex, statey, statez, k(1), k(2), k(3), idir)
         case (2)
            rhs_du(term) = ads_data%dUn13(statex, statey, statez, k(1), k(2), k(3), idir)
         case (3)
            rhs_du(term) = ads_data%dUn23(statex, statey, statez, k(1), k(2), k(3), idir)
         case default
            stop "wrong RHS derivative state"
         end select
      end do

      source = forcing_cb(u, du, X)
      ret = J*W*(u*v + ads%tau*( &
         alpha(1)*rhs_du(1)*dv(1) + &
         alpha(3)*rhs_du(3)*dv(2) + &
         alpha(5)*rhs_du(5)*dv(3) + &
         alpha(7)*source*v))

   end subroutine heat_igrm_rhs_point

end module RHS_fun
