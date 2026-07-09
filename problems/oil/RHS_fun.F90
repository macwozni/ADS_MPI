!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise forcing callback for the oil transport example.
!>
!> @details
!> This module provides callbacks passed to the current ADS \ref Step API.
!> The scalar \ref forcing callback is retained for source-like terms, while
!> \ref oil_rhs_point restores the oil problem's full nonlinear weak form.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pointwise source term for the current ADS API.
!>
!> @details
!> At the initial pseudo-step the callback returns the initial state. For
!> later time steps it returns pump injection minus the nonnegative drain
!> contribution.
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
!> Pointwise source value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(ret)
      use input_data, ONLY: t, pumping, draining, initial_state
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: ret

      if (t > 0.d0) then
         ret = pumping(X(1), X(2), X(3)) - &
               max(0.d0, draining(un, X(1), X(2), X(3)))
      else
         ret = initial_state(X(1), X(2), X(3))
      endif

   end function forcing

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the full pointwise oil weak-form contribution.
!>
!> @details
!> The oil model cannot be represented by the scalar source callback alone:
!> its diffusion term depends on the cached permeability, the reconstructed
!> solution, and the gradient of the test function. This callback therefore
!> restores the legacy contribution
!> \f$-K_q \exp(\mu u)\nabla u\cdot\nabla v\f$ and updates the drained
!> diagnostic at the same quadrature point.
!
! Input:
! ------
!> @param[in] ads
!> ADS setup structure with basis and quadrature data.
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
!> Reconstructed gradient of the previous solution.
!>
!> @param[in] un11
!> Previous-step solution value at the quadrature point.
!
! Output:
! -------
!> @param[out] ret
!> Weighted pointwise contribution accumulated into the RHS.
!
!---------------------------------------------------------------------------
   subroutine oil_rhs_point( &
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
      use input_data, ONLY: t, mi, Kqvals, drained, initial_state, draining
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
      real(kind=8) :: v, dvx, dvy, dvz
      real(kind=8) :: source, vdrain, kqval, grad_term

      v = ads%NNx(0, a(1), k(1), e(1))* &
          ads%NNy(0, a(2), k(2), e(2))* &
          ads%NNz(0, a(3), k(3), e(3))
      dvx = ads%NNx(1, a(1), k(1), e(1))* &
            ads%NNy(0, a(2), k(2), e(2))* &
            ads%NNz(0, a(3), k(3), e(3))
      dvy = ads%NNx(0, a(1), k(1), e(1))* &
            ads%NNy(1, a(2), k(2), e(2))* &
            ads%NNz(0, a(3), k(3), e(3))
      dvz = ads%NNx(0, a(1), k(1), e(1))* &
            ads%NNy(0, a(2), k(2), e(2))* &
            ads%NNz(1, a(3), k(3), e(3))

      if (t > 0.d0) then
         kqval = Kqvals(k(1), k(2), k(3), &
                        e(1) - ads%mine(1) + 1, &
                        e(2) - ads%mine(2) + 1, &
                        e(3) - ads%mine(3) + 1)
         source = forcing_cb(un11, du, X)
         vdrain = max(0.d0, draining(un11, X(1), X(2), X(3)))
         grad_term = du(1)*dvx + du(2)*dvy + du(3)*dvz
         ret = J*W*(v*un11 + ads%tau*(-kqval*exp(mi*un11)*grad_term + v*source))
         drained = drained + J*W*v*ads%tau*vdrain
      else
         ret = J*W*v*initial_state(X(1), X(2), X(3))
      end if

   end subroutine oil_rhs_point

end module RHS_fun
