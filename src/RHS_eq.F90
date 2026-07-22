!------------------------------------------------------------------------------
!
! MODULE: RHS_eq
!
! DESCRIPTION:
!> @file RHS_eq.F90
!> @brief Module providing pointwise right-hand-side evaluation for the
!> ADS assembly workflow.
!>
!> @details
!> This module contains the local integrand construction used during
!> right-hand-side assembly in the ADS scheme.
!>
!> The provided functionality includes:
!> - evaluation of a test-function value and its directional derivatives
!>   at one quadrature point,
!> - incorporation of previously reconstructed solution values and
!>   derivative data,
!> - inclusion of an external forcing term supplied through a callback,
!> - formation of the pointwise contribution accumulated later by the
!>   assembly layer.
!>
!> Within the project, this module is called from
!> \ref projection_engine via \ref ComputePointForRHS and therefore sits
!> on the critical path between:
!> - quadrature-point reconstruction of state variables,
!> - pointwise PDE residual evaluation,
!> - elementwise accumulation into ADS right-hand-side buffers.
!
!------------------------------------------------------------------------------
module RHS_eq

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Computes the pointwise contribution to the assembled
!> right-hand side at one quadrature point.
!>
!> @details
!> This routine evaluates the contribution associated with a single
!> quadrature point and a single local tensor-product basis function.
!>
!> The procedure:
!> - selects the active time-substep coefficient vector from
!>   \p alpha_step,
!> - chooses the appropriate intermediate solution state depending on
!>   \p substep,
!> - evaluates the tensor-product test function and its directional
!>   derivatives at the specified quadrature point,
!> - evaluates the forcing term through the callback \p forcing,
!> - forms a weighted combination of previous-solution values, derivative
!>   contributions, forcing, and time-step scaling.
!>
!> The test-function value is denoted by
!>
!> \f[
!> v = N_x N_y N_z,
!> \f]
!>
!> while the directional derivatives are assembled as
!>
!> \f[
!> \partial_x v,\quad \partial_y v,\quad \partial_z v.
!> \f]
!>
!> These quantities are then combined with the derivative vector
!> \p du, the forcing value, the time-step length \p ads%tau, and the
!> selected intermediate solution value.
!
! Input:
! ------
!> @param[in] ads
!> ADS setup structure containing basis values, derivative tables, and
!> time-step information.
!>
!> @param[in] X
!> Physical coordinates of the current quadrature point.
!>
!> @param[in] k
!> Quadrature-point indices in the three coordinate directions.
!>
!> @param[in] e
!> Element indices in the three coordinate directions.
!>
!> @param[in] a
!> Local basis-function indices in the three coordinate directions.
!>
!> @param[in] du
!> Reconstructed directional derivatives of the previous solution at the
!> current quadrature point.
!>
!> @param[in] n
!> Time-history or step index passed through the assembly pipeline.
!>
!> @param[in] un11
!> Reference previous-step solution value at the current quadrature
!> point.
!>
!> @param[in] un13
!> Intermediate solution value associated with the first fractional
!> stage.
!>
!> @param[in] un23
!> Intermediate solution value associated with the second fractional
!> stage.
!>
!> @param[in] ads_data
!> Runtime ADS data structure passed through the interface.
!>
!> @param[in] J
!> Jacobian factor of the current element mapping.
!>
!> @param[in] W
!> Product of quadrature weights associated with the current quadrature
!> point.
!>
!> @param[in] direction
!> Directional enrichment selector of the current substep.
!>
!> @param[in] substep
!> Number of the current substep.
!>
!> @param[in] alpha_step
!> Table of coefficients used by all substeps; the active column is
!> selected according to \p substep.
!>
!> @param[in] forcing
!> Callback returning the forcing value for the current state and point.
!
! Output:
! -------
!> @param[out] ret
!> Pointwise contribution to be accumulated into the assembled
!> right-hand side.
!
! Notes:
! ------
!> @note
!> The active coefficient vector is extracted as `alpha_step(:,substep)`.
!
!> @note
!> The arguments \p n, \p ads_data, \p J, \p W, and \p direction are
!> presently carried through the interface for consistency with the
!> surrounding workflow, although not all of them are explicitly used in
!> the current implementation.
!
!> @warning
!> The procedure assumes that \p substep belongs to the set
!> \f$\{1,2,3\}\f$.
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
   un11, &
   un13, &
   un23, &
   ads_data, J, W, direction, substep, &
   alpha_step, &
   forcing, &
   ret)
   use Setup, ONLY: ADS_Setup, ADS_compute_data
   use Interfaces, ONLY: forcing_fun
   use ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
   implicit none
!> @brief ADS setup structure with basis tables and time-step data.
   type(ADS_setup), intent(in) :: ads
!> @brief Physical coordinates of the current quadrature point.
   real(kind=8), intent(in), dimension(3)  :: X
!> @brief Quadrature-point indices in the three directions.
   integer(kind=4), intent(in), dimension(3)  :: k
!> @brief Element indices in the three directions.
   integer(kind=4), intent(in), dimension(3)  :: e
!> @brief Local basis-function indices in the three directions.
   integer(kind=4), intent(in), dimension(3)  :: a
!> @brief Directional derivatives of the reconstructed solution.
   real(kind=8), intent(in), dimension(3)  :: du
!> @brief Time-history or step index.
   integer(kind=4), intent(in) :: n
!> @brief Previous-step and intermediate solution values.
   real(kind=8), intent(in) :: un11, un13, un23
!> @brief Runtime ADS data structure passed through the interface.
   type(ADS_compute_data), intent(in) :: ads_data
!> @brief Element Jacobian and quadrature-weight product.
   real(kind=8), intent(in)  :: J, W
!> @brief Directional enrichment selector.
   integer(kind=4), dimension(3), intent(in) :: direction
!> @brief Number of the active substep.
   integer(kind=4), intent(in) :: substep
!> @brief Forcing callback used in pointwise evaluation.
   procedure(forcing_fun) :: forcing
!> @brief Coefficient table for all substeps.
   real(kind=8), intent(in), dimension(7, 3) :: alpha_step
!> @brief Pointwise contribution returned to the assembly routine.
   real(kind=8), intent(out) :: ret
!> @brief Temporary accumulator for the pointwise value.
   real(kind=8) :: fval
!> @brief Directional derivatives of the local test function.
   real(kind=8), dimension(3) :: dv
!> @brief Forcing value, test-function value, and selected state value.
   real(kind=8) :: rhs, v, u
!> @brief Active coefficient vector extracted for the current substep.
   real(kind=8), dimension(7) :: alpha
!> @brief Directional derivatives selected independently for RHS terms.
   real(kind=8), dimension(6) :: rhs_du
!> @brief Local indices into the stored derivative-state buffers.
   integer(kind=4) :: statex, statey, statez
!> @brief Loop and direction selectors for derivative-state lookup.
   integer(kind=4) :: term, idir

   select case (substep)
   case (1)
      u = un11
   case (2)
      u = un13
   case (3)
      u = un23
   case default
      write (ERROR_UNIT, '(A,I0)') 'wrong substep: ', substep
      flush (ERROR_UNIT)
      STOP 1
   end select
   alpha = alpha_step(:, substep)

   v = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
   dv(1) = ads%NNx(1, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
   dv(2) = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(1, a(2), k(2), e(2))*ads%NNz(0, a(3), k(3), e(3))
   dv(3) = ads%NNx(0, a(1), k(1), e(1))*ads%NNy(0, a(2), k(2), e(2))*ads%NNz(1, a(3), k(3), e(3))

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
         write (ERROR_UNIT, '(A,I0,A,I0,A,I0)') &
            'wrong RHS derivative state: ', ads_data%rhs_du_state(term, substep), &
            ', term: ', term, ', substep: ', substep
         flush (ERROR_UNIT)
         STOP 1
      end select
   end do

   rhs = forcing(u, du, X)

   fval = u*v
   fval = fval + ads%tau*( &
      alpha(1)*rhs_du(1)*dv(1) + alpha(2)*rhs_du(2)*v + &
      alpha(3)*rhs_du(3)*dv(2) + alpha(4)*rhs_du(4)*v + &
      alpha(5)*rhs_du(5)*dv(3) + alpha(6)*rhs_du(6)*v + &
      alpha(7)*rhs*v)

   ret = J*W*fval

end subroutine ComputePointForRHS

end module RHS_eq
