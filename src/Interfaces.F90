!------------------------------------------------------------------------------
!
! MODULE: Interfaces
!
! DESCRIPTION:
!> @file Interfaces.F90
!> @brief Module defining abstract procedural interfaces used across the
!> project.
!>
!> @details
!> This module collects interface blocks that formalize callback
!> signatures passed between higher-level computational layers.
!>
!> In the current version, the module defines:
!> - \ref forcing_fun, an abstract scalar forcing-function interface
!>   used during pointwise right-hand-side evaluation,
!> - \ref rhs_point_fun, an optional full point-integrand interface used by
!>   problems whose weak form cannot be represented by a scalar source alone.
!>
!> This interface is consumed by modules such as:
!> - \ref RHS_eq,
!> - \ref projection_engine,
!> - \ref ADSS,
!>
!> where user-provided or problem-specific forcing terms are injected
!> into the ADS assembly workflow.
!
!------------------------------------------------------------------------------
module Interfaces

   implicit none

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Abstract interface for scalar forcing functions used in the
!> right-hand-side assembly process.
!>
!> @details
!> This interface specifies the signature of callback procedures that
!> evaluate a scalar forcing term at a given point, using:
!> - the reconstructed scalar solution value \p un,
!> - the reconstructed directional derivative vector \p du,
!> - the physical coordinates \p X of the evaluation point.
!>
!> The returned value is a scalar contribution later used by the ADS
!> pointwise assembly layer.
!
! Input:
! ------
!> @param[in] un
!> Reconstructed scalar solution value at the current quadrature point.
!>
!> @param[in] du
!> Reconstructed directional derivatives of the solution.
!>
!> @param[in] X
!> Physical coordinates of the current quadrature point.
!
! Output:
! -------
!> @return ret
!> Scalar forcing value evaluated at the supplied state and point.
!
! Notes:
! ------
!> @note
!> The interface is abstract and does not provide an implementation.
!> Concrete problem-dependent forcing functions must conform to this
!> signature.
!
!---------------------------------------------------------------------------
abstract interface
   function forcing_fun(un, du, X) result(ret)
      implicit none
!> @brief Reconstructed scalar solution value at the evaluation point.
      real(kind=8), intent(in) :: un
!> @brief Directional derivatives of the solution.
      real(kind=8), intent(in), dimension(3) :: du
!> @brief Physical coordinates of the evaluation point.
      real(kind=8), intent(in), dimension(3) :: X
!> @brief Scalar forcing value returned by the callback.
      real(kind=8) :: ret
   end function

end interface

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Abstract interface for a full pointwise RHS integrand callback.
!>
!> @details
!> This callback receives the same quadrature-point context as the default
!> \ref RHS_eq implementation and returns the already weighted contribution
!> accumulated into the element-local right-hand-side array. It is intended
!> for problem models with spatially varying or nonlinear weak-form terms.
!
!---------------------------------------------------------------------------
abstract interface
   subroutine rhs_point_fun( &
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
      import :: forcing_fun
      implicit none
!> @brief ADS setup structure with basis tables and time-step data.
      type(ADS_setup), intent(in) :: ads
!> @brief Physical coordinates of the current quadrature point.
      real(kind=8), intent(in), dimension(3) :: X
!> @brief Quadrature-point indices in the three directions.
      integer(kind=4), intent(in), dimension(3) :: k
!> @brief Element indices in the three directions.
      integer(kind=4), intent(in), dimension(3) :: e
!> @brief Local basis-function indices in the three directions.
      integer(kind=4), intent(in), dimension(3) :: a
!> @brief Reconstructed directional derivatives of the solution.
      real(kind=8), intent(in), dimension(3) :: du
!> @brief Time-history or step index.
      integer(kind=4), intent(in) :: n
!> @brief Previous-step and intermediate solution values.
      real(kind=8), intent(in) :: un11, un13, un23
!> @brief Runtime ADS data structure.
      type(ADS_compute_data), intent(in) :: ads_data
!> @brief Element Jacobian and quadrature-weight product.
      real(kind=8), intent(in) :: J, W
!> @brief Directional enrichment selector.
      integer(kind=4), intent(in), dimension(3) :: direction
!> @brief Number of the active substep.
      integer(kind=4), intent(in) :: substep
!> @brief Coefficient table for all substeps.
      real(kind=8), intent(in), dimension(7, 3) :: alpha_step
!> @brief Scalar forcing callback retained for source-like terms.
      procedure(forcing_fun) :: forcing
!> @brief Pointwise contribution returned to the assembly routine.
      real(kind=8), intent(out) :: ret
   end subroutine rhs_point_fun

end interface

end module Interfaces
