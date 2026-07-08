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
!>   used during pointwise right-hand-side evaluation.
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
interface
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

end module Interfaces