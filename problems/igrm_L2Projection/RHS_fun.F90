!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise forcing callback for the iGRM L2 projection example.
!>
!> @details
!> This module provides the forcing function passed to \ref MultiStep.
!> The commented block above the active function preserves the older
!> quadrature-level RHS interface for reference.
!
!------------------------------------------------------------------------------
module RHS_fun

implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pointwise forcing used by the iGRM L2 projection.
!>
!> @details
!> The callback returns the product of the physical coordinates multiplied
!> by a constant source value. The previous solution and gradient are
!> accepted to satisfy the current ADS forcing interface.
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
function forcing (un,du,X) result (ret)
    implicit none
    real (kind = 8), intent(in) :: un
    real (kind = 8), intent(in), dimension(3) :: du
    real (kind = 8), intent(in), dimension(3) :: X
    real (kind = 8) :: ret

    real (kind=8) :: fval
    real (kind=8) :: v

    fval = 1.d0 !initial_state(X(1),X(2),X(3))
    v = X(1)*X(2)*X(3)
    ret = v*fval
 end function



end module RHS_fun
