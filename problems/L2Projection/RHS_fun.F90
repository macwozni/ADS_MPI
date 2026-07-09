!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise forcing callback for the scalar L2 projection problem.
!>
!> @details
!> This module follows the current ADS API: the problem driver passes
!> \ref forcing to \ref Step, and the common RHS assembly in \ref src/RHS_eq.F90
!> evaluates it at quadrature points.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Returns the constant source term used by the L2 projection test.
!>
!> @details
!> The current test projects the constant value one. The previous solution,
!> its gradient, and the physical point are accepted only to match the ADS
!> pointwise forcing interface.
!
! Input:
! ------
!> @param[in] un
!> Solution value from the previous state at the quadrature point.
!>
!> @param[in] du
!> Gradient of the previous solution at the quadrature point.
!>
!> @param[in] X
!> Physical coordinates of the quadrature point.
!
! Output:
! -------
!> @return ret
!> Constant forcing value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(ret)
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: ret

      ret = 1.d0

   end function forcing

end module RHS_fun
