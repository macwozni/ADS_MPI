!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise forcing callback for the pure-diffusion iGRM example.
!>
!> @details
!> This module provides the callback passed to the current ADS \ref Step API.
!> The migrated pure-diffusion example currently uses zero volume forcing.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Evaluates the pure-diffusion volume forcing.
!>
!> @details
!> Arguments are accepted to match the current ADS pointwise forcing
!> interface.
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
!> @return fval
!> Volume forcing value.
!
!---------------------------------------------------------------------------
   function forcing(un, du, X) result(fval)
      implicit none
      real(kind = 8), intent(in) :: un
      real(kind = 8), intent(in), dimension(3) :: du
      real(kind = 8), intent(in), dimension(3) :: X
      real(kind = 8) :: fval

      fval = 0.d0

   end function forcing

end module RHS_fun
