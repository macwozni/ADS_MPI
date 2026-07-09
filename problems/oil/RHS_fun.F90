!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Pointwise forcing callback for the oil transport example.
!>
!> @details
!> This module provides the callback passed to the current ADS \ref Step API.
!> Problem-specific pump, drain, and initial-state functions are stored in
!> \ref input_data.
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

end module RHS_fun
