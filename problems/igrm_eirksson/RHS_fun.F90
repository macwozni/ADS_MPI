!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Volume-forcing callback for the stationary 3D Eriksson problem.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none
   private

   public :: forcing

contains

!> @brief Evaluate the manufactured source at one physical point.
!>
!> The previous solution and its gradient are accepted only to satisfy the
!> common ADS forcing callback interface; this stationary linear source does
!> not depend on either value.
   function forcing(un, du, X) result(value)
      use input_data, only: exact_forcing
      implicit none
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: du(3)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value

      value = exact_forcing(X)

   end function forcing

end module RHS_fun
