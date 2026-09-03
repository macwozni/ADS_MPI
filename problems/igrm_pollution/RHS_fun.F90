!------------------------------------------------------------------------------
!> @file RHS_fun.F90
!> @brief Source callback for the corrected 3D pollution DPG problem.
!------------------------------------------------------------------------------
module RHS_fun

   implicit none
   private

   public :: forcing

contains

!------------------------------------------------------------------------------
!> Return the compact emission source; the previous state is intentionally
!> ignored because the source is prescribed in physical space.
!------------------------------------------------------------------------------
function forcing(un, du, X) result(value)
   use input_data, only: emission
   implicit none

   real(kind=8), intent(in) :: un
   real(kind=8), intent(in) :: du(3)
   real(kind=8), intent(in) :: X(3)
   real(kind=8) :: value

   value = emission(X)

end function forcing

end module RHS_fun
