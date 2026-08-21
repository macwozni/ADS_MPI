!------------------------------------------------------------------------------
!
! MODULE: RHS_fun
!
! DESCRIPTION:
!> @file RHS_fun.F90
!> @brief Component forcing callbacks for the stationary 3D Stokes problem.
!
!------------------------------------------------------------------------------
module RHS_fun

   implicit none
   private

   public :: forcing_x
   public :: forcing_y
   public :: forcing_z

contains

!> @brief Evaluate the x component of the manufactured body force.
   function forcing_x(un, du, X) result(value)
      use input_data, only: body_force
      implicit none
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: du(3)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: force(3)

      force = body_force(X)
      value = force(1)

   end function forcing_x

!> @brief Evaluate the y component of the manufactured body force.
   function forcing_y(un, du, X) result(value)
      use input_data, only: body_force
      implicit none
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: du(3)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: force(3)

      force = body_force(X)
      value = force(2)

   end function forcing_y

!> @brief Evaluate the z component of the manufactured body force.
   function forcing_z(un, du, X) result(value)
      use input_data, only: body_force
      implicit none
      real(kind=8), intent(in) :: un
      real(kind=8), intent(in) :: du(3)
      real(kind=8), intent(in) :: X(3)
      real(kind=8) :: value
      real(kind=8) :: force(3)

      force = body_force(X)
      value = force(3)

   end function forcing_z

end module RHS_fun
