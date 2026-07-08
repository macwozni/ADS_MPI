!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data and pointwise forcing for the scalar L2 projection
!> problem.
!>
!> @details
!> This module stores command-line parameters describing the tensor-product
!> grid, spline order, and MPI process grid used by the L2 projection test.
!> It also provides the forcing callback consumed by the current ADS API.
!
!------------------------------------------------------------------------------
module input_data

   implicit none

!> @brief Polynomial order of the approximation space.
   integer(kind = 4) :: order

!> @brief Numbers of elements in the three parametric directions.
   integer(kind = 4) :: isizex, isizey, isizez

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads the L2 projection parameters from command-line arguments.
!>
!> @details
!> The expected argument list is:
!> `<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>`.
!> The values are stored in module variables used by the program driver.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      character(100) :: input
      integer(kind = 4) :: length
      integer(kind = 4) :: status

      ! ./l2 <size> <order> <procx> <procy> <procz>

      if (COMMAND_ARGUMENT_COUNT() .NE. 7) then
         write(*,*) "proper usage with arguments: ", &
         "<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>"
         STOP 5
      end if
      
      
      call GET_COMMAND_ARGUMENT(1, input, length, status)
      read(input, *) isizex
      call GET_COMMAND_ARGUMENT(2, input, length, status)
      read(input, *) isizey
      call GET_COMMAND_ARGUMENT(3, input, length, status)
      read(input, *) isizez
      call GET_COMMAND_ARGUMENT(4, input, length, status)
      read(input, *) order
      call GET_COMMAND_ARGUMENT(5, input, length, status)
      read(input, *) procx
      call GET_COMMAND_ARGUMENT(6, input, length, status)
      read(input, *) procy
      call GET_COMMAND_ARGUMENT(7, input, length, status)
      read(input, *) procz
      
   end subroutine InitializeParameters


!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Returns the constant source term used by the L2 projection test.
!>
!> @details
!> The callback follows the pointwise forcing interface expected by
!> \ref ADSS. The current test projects the constant value one, so the
!> previous solution, its gradient, and the physical point are accepted
!> only to match the interface.
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

end module input_data
