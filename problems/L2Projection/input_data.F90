!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input data for the scalar L2 projection problem.
!>
!> @details
!> This module stores command-line parameters describing the tensor-product
!> grid, spline order, and MPI process grid used by the L2 projection test.
!> The forcing callback consumed by the current ADS API is defined in
!> \ref RHS_fun.
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

      ! ./l2 <size> <order> <procx> <procy> <procz>

      if (COMMAND_ARGUMENT_COUNT() .NE. 7) then
         write(*,*) "proper usage with arguments: ", &
         "<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>"
         STOP 5
      end if
      
      
      call ReadIntegerArgument(1, isizex)
      call ReadIntegerArgument(2, isizey)
      call ReadIntegerArgument(3, isizez)
      call ReadIntegerArgument(4, order)
      call ReadIntegerArgument(5, procx)
      call ReadIntegerArgument(6, procy)
      call ReadIntegerArgument(7, procz)

   contains

      subroutine ReadArgument(arg, input)
         implicit none
         integer(kind = 4), intent(in) :: arg
         character(*), intent(out) :: input
         integer(kind = 4) :: length
         integer(kind = 4) :: status

         call GET_COMMAND_ARGUMENT(arg, input, length, status)
         if (status /= 0) then
            write(*,*) "invalid command argument: ", arg
            STOP 5
         end if
      end subroutine ReadArgument

      subroutine ReadIntegerArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         integer(kind = 4), intent(out) :: value
         character(100) :: input
         integer(kind = 4) :: read_status

         call ReadArgument(arg, input)
         read(input, *, iostat = read_status) value
         if (read_status /= 0) then
            write(*,*) "invalid integer argument: ", arg
            STOP 5
         end if
      end subroutine ReadIntegerArgument
      
   end subroutine InitializeParameters

end module input_data
