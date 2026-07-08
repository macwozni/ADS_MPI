!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input parameters for the iGRM L2 projection example.
!>
!> @details
!> The module stores the command-line grid dimensions, polynomial order,
!> and MPI process-grid dimensions used by the iGRM L2 projection driver.
!> The actual forcing callback is defined in \ref RHS_fun for this example.
!
!------------------------------------------------------------------------------
module input_data

   implicit none

!> @brief Polynomial order requested from the command line.
   integer(kind = 4) :: order

!> @brief Numbers of elements in the three parametric directions.
   integer(kind = 4) :: isizex, isizey, isizez

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads iGRM L2 projection parameters from command-line arguments.
!>
!> @details
!> The expected argument list is:
!> `<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>`.
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



end module input_data
