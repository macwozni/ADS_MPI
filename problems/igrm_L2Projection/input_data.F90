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
!> MPI process-grid dimensions, and selected iGRM time scheme used by the
!> iGRM L2 projection driver. The actual forcing callback is defined in
!> \ref RHS_fun for this example.
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

!> @brief Time-step length used to configure the selected iGRM scheme.
   real(kind = 8) :: scheme_tau = 1.d0

!> @brief Time scheme selected for the iGRM multi-step driver.
   character(len = 32) :: time_scheme = "dg"


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads iGRM L2 projection parameters from command-line arguments.
!>
!> @details
!> The accepted argument lists are:
!> `<isizex> <isizey> <isizez> <order> <procx> <procy> <procz>`,
!> `<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> <scheme>`,
!> or
!> `<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> <tau> <scheme>`.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      integer(kind = 4) :: nargs

      nargs = COMMAND_ARGUMENT_COUNT()

      if (nargs .NE. 7 .AND. nargs .NE. 8 .AND. nargs .NE. 9) then
         write(*,*) "proper usage with arguments: ", &
         "<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> [scheme]"
         write(*,*) "or: ", &
         "<isizex> <isizey> <isizez> <order> <procx> <procy> <procz> <tau> <scheme>"
         STOP 5
      end if

      call ReadIntegerArgument(1, isizex)
      call ReadIntegerArgument(2, isizey)
      call ReadIntegerArgument(3, isizez)
      call ReadIntegerArgument(4, order)
      call ReadIntegerArgument(5, procx)
      call ReadIntegerArgument(6, procy)
      call ReadIntegerArgument(7, procz)
      scheme_tau = 1.d0
      time_scheme = "dg"
      if (nargs .EQ. 8) call ReadStringArgument(8, time_scheme)
      if (nargs .EQ. 9) then
         call ReadRealArgument(8, scheme_tau)
         call ReadStringArgument(9, time_scheme)
      end if

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

      subroutine ReadRealArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         real(kind = 8), intent(out) :: value
         character(100) :: input
         integer(kind = 4) :: read_status

         call ReadArgument(arg, input)
         read(input, *, iostat = read_status) value
         if (read_status /= 0) then
            write(*,*) "invalid real argument: ", arg
            STOP 5
         end if
      end subroutine ReadRealArgument

      subroutine ReadStringArgument(arg, value)
         implicit none
         integer(kind = 4), intent(in) :: arg
         character(*), intent(out) :: value
         character(100) :: input

         call ReadArgument(arg, input)
         value = adjustl(input)
         call Lowercase(value)
      end subroutine ReadStringArgument

      subroutine Lowercase(value)
         implicit none
         character(*), intent(inout) :: value
         integer(kind = 4) :: i, code

         do i = 1, len_trim(value)
            code = iachar(value(i:i))
            if (code >= iachar('A') .AND. code <= iachar('Z')) value(i:i) = achar(code + 32)
         end do
      end subroutine Lowercase

   end subroutine InitializeParameters



end module input_data
