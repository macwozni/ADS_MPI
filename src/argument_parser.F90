!------------------------------------------------------------------------------
!
! MODULE: argument_parser
!
! DESCRIPTION:
!> @file argument_parser.F90
!> @brief Shared command-line argument parsing helpers.
!>
!> @details
!> This module centralizes the small typed readers used by problem-specific
!> input_data modules. It keeps argument fetching, conversion checks, and
!> simple string normalization in one place.
!
!------------------------------------------------------------------------------
module argument_parser

   implicit none

   private

   public :: ReadArgument
   public :: ReadIntegerArgument
   public :: ReadRealArgument
   public :: ReadStringArgument
   public :: ReadTimeSchemeArgument
   public :: SelectTimeScheme
   public :: Lowercase
   public :: TIME_SCHEME_DG
   public :: TIME_SCHEME_PR
   public :: TIME_SCHEME_BE

!> @brief Supported iGRM time-scheme identifiers.
   integer(kind = 4), parameter :: TIME_SCHEME_DG = 1
   integer(kind = 4), parameter :: TIME_SCHEME_PR = 2
   integer(kind = 4), parameter :: TIME_SCHEME_BE = 3

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads one raw command-line argument.
!
! Input:
! ------
!> @param[in] arg
!> One-based command-line argument index.
!
! Output:
! -------
!> @param[out] input
!> Raw command-line argument value.
!
!---------------------------------------------------------------------------
subroutine ReadArgument(arg, input)
      implicit none
!> @brief One-based command-line argument index.
      integer(kind = 4), intent(in) :: arg
!> @brief Raw command-line argument value.
      character(*), intent(out) :: input
!> @brief Argument length and retrieval status.
      integer(kind = 4) :: length, status

      call GET_COMMAND_ARGUMENT(arg, input, length, status)
      if (status /= 0) then
         write(*,*) "invalid command argument: ", arg
         STOP 5
      end if

end subroutine ReadArgument

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads one integer command-line argument.
!
! Input:
! ------
!> @param[in] arg
!> One-based command-line argument index.
!
! Output:
! -------
!> @param[out] value
!> Parsed integer value.
!
!---------------------------------------------------------------------------
subroutine ReadIntegerArgument(arg, value)
      implicit none
!> @brief One-based command-line argument index.
      integer(kind = 4), intent(in) :: arg
!> @brief Parsed integer value.
      integer(kind = 4), intent(out) :: value
!> @brief Raw argument text.
      character(100) :: input
!> @brief Numeric conversion status.
      integer(kind = 4) :: read_status

      call ReadArgument(arg, input)
      read(input, *, iostat = read_status) value
      if (read_status /= 0) then
         write(*,*) "invalid integer argument: ", arg
         STOP 5
      end if

end subroutine ReadIntegerArgument

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads one real-valued command-line argument.
!
! Input:
! ------
!> @param[in] arg
!> One-based command-line argument index.
!
! Output:
! -------
!> @param[out] value
!> Parsed real value.
!
!---------------------------------------------------------------------------
subroutine ReadRealArgument(arg, value)
      implicit none
!> @brief One-based command-line argument index.
      integer(kind = 4), intent(in) :: arg
!> @brief Parsed real value.
      real(kind = 8), intent(out) :: value
!> @brief Raw argument text.
      character(100) :: input
!> @brief Numeric conversion status.
      integer(kind = 4) :: read_status

      call ReadArgument(arg, input)
      read(input, *, iostat = read_status) value
      if (read_status /= 0) then
         write(*,*) "invalid real argument: ", arg
         STOP 5
      end if

end subroutine ReadRealArgument

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads one string command-line argument and lowercases it.
!
! Input:
! ------
!> @param[in] arg
!> One-based command-line argument index.
!
! Output:
! -------
!> @param[out] value
!> Parsed, left-adjusted, lowercase string value.
!
!---------------------------------------------------------------------------
subroutine ReadStringArgument(arg, value)
      implicit none
!> @brief One-based command-line argument index.
      integer(kind = 4), intent(in) :: arg
!> @brief Parsed, normalized string value.
      character(*), intent(out) :: value
!> @brief Raw argument text.
      character(100) :: input

      call ReadArgument(arg, input)
      value = adjustl(input)
      call Lowercase(value)

end subroutine ReadStringArgument

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads and normalizes one iGRM time-scheme argument.
!
! Input:
! ------
!> @param[in] arg
!> One-based command-line argument index.
!
! Output:
! -------
!> @param[out] selected_time_scheme
!> Normalized time-scheme identifier.
!>
!> @param[out] time_scheme
!> Canonical short time-scheme name.
!
!---------------------------------------------------------------------------
subroutine ReadTimeSchemeArgument(arg, selected_time_scheme, time_scheme)
      implicit none
!> @brief One-based command-line argument index.
      integer(kind = 4), intent(in) :: arg
!> @brief Normalized time-scheme identifier.
      integer(kind = 4), intent(out) :: selected_time_scheme
!> @brief Canonical short time-scheme name.
      character(*), intent(out) :: time_scheme
!> @brief Raw argument text.
      character(100) :: input

      call ReadArgument(arg, input)
      call SelectTimeScheme(input, selected_time_scheme, time_scheme)

end subroutine ReadTimeSchemeArgument

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Maps a time-scheme name or alias to its canonical identifier.
!
! Input:
! ------
!> @param[in] value
!> Time-scheme name or alias.
!
! Output:
! -------
!> @param[out] selected_time_scheme
!> Normalized time-scheme identifier.
!>
!> @param[out] time_scheme
!> Canonical short time-scheme name.
!
!---------------------------------------------------------------------------
subroutine SelectTimeScheme(value, selected_time_scheme, time_scheme)
      implicit none
!> @brief Time-scheme name or alias.
      character(*), intent(in) :: value
!> @brief Normalized time-scheme identifier.
      integer(kind = 4), intent(out) :: selected_time_scheme
!> @brief Canonical short time-scheme name.
      character(*), intent(out) :: time_scheme
!> @brief Lowercase working copy.
      character(len = 32) :: normalized

      normalized = adjustl(value)
      call Lowercase(normalized)
      select case (trim(normalized))
      case ("dg", "douglas-gunn", "douglas_gunn")
         selected_time_scheme = TIME_SCHEME_DG
         time_scheme = "dg"
      case ("pr", "peaceman-rachford", "peaceman_rachford")
         selected_time_scheme = TIME_SCHEME_PR
         time_scheme = "pr"
      case ("be", "backward-euler", "backward_euler", "backwardeuler")
         selected_time_scheme = TIME_SCHEME_BE
         time_scheme = "be"
      case ("fe", "forward-euler", "forward_euler", "forwardeuler")
         write(*, *) "forward euler is not an iGRM time-scheme option"
         STOP 5
      case default
         write(*, *) "unknown time scheme: ", trim(value)
         STOP 5
      end select

end subroutine SelectTimeScheme

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Converts ASCII letters in a string to lowercase in place.
!
! Input/Output:
! -------------
!> @param[inout] value
!> String normalized in place.
!
!---------------------------------------------------------------------------
subroutine Lowercase(value)
      implicit none
!> @brief String normalized in place.
      character(*), intent(inout) :: value
!> @brief Character index and ASCII code.
      integer(kind = 4) :: i, code

      do i = 1, len_trim(value)
         code = iachar(value(i:i))
         if (code >= iachar('A') .AND. code <= iachar('Z')) value(i:i) = achar(code + 32)
      end do

end subroutine Lowercase

end module argument_parser
