program argument_parser_probe
   use argument_parser, only: ReadArgument, ReadIntegerArgument, ReadRealArgument, &
                              ReadStringArgument, ReadTimeSchemeArgument
   implicit none

   character(len=32) :: mode, expectation
   character(len=100) :: value, scheme
   integer(kind=4) :: integer_value, selected_scheme
   real(kind=8) :: real_value

   if (command_argument_count() < 1) stop 90
   call get_command_argument(1, mode)

   select case (trim(mode))
   case ("raw")
      call ReadArgument(2, value)
      write(*, '(A)') "raw:" // trim(value)
   case ("integer")
      call ReadIntegerArgument(2, integer_value)
      write(*, '("integer:",I0)') integer_value
   case ("real")
      call ReadRealArgument(2, real_value)
      if (command_argument_count() < 3) stop 92
      call get_command_argument(3, expectation)
      select case (trim(expectation))
      case ("negative-decimal")
         if (abs(real_value + 12.5d0) > 1.0d-14) stop 93
      case ("d-exponent")
         if (abs(real_value - 6.25d-3) > 1.0d-14) stop 93
      case default
         stop 94
      end select
      write(*, '(A)') "real:ok"
   case ("string")
      call ReadStringArgument(2, value)
      write(*, '(A)') "string:" // trim(value)
   case ("scheme")
      call ReadTimeSchemeArgument(2, selected_scheme, scheme)
      write(*, '("scheme:",I0,":",A)') selected_scheme, trim(scheme)
   case default
      stop 91
   end select

end program argument_parser_probe
