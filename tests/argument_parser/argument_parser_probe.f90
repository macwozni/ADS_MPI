program argument_parser_probe
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_value
   use argument_parser, only: ReadArgument, ReadIntegerArgument, ReadRealArgument, &
                              ReadStringArgument, ReadTimeSchemeArgument, &
                              RequireNonnegativeInteger, RequirePositiveInteger, RequirePositiveReal, &
                              RequireSafeSplineDimensions
   implicit none

   character(len=32) :: mode, expectation
   character(len=100) :: value, scheme
   integer(kind=4) :: degree_increment, integer_value, selected_scheme
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
      case ("leading-dot")
         if (abs(real_value - 0.5d0) > 1.0d-14) stop 93
      case ("trailing-dot")
         if (abs(real_value - 5.d0) > 1.0d-14) stop 93
      case ("integer-form")
         if (abs(real_value - 1.d0) > 1.0d-14) stop 93
      case ("e-exponent")
         if (abs(real_value - 1.d3) > 1.0d-11) stop 93
      case default
         stop 94
      end select
      write(*, '(A)') "real:ok"
   case ("positive-integer")
      call ReadIntegerArgument(2, integer_value)
      call RequirePositiveInteger(integer_value, "probe integer")
      write(*, '(A)') "positive-integer:ok"
   case ("nonnegative-integer")
      call ReadIntegerArgument(2, integer_value)
      call RequireNonnegativeInteger(integer_value, "probe integer")
      write(*, '(A)') "nonnegative-integer:ok"
   case ("positive-real")
      call ReadRealArgument(2, real_value)
      call RequirePositiveReal(real_value, "probe real")
      write(*, '(A)') "positive-real:ok"
   case ("positive-real-nonfinite")
      real_value = ieee_value(0.0d0, ieee_positive_inf)
      call RequirePositiveReal(real_value, "probe real")
      write(*, '(A)') "positive-real-nonfinite:unexpected-success"
   case ("safe-spline")
      call ReadIntegerArgument(2, integer_value)
      call ReadIntegerArgument(3, selected_scheme)
      if (command_argument_count() == 4) then
         call ReadIntegerArgument(4, degree_increment)
         call RequireSafeSplineDimensions(integer_value, selected_scheme, degree_increment)
      else
         call RequireSafeSplineDimensions(integer_value, selected_scheme)
      end if
      write(*, '(A)') "safe-spline:ok"
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
